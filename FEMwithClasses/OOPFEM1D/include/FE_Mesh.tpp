#include "FE_Mesh.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

// TODO:
// - account for extra DoFs in higher order elements

// constructors
template <typename NodeType>
FE_Mesh<NodeType>::FE_Mesh() // maybe initialise all to 0
{
}

template <typename NodeType>
FE_Mesh<NodeType>::FE_Mesh(int n, int p)
 : n(n), p(p), nodes(p*n+1, 0.0), elements(n), stiffness(), load(p*n+1, 0.0)
{
}

// destructor
template <typename NodeType>
FE_Mesh<NodeType>::~FE_Mesh()
{
}

// load mesh from Triangle generated files
template <typename NodeType>
void FE_Mesh<NodeType>::parseMesh(std::string filename)
{
	std::ifstream eleFile(filename + ".ele");
	std::ifstream nodeFile(filename + ".node");

	if (!eleFile.good() || !nodeFile.good())
	{
		throw std::runtime_error("Mesh files not found");
	}

	std::string line;

	// read node header
	std::getline(nodeFile, line);
	std::istringstream nodeHeader(line);
	int noVertices, dimension, noAttributesNode, boundary;
	nodeHeader >> noVertices >> dimension >> noAttributesNode >> boundary;

	nodes.resize(noVertices);

	// read nodes
	while (std::getline(nodeFile, line))
	{
		std::istringstream iss(line);
		Point2D node;
		iss >> node.id >> node.x >> node.y;
		node.id--; // 0-indexed
		nodes[node.id] = node;
	}

	// read element header
	std::getline(eleFile, line);
	std::istringstream eleHeader(line);
	int noElems, noNodesElem, noAttributes;
	eleHeader >> noElems >> noNodesElem >> noAttributes;

	elements.resize(noElems);

	// read elements
	while(std::getline(eleFile, line))
	{
		std::istringstream iss(line);
		Element elem;
		elem.p = p;
		iss >> elem.id;
		elem.id--; // 0-indexed
		elem.local_DoF.resize(noNodesElem);
		for (int i=0; i<noNodesElem; i++)
		{
			iss >> elem.local_DoF[i];
			elem.local_DoF[i]--; // 0-indexed
		}
		elements[elem.id] = elem;
	}
}

// construct mesh of n elements and n+1 nodes
template <typename NodeType>
void FE_Mesh<NodeType>::constructMesh()
{
	const double h = 1.0 / n; // element size
	int int_i = n + 1; // start of extra nodes indexing (find proper word for this)

	// boundary nodes
	nodes.resize(n+1); // was p*n+1, haven't tested n+1 yet
	nodes.at(0) = 0.0;
	for (int i=1; i<n+1; i++)
	{
		nodes.at(i) = i * h;
	}

	// elements
	std::vector<int> local_dofs(p+1);
	std::vector<double> element_nodes(p+1);
	for (int i=0; i<n; i++)
	{
		// element dofs
		local_dofs.at(0) = i;
		local_dofs.at(1) = i+1;
		if (p > 1)
		{
			for (int j=2; j<p+1; j++)
			{
				local_dofs.at(j) = int_i++;
			}
		}
		// local_dofs.at(p) = i+1;

		elements.at(i) = Element({nodes.at(i), nodes.at(i+1)}, local_dofs, p);
	}
}

// algorithm to allocate space for CSR stiffness matrix based on where
// non-zero entries will be
template <typename NodeType>
void FE_Mesh<NodeType>::allocateStiffness()
{
	int no_nodes = p*n + 1;
	std::vector<int> nnz_per_row(no_nodes, 0); // not so sure about p*n + 1 but its excess
	// std::vector<int> DoF(p+1); // hp-fem - change !!

	for (int k=0; k<n; k++) // for each element
	{
		const std::vector<int>& DoF = elements.at(k).local_DoF;
		for (int j=0; j<DoF.size(); j++) // for each connection
		{
			nnz_per_row[DoF[j]] += DoF.size();
		}
	}

	stiffness.row_start.resize(no_nodes + 1, 0);
	for (int i=1; i<no_nodes + 1; i++) // for each node (row)
	{
		stiffness.row_start[i] = stiffness.row_start[i-1] + nnz_per_row[i-1];
	}

	int nnz = stiffness.row_start.back();
	stiffness.entries.resize(nnz, 0.0);
	stiffness.col_no.resize(nnz, 0);

	std::vector<int> row_start_copy = stiffness.row_start; // should use reference? maybe?

	for (int k=0; k<n; k++)
	{
		std::vector<int>& DoF = elements.at(k).local_DoF;
		for (int i=0; i<DoF.size(); i++)
		{
			for (int j=0; j<DoF.size(); j++)
			{
				// std::cout << "i=" << i << ", j=" << j << "\nDoF.at(i)=" << DoF.at(i) << ", DoF.at(j)=" << DoF.at(j) << "\nrow_start_copy.at(DoF.at(i))=" << row_start_copy[DoF.at(i)] << std::endl;
				stiffness.col_no[row_start_copy[DoF[i]]++] = DoF[j];
			}
		}
	}
}

// find each element's local stiffness and construct global stiffness using each
// element's DoFs
template <typename NodeType>
void FE_Mesh<NodeType>::assembleStiffnessMatrix()
{
	std::vector<std::vector<double>> local_stiffness(p+1, std::vector<double>(p+1, 0)); // changes needed for hp-fem
	std::vector<int> globalDoF(p+1);

	for (int k=0; k<n; k++)
	{
		local_stiffness = elements.at(k).getLocalStiffness();
		globalDoF = elements.at(k).local_DoF;
		// handle DoF
		for (int i=0; i<p+1; i++)
		{
			for (int j=0; j<p+1; j++)
			{
				stiffness(globalDoF.at(i), globalDoF.at(j)) += local_stiffness.at(i).at(j);
			}
		}
	}
}

// find each element's local load and construct global load using each
// element's DoFs
template <typename NodeType>
void FE_Mesh<NodeType>::assembleLoadVector(double (*f)(double))
{
	std::vector<double> local_load(p+1); // changes needed for hp-fem
	std::vector<int> globalDoF(p+1);

	for (int k=0; k<n; k++)
	{
		local_load = elements.at(k).getLocalLoad(f);
		globalDoF = elements.at(k).local_DoF;
		// handle DoF
		for (int i=0; i<p+1; i++)
		{
			load.at(globalDoF.at(i)) += local_load.at(i);
		}
	}
}

// apply left and right boundary conditions
template <typename NodeType>
void FE_Mesh<NodeType>::applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1)
{
	if (boundary_u0)
	{
		for (int i=stiffness.row_start.at(0); i<stiffness.row_start.at(1); i++)
		{
			stiffness.entries.at(i) = 0.0;
		}

		for (int i=1; i<stiffness.col_no.size(); i++) // probably inefficient, revisit
		{
			if (stiffness.col_no.at(i) == 0)
			{
				stiffness.entries.at(i) = 0.0;
			}
		}

		stiffness(0, 0) = 1.0;
		load.at(0) = u0;
	}

	if (boundary_u1)
	{
		for (int i=stiffness.row_start.at(n); i<stiffness.row_start.at(n+1); i++)
		{
			stiffness.entries.at(i) = 0.0;
		}

		for (int i=0; i<stiffness.col_no.size(); i++) // probably inefficient, revisit
		{
			if (stiffness.col_no.at(i) == n)
			{
				stiffness.entries.at(i) = 0.0;
			}
		}

		stiffness(n, n) = 1.0;
		load.at(n) = u1;
	}
}