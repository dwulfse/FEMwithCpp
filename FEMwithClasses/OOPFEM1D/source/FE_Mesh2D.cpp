#include "FE_Mesh2D.hpp"
#include "Element2D.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

// TODO:
// - account for extra DoFs in higher order elements

// constructors
FE_Mesh2D::FE_Mesh2D() // maybe initialise all to 0
{
}

FE_Mesh2D::FE_Mesh2D(int n, int p, int d)
 : FE_Mesh(n, p, d), nodes(p*n+1, {0, 0.0, 0.0}), boundary_nodes(0)
{
}

// destructor
FE_Mesh2D::~FE_Mesh2D()
{
}

// get node i coordinates
std::vector<double> FE_Mesh2D::getNode(int i)
{
	return { nodes[i].x, nodes[i].y };
}

// get number of nodes
int FE_Mesh2D::getNoNodes()
{
	return nodes.size();	// TODO: update when p>1 implemented
}

// load mesh from Triangle generated files
void FE_Mesh2D::constructMesh(std::string filename)
{
	std::ifstream eleFile("C:/Users/dylan/OneDrive/Documents/uni_work/FEMwithCPP/Code/FEMwithClasses/OOPFEM1D/main/domain.1.ele");
	std::ifstream nodeFile("C:/Users/dylan/OneDrive/Documents/uni_work/FEMwithCPP/Code/FEMwithClasses/OOPFEM1D/main/domain.1.node");

	if (!eleFile.good() || !nodeFile.good())
	{
		throw std::runtime_error("Mesh files not found");
	}

	std::string line;

	// read node header
	std::getline(nodeFile, line);
	if (line.empty())
	{
		throw std::runtime_error("Empty node file header");
	}
	std::istringstream nodeHeader(line);
	int noVertices, dimension, noAttributesNode, boundary;
	nodeHeader >> noVertices >> dimension >> noAttributesNode >> boundary;

	nodes.resize(noVertices);
	load.resize(noVertices, 0.0);

	// read nodes
	while (std::getline(nodeFile, line))
	{
		if (line.empty() || line[0] == '#') continue;
		std::istringstream iss(line);
		Point2D node;
		int boundary_flag;
		if (!(iss >> node.id >> node.x >> node.y >> boundary_flag))
		{
			std::cerr << "Couldn't read node file line: " << line << std::endl;
			continue;
		}
		node.id--; // 0-indexed
		if (node.id < 0 || node.id >= noVertices)
		{
			std::cerr << "Node id " << node.id << " out of range" << std::endl;
			continue;
		}
		nodes[node.id] = node;
		if (boundary_flag == 1)
		{
			boundary_nodes.push_back(node.id);
		}
	}

	std::cout << "Nodes: " << std::endl;
	for (int i=0; i<noVertices; i++)
	{
		std::cout << nodes[i].id << " " << nodes[i].x << " " << nodes[i].y << std::endl;
	}

	std::cout << "boundary_nodes: " << std::endl;
	for (int i=0; i<boundary_nodes.size(); i++)
	{
		std::cout << boundary_nodes[i] << " ";
	}
	std::cout << std::endl;

	// read element header
	std::getline(eleFile, line);
	if (line.empty())
	{
		throw std::runtime_error("Empty element file header");
	}
	std::istringstream eleHeader(line);
	int noElems, noNodesElem, noAttributes;
	eleHeader >> noElems >> noNodesElem >> noAttributes;
	if (noNodesElem != 3)
	{
		throw std::runtime_error("Only triangles supported");
	}

	n = noElems;
	elements.resize(noElems);

	// read elements
	std::vector<int> local_dofs(noNodesElem);
	std::vector<Point2D> element_nodes(noNodesElem);
	int id;
	while(std::getline(eleFile, line))
	{
		if (line.empty() || line[0] == '#') continue;
		std::istringstream iss(line);
		if (!(iss >> id))
		{
			std::cerr << "Couldn't read element id from line" << line << std::endl;
			continue;
		}
		id--; // 0-indexed
		if (id < 0 || id >= noElems)
		{
			std::cerr << "Element id " << id << " out of range" << std::endl;
			continue;
		}
		for (int i=0; i<noNodesElem; i++)
		{
			if (!(iss >> local_dofs[i]))
			{
				std::cerr << "Couldn't read local dof for element " << id+1 << " from line" << line << std::endl;
				break;
			}
			local_dofs[i]--; // 0-indexed
			if (local_dofs[i] < 0 || local_dofs[i] >= noVertices)
			{
				std::cerr << "Local dof " << local_dofs[i] << " out of range" << std::endl;
				break;
			}
			element_nodes[i] = nodes[local_dofs[i]];
		}
		elements[id] = std::make_unique<Element2D>(id, p, local_dofs, element_nodes);
	}

	std::cout << "Elements: " << std::endl;
	for (int i=0; i<noElems; i++)
	{
		std::cout << "id: " << elements[i]->id << " \tDoF: ";
		for (int j=0; j<noNodesElem; j++)
		{
			std::cout << elements[i]->local_DoF[j] << " ";
		}
		std::cout << std::endl;
	}

	eleFile.close();
	nodeFile.close();
}

// evaluate solution at x
double FE_Mesh2D::evaluateSolution(std::vector<double> x, std::vector<double> solution)
{
	for (int i=0; i<elements.size(); i++)
	{
		std::vector<int> elem_dof = elements[i]->local_DoF;
		Point2D P = nodes[elem_dof[0]], Q = nodes[elem_dof[1]], R = nodes[elem_dof[2]];
		double det = (Q.x - P.x) * (R.y - P.y) - (R.x - P.x) * (Q.y - P.y);

		std::vector<double> xi(3, 0.0);
		xi[0] = ((Q.x - x[0]) * (R.y - x[1]) - (R.x - x[0]) * (Q.y - x[1])) / det;
		xi[1] = ((R.x - x[0]) * (P.y - x[1]) - (P.x - x[0]) * (R.y - x[1])) / det;
		xi[2] = 1.0 - xi[0] - xi[1];

		// check if x in triangle
		if (xi[0] >= -1e-12 && xi[1] >= -1e-12 && xi[2] >= -1e-12)
		{
			double u_val = 0.0;
			for (int j=0; j<3; j++)
			{
				u_val += solution[elem_dof[j]] * xi[j];
			}
			return u_val;
		}
	}
	throw std::runtime_error("Point not in any element on mesh");
}

// apply boundary conditions
void FE_Mesh2D::applyBoundaryConditions(double u_val, double /*unused*/, bool apply_boundary, bool /*unused*/)
{
	if (!apply_boundary) return;

	for (int k=0; k<boundary_nodes.size(); k++)
	{
		// zero out row
		int index = boundary_nodes[k];
		for (int i=stiffness.row_start[index]; i<stiffness.row_start[index+1]; i++)
		{
			stiffness.entries[i] = 0.0;
		}

		// zero out column
		for (int i=0; i<stiffness.row_start.size(); i++)
		{
			for (int j=stiffness.row_start[i]; j<stiffness.row_start[i+1]; j++)
			{
				if (stiffness.col_no[j] == index)
				{
					stiffness.entries[j] = 0.0;
				}
			}
		}

		// set diagonal to 1
		stiffness(index, index) = 1.0;
		// set load vector to boundary value
		load[index] = u_val;
	}
}

// // algorithm to allocate space for CSR stiffness matrix based on where
// // non-zero entries will be
// void FE_Mesh2D::allocateStiffness()
// {
// 	int no_nodes = p*n + 1;
// 	std::vector<int> nnz_per_row(no_nodes, 0); // not so sure about p*n + 1 but its excess
// 	// std::vector<int> DoF(p+1); // hp-fem - change !!

// 	for (int k=0; k<n; k++) // for each element
// 	{
// 		const std::vector<int>& DoF = elements[k].local_DoF;
// 		for (int j=0; j<DoF.size(); j++) // for each connection
// 		{
// 			nnz_per_row[DoF[j]] += DoF.size();
// 		}
// 	}

// 	stiffness.row_start.resize(no_nodes + 1, 0);
// 	for (int i=1; i<no_nodes + 1; i++) // for each node (row)
// 	{
// 		stiffness.row_start[i] = stiffness.row_start[i-1] + nnz_per_row[i-1];
// 	}

// 	int nnz = stiffness.row_start.back();
// 	stiffness.entries.resize(nnz, 0.0);
// 	stiffness.col_no.resize(nnz, 0);

// 	std::vector<int> row_start_copy = stiffness.row_start; // should use reference? maybe?

// 	for (int k=0; k<n; k++)
// 	{
// 		std::vector<int>& DoF = elements[k].local_DoF;
// 		for (int i=0; i<DoF.size(); i++)
// 		{
// 			for (int j=0; j<DoF.size(); j++)
// 			{
// 				// std::cout << "i=" << i << ", j=" << j << "/nDoF[i]=" << DoF[i] << ", DoF[j]=" << DoF[j] << "\nrow_start_copy[DoF[i]]=" << row_start_copy[DoF[i]] << std::endl;
// 				stiffness.col_no[row_start_copy[DoF[i]]++] = DoF[j];
// 			}
// 		}
// 	}
// }

// // find each element's local stiffness and construct global stiffness using each
// // element's DoFs
// void FE_Mesh2D::assembleStiffnessMatrix()
// {
// 	std::vector<std::vector<double>> local_stiffness(p+1, std::vector<double>(p+1, 0)); // changes needed for hp-fem
// 	std::vector<int> globalDoF(p+1);

// 	for (int k=0; k<n; k++)
// 	{
// 		local_stiffness = elements[k].getLocalStiffness();
// 		globalDoF = elements[k].local_DoF;
// 		// handle DoF
// 		for (int i=0; i<p+1; i++)
// 		{
// 			for (int j=0; j<p+1; j++)
// 			{
// 				stiffness(globalDoF[i], globalDoF[j]) += local_stiffness[i][j];
// 			}
// 		}
// 	}
// }

// // find each element's local load and construct global load using each
// // element's DoFs
// void FE_Mesh2D::assembleLoadVector(double (*f)(double))
// {
// 	std::vector<double> local_load(p+1); // changes needed for hp-fem
// 	std::vector<int> globalDoF(p+1);

// 	for (int k=0; k<n; k++)
// 	{
// 		local_load = elements[k].getLocalLoad(f);
// 		globalDoF = elements[k].local_DoF;
// 		// handle DoF
// 		for (int i=0; i<p+1; i++)
// 		{
// 			load[globalDoF.at(i]) += local_load[i];
// 		}
// 	}
// }

// // apply left and right boundary conditions
// void FE_Mesh2D::applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1)
// {
// 	if (boundary_u0)
// 	{
// 		for (int i=stiffness.row_start[0]; i<stiffness.row_start[1]; i++)
// 		{
// 			stiffness.entries[i] = 0.0;
// 		}

// 		for (int i=1; i<stiffness.col_no.size(); i++) // probably inefficient, revisit
// 		{
// 			if (stiffness.col_no[i] == 0)
// 			{
// 				stiffness.entries[i] = 0.0;
// 			}
// 		}

// 		stiffness(0, 0) = 1.0;
// 		load[0] = u0;
// 	}

// 	if (boundary_u1)
// 	{
// 		for (int i=stiffness.row_start[n]; i<stiffness.row_start[n+1]; i++)
// 		{
// 			stiffness.entries[i] = 0.0;
// 		}

// 		for (int i=0; i<stiffness.col_no.size(); i++) // probably inefficient, revisit
// 		{
// 			if (stiffness.col_no[i] == n)
// 			{
// 				stiffness.entries[i] = 0.0;
// 			}
// 		}

// 		stiffness(n, n) = 1.0;
// 		load[n] = u1;
// 	}
// }