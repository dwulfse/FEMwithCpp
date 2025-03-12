#include "FE_Mesh.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

// TODO:
// - account for extra DoFs in higher order elements

// algorithm to allocate space for CSR stiffness matrix based on where
// non-zero entries will be
void FE_Mesh::allocateStiffness()
{
	int no_nodes = getNoNodes();
	std::vector<int> nnz_per_row(no_nodes, 0); // not so sure about p*n + 1 but its excess
	// std::vector<int> DoF(p+1); // hp-fem - change !!

	for (int k=0; k<elements.size(); k++) // for each element
	{
		const std::vector<int>& DoF = elements[k]->local_DoF;
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

	for (int k=0; k<elements.size(); k++)
	{
		std::vector<int>& DoF = elements[k]->local_DoF;
		for (int i=0; i<DoF.size(); i++)
		{
			for (int j=0; j<DoF.size(); j++)
			{
				int row = DoF[i];
				int pos = row_start_copy[row];
				if (pos >= stiffness.row_start[row+1])
				{
					std::cerr << "Error: row " << row << " insertion index " << pos
							  << " exceeds allocated range ["
							  << stiffness.row_start[row] << ", "
							  << stiffness.row_start[row+1] << ")\n";
					throw std::runtime_error("CSR assembly index out of range");
				}
				stiffness.col_no[pos] = DoF[j];
				row_start_copy[row]++;
			}
		}
	}

	std::cout << "Stiffness matrix: " << std::endl;
	stiffness.print(true);
}

// find each element's local stiffness and construct global stiffness using each
// element's DoFs
void FE_Mesh::assembleStiffnessMatrix()
{
	int noLocalDoFs = d == 1 ? p+1 : (p+1)*(p+2)/2;
	std::vector<std::vector<double>> local_stiffness(noLocalDoFs, std::vector<double>(noLocalDoFs, 0)); // changes needed for hp-fem
	std::vector<int> globalDoF(noLocalDoFs);

	for (int k=0; k<n; k++)
	{
		local_stiffness = elements[k]->getLocalStiffness();
		globalDoF = elements[k]->local_DoF;
		// handle DoF
		for (int i=0; i<noLocalDoFs; i++)
		{
			for (int j=0; j<noLocalDoFs; j++)
			{
				stiffness(globalDoF[i], globalDoF[j]) += local_stiffness[i][j];
			}
		}
	}

	std::cout << "Stiffness matrix: " << std::endl;
	stiffness.print(false);
}

// find each element's local load and construct global load using each
// element's DoFs
void FE_Mesh::assembleLoadVector(double (*f)(const std::vector<double>&))
{
	int noLocalDoFs = d == 1 ? p+1 : (p+1)*(p+2)/2;
	std::vector<double> local_load(noLocalDoFs); // changes needed for hp-fem
	std::vector<int> globalDoF(noLocalDoFs);

	for (int k=0; k<n; k++)
	{
		local_load = elements[k]->getLocalLoad(f);
		globalDoF = elements[k]->local_DoF;
		// handle DoF
		for (int i=0; i<noLocalDoFs; i++)
		{
			load[globalDoF[i]] += local_load[i];
		}
	}

	std::cout << "Load vector: " << std::endl;
	for (int i=0; i<load.size(); i++)
	{
		std::cout << load[i] << std::endl;
	}
}