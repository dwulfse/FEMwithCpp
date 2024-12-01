#include "FE_Mesh.hpp"
#include <iostream>

FE_Mesh::FE_Mesh(int n)
 : n(n), nodes(n+1), elements(n), stiffness(), load(n+1)
{
}

FE_Mesh::FE_Mesh() {} // maybe initialise all to 0

void FE_Mesh::constructMesh()
{
	const double h = 1.0/n;

	nodes[0] = 0.0;
	for (int i=1; i<n+1; i++)
	{
		nodes[i] = i * h;
		Element elem({nodes[i-1], nodes[i]}, {i-1, i});
		elements[i-1] = elem;
	}
}

void FE_Mesh::allocateStiffness()
{
	std::vector<int> nnz_per_row(n+1);
	std::vector<int> DoF(elements[0].node_indicies.size());

	for (int i=0; i<n; i++) // for each element
	{
		DoF = elements[i].node_indicies;
		for (int j=0; j<DoF.size(); j++) // for each connection
		{
			nnz_per_row[DoF[j]] += DoF.size();
		}
	}

	stiffness.row_start.resize(n+2);
	for (int i=1; i<n+2; i++) // for each node (row)
	{
		stiffness.row_start[i] = stiffness.row_start[i-1] + nnz_per_row[i-1];
	}

	int nnz = stiffness.row_start[n+1];
	stiffness.entries.resize(nnz);
	stiffness.col_no.resize(nnz);

	std::vector<int> row_start_copy = stiffness.row_start;

	for (int k=0; k<n; k++)
	{
		DoF = elements[k].node_indicies;
		for (int i=0; i<DoF.size(); i++)
		{
			for (int j=0; j<DoF.size(); j++)
			{
				stiffness.col_no[row_start_copy[DoF[i]]] = DoF[j];
				row_start_copy[DoF[i]]++;
			}
		}
	}
}

CSRMatrix FE_Mesh::assembleStiffnessMatrix()
{
	std::vector<std::vector<double>> local_stiffness(2, std::vector<double>(2));
	std::vector<int> globalDoF;

	for (int k=0; k<n; k++)
	{
		local_stiffness = elements[k].getLocalStiffness();
		std::cout << "Local Stiffness Matrix for element " << k << ":" << std::endl;
		for (int i=0; i<2; i++)
		{
			for (int j=0; j<2; j++)
			{
				std::cout << local_stiffness[i][j] << " ";
			}
			std::cout << std::endl;
		}
		globalDoF = elements[k].node_indicies;
		// handle DoF
		for (int i=0; i<globalDoF.size(); i++)
		{
			for (int j=0; j<globalDoF.size(); j++)
			{
				stiffness(globalDoF[i], globalDoF[j]) += local_stiffness[i][j];
			}
		}
	}

	std::cout << "Stiffness Matrix before BCs:" << std::endl;
	for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
					try{
						std::cout << stiffness(i, j) << "		";
					}
					catch (std::out_of_range& e)
					{
						std::cout << "0		";
					}
			}
			std::cout << std::endl;
	}

	std::cout << "entries: ";
	for (int i=0; i<stiffness.entries.size(); i++)
	{
		std::cout << stiffness.entries[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "col_no: ";
	for (int i=0; i<stiffness.col_no.size(); i++)
	{
		std::cout << stiffness.col_no[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "row_start: ";
	for (int i=0; i<stiffness.row_start.size(); i++)
	{
		std::cout << stiffness.row_start[i] << " ";
	}
	std::cout << std::endl;

	return stiffness;
}

std::vector<double> FE_Mesh::assembleLoadVector(double (*f)(double))
{
	std::vector<double> local_load(2);
	std::vector<int> globalDoF;

	for (int k=0; k<n; k++)
	{
		local_load = elements[k].getLocalLoad(f);
		globalDoF = elements[k].node_indicies;
		// handle DoF
		for (int i=0; i<globalDoF.size(); i++)
		{
			load[globalDoF[i]] += local_load[i];
		}
	}

	return load;
}

void FE_Mesh::applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1)
{
	if (boundary_u0)
	{
		for (int i=stiffness.row_start[0]; i<stiffness.row_start[1]; i++)
		{
			stiffness.entries[i] = 0.0;
		}
		stiffness(0, 0) = 1.0;
		load[0] = u0;
	}

	if (boundary_u1)
	{
		for (int i=stiffness.row_start[n]; i<stiffness.row_start[n+1]; i++)
		{
			stiffness.entries[i] = 0.0;
		}
		stiffness(n, n) = 1.0;
		load[n] = u1;
	}

	std::cout << "Stiffness Matrix w/ BCs:" << std::endl;
	for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
					try {
						std::cout << stiffness(i, j) << "		";
					}
					catch (std::out_of_range& e)
					{
						std::cout << "0		";
					}
			}
			std::cout << std::endl;
	}

	std::cout << "entries: ";
	for (int i=0; i<stiffness.entries.size(); i++)
	{
		std::cout << stiffness.entries[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "col_no: ";
	for (int i=0; i<stiffness.col_no.size(); i++)
	{
		std::cout << stiffness.col_no[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "row_start: ";
	for (int i=0; i<stiffness.row_start.size(); i++)
	{
		std::cout << stiffness.row_start[i] << " ";
	}
	std::cout << std::endl;
}