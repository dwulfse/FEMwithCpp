#include "Solver.hpp"
#include <iostream>

// constructor
Solver::Solver(int n, int p, CSRMatrix A, std::vector<double> b, std::vector<double> u)
 : n(n), p(p), A(A), b(b), u(u)
{
}

// destructor
Solver::~Solver()
{
}

void Solver::setupEigen()
{
	int noNodes = b.size(); 	// was p*n + 1
	EigenA.resize(noNodes, noNodes);
	Eigenb.resize(noNodes);

	std::vector<Eigen::Triplet<double>> triplets;
	for (int i=0; i<noNodes; i++) ///////////
	{
		for (int j=A.row_start[i]; j<A.row_start[i+1]; j++)
		{
			triplets.emplace_back(Eigen::Triplet<double>(i, A.col_no[j], A.entries[j]));
		}
	}

	// log triplets
	// for (int i=0; i<triplets.size(); i++)
	// {
	// 	std::cout << "triplets[" << i << "] = (" << triplets[i].row() << ", " << triplets[i].col() << ", " << triplets[i].value() << ")" << std::endl;
	// }

	EigenA.setFromTriplets(triplets.begin(), triplets.end()); // problem for now
	EigenA.makeCompressed();

	std::cout << "Eigen compressed matrix: " << std::endl;
	std::cout << EigenA << std::endl;

	for (int i=0; i<noNodes; i++)
	{
		Eigenb(i) = b[i];
	}
}

std::vector<double> Solver::solveEigen()
{
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(EigenA);

	if (solver.info() != Eigen::Success)
	{
		throw std::runtime_error("Eigen factorisation failed");
	}

	Eigen::VectorXd solution = solver.solve(Eigenb);

	// convert Eigen::VectorXd to std::vector<double> u
	for (int i=0; i<b.size(); i++)
	{
		u[i] = solution(i);
	}

	return u;
}