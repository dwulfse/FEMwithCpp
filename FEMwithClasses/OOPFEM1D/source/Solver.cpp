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
	EigenA.resize(p*n+1, p*n+1);
	Eigenb.resize(p*n+1);

	std::vector<Eigen::Triplet<double>> triplets;
	for (int i=0; i<p*n+1; i++) ///////////
	{
		for (int j=A.row_start.at(i); j<A.row_start[i+1]; j++)
		{
			triplets.emplace_back(Eigen::Triplet<double>(i, A.col_no.at(j), A.entries.at(j)));
		}
	}

	// log triplets
	// for (int i=0; i<triplets.size(); i++)
	// {
	// 	std::cout << "triplets[" << i << "] = (" << triplets.at(i).row() << ", " << triplets.at(i).col() << ", " << triplets.at(i).value() << ")" << std::endl;
	// }

	EigenA.setFromTriplets(triplets.begin(), triplets.end()); // problem for now
	EigenA.makeCompressed();

	for (int i=0; i<p*n+1; i++)
	{
		Eigenb(i) = b.at(i);
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
	for (int i=0; i<p*n+1; i++)
	{
		u.at(i) = solution(i);
	}

	return u;
}