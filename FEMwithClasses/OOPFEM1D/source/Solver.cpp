#include "Solver.hpp"

Solver::Solver(int n, CSRMatrix A, std::vector<double> b, std::vector<double> u)
 : n(n), A(A), b(b), u(u)
{
}

void Solver::setupEigen()
{
	EigenA.resize(n+1, n+1);
	Eigenb.resize(n+1);

	std::vector<Eigen::Triplet<double>> triplets;
	for (int i=0; i<n+1; i++)
	{
		for (int j=A.row_start[i]; j<A.row_start[i+1]; j++)
		{
			triplets.emplace_back(Eigen::Triplet<double>(i, A.col_no[j], A.entries[j]));
		}
	}

	EigenA.setFromTriplets(triplets.begin(), triplets.end()); // problem for now
	EigenA.makeCompressed();

	for (int i=0; i<n+1; i++)
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
	for (int i=0; i<n+1; i++)
	{
		u[i] = solution(i);
	}

	return u;
}