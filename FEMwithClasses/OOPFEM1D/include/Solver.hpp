#ifndef SOLVERHEADERDEF
#define SOLVERHEADERDEF

#include "CSR_Matrix.hpp"
#include "../thirdparty/Eigen/Sparse"
#include "../thirdparty/Eigen/SparseLU"
#include <vector>

class Solver
{
	public:
		int n;
		CSRMatrix A;
		std::vector<double> b;
		std::vector<double> u;
		Eigen::SparseMatrix<double> EigenA;
		Eigen::VectorXd Eigenb;

		Solver(int n, CSRMatrix A, std::vector<double> b, std::vector<double> u);

		void setupEigen();
		std::vector<double> solveEigen();
};

#endif