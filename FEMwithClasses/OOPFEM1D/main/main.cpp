#include "FE_Solution.hpp"
#include "Polynomial.hpp"
#include <iostream>
#include <string>

//TODO:
// - polynomial spaces per element
// - more general boundary conditions
// - 2D
// - non-linear problems
// - consider using Eigen for the CSR matrix throughout
// - memoize reference stiffness
// - CSR_Matrix constructor

int main()
{
	const int n = 4; // number of elements
	const int polynomialDegree = 1; // polynomial degree
	const int dimension = 2; // dimension

	std::vector<double> solution;
	std::string filenameNoExt = "domain.1";

	// boundary conditions
	double u0 = 0.0;
	double u1 = 0.0;
	bool boundary_u0 = true;
	bool boundary_u1 = true;

	FE_Solution FEM(n, polynomialDegree, dimension);
	if (dimension == 1)
	{
		auto f = [](const std::vector<double>& x) { return M_PI * M_PI * sin(M_PI * x[0]); }; // RHS function
		auto f_analytic = [](const std::vector<double>& x) { return sin(M_PI * x[0]); }; // analytic solution
		solution = FEM.solve(f);

		std::cout << "u vector: " << std::endl;
		for (int i=0; i<solution.size(); i++)
		{
			std::cout << "u" << i << " = " << solution[i] << std::endl;
		}
	
		std::cout << "Solution:" << std::endl;
		for (int i=0; i<5; i++)
		{
			double x = i / 4.0;
			std::cout << "u_h(" << x << ") = " << FEM.evaluateSolution({x}) << std::endl;
		}
		// should be [0, 0.09375, 0.125, 0.09375, 0]
	}
	else if (dimension == 2)
	{
		auto f2 = [](const std::vector<double>& x) { return 2.0 * M_PI * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]); }; // RHS function
		auto f2_analytic = [](const std::vector<double>& x) { return sin(M_PI * x[0]) * sin(M_PI * x[1]); }; // analytic solution
		solution = FEM.solve(f2, filenameNoExt);

		std::cout << "u vector: " << std::endl;
		for (int i=0; i<solution.size(); i++)
		{
			std::cout << "u" << i << " = " << solution[i] << std::endl;
		}
	
		std::cout << "Solution:" << std::endl;
		int nGridPoints = 10;
		for (int i=0; i<nGridPoints+1; i++)
		{
			for (int j=0; j<nGridPoints+1; j++)
			{
				double x = i / (double)nGridPoints;
				double y = j / (double)nGridPoints;
				std::cout << x << ", " << y << ", " << FEM.evaluateSolution({x, y}) << std::endl;
			}
		}
	}

	// convergence analysis
	// std::ofstream file("hp_error.csv");
	// file << "p,h,error\n";

	// for (int p=1; p<=3; p++)
	// {
	// 	for (int n=3; n<=20; n++)
	// 	{
	// 		FE_Solution FEM(n, p, dimension);
	// 		std::vector<double> solution = FEM.solve(f);
	// 		double error = FEM.getL2Error(f_analytic);
	// 		double h = 1.0 / n;
	// 		file << p << "," << h << "," << error << "\n";
	// 	}
	// }

	// file.close();

	return 0;
}