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
	// std::cout << "testing lobatto:" << std::endl;
	// for (int i=0; i<10; i++)
	// {
	// 	std::cout << "lobatto(" << i << ", 0.0) = " << lobatto(i, 0.0) << std::endl;
	// }

	const int n = 4; // number of elements
	const int polynomialDegree = 2; // polynomial degree
	const int dimension = 1; // dimension
	auto f = [](double x) { return M_PI * M_PI * sin(M_PI * x); }; // RHS function
	auto f_analytic = [](double x) { return sin(M_PI * x); }; // analytic solution

	std::string filenameNoExt = "domain.3";

	// boundary conditions
	double u0 = 0.0;
	double u1 = 0.0;
	bool boundary_u0 = true;
	bool boundary_u1 = true;

	FE_Solution FEM(n, polynomialDegree, dimension);
	if (dimension == 1)
	{
		std::vector<double> solution = FEM.solve(f);
	}
	else if (dimension == 2)
	{
		std::vector<double> solution = FEM.solve(f, filenameNoExt);
	}

	std::cout << "u vector: " << std::endl;
	for (int i=0; i<polynomialDegree*n+1; i++)
	{
		std::cout << "u" << i << " = " << solution.at(i) << std::endl;
	}

	std::cout << "Solution:" << std::endl;
	for (int i=0; i<5; i++)
	{
		double x = i / 4.0;
		std::cout << "u_h(" << x << ") = " << FEM.evaluateSolution(x) << std::endl;
	}
	// should be [0, 0.09375, 0.125, 0.09375, 0]

	std::ofstream file("hp_error.csv");
	file << "p,h,error\n";

	for (int p=1; p<=3; p++)
	{
		for (int n=3; n<=20; n++)
		{
			FE_Solution FEM(n, p, dimension);
			std::vector<double> solution = FEM.solve(f);
			double error = FEM.getL2Error(f_analytic);
			double h = 1.0 / n;
			file << p << "," << h << "," << error << "\n";
		}
	}

	file.close();

	return 0;
}