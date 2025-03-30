#include "FE_Solution.hpp"
#include "Polynomial.hpp"
#include <iostream>
#include <string>

#include "Point2D.hpp"
#include "Element2D.hpp"
#include "PolynomialSpace.hpp"

//TODO:
// - polynomial spaces per element?
// - more general boundary conditions
// 		- just non-zero
// 		- inhomogeneous
// - non-linear problems
// 		- non-linear in general
// 		- Stoke's flow
// 		- Navier-Stokes
//  	- non-linear solver
// - consider using Eigen for the CSR matrix throughout
// - memoize reference stiffness
// - CSR_Matrix constructor

int main()
{
	// setup problem
	const int n = 4; // number of elements - UNUSED FOR >1D
	const int polynomialDegree = 1; // polynomial degree
	const int dimension = 2; // dimension

	std::vector<double> solution;
	std::string filenameNoExt = "domain.4";

	// boundary conditions
	double u0 = 0.0;
	double u1 = 0.0;
	bool boundary_u0 = true;
	bool boundary_u1 = true;

	if (dimension == 1)
	{
		// -u'' = f
		// f(x) = pi^2 sin(pi x)
		// u(x) = sin(pi x)
		auto f = [](const std::vector<double>& x) { return M_PI * M_PI * sin(M_PI * x[0]); }; // RHS function
		auto f_analytic = [](const std::vector<double>& x) { return sin(M_PI * x[0]); }; // analytic solution
		// create and solve problem
		FE_Solution FEM(n, polynomialDegree, dimension);
		solution = FEM.solve(f);

		// output u vector
		std::cout << "u vector: " << std::endl;
		for (int i=0; i<solution.size(); i++)
		{
			std::cout << "u" << i << " = " << solution[i] << std::endl;
		}
		
		// output solution
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
		// -grad^2 u = f
		// f(x, y) = 2pi^2 sin(pi x) sin(pi y)
		// u(x, y) = sin(pi x) sin(pi y)
		auto f2 = [](const std::vector<double>& x) { return 2.0 * M_PI * M_PI * sin(M_PI * x[0]) * sin(M_PI * x[1]); }; // RHS function
		auto f2_analytic = [](const std::vector<double>& x) { return sin(M_PI * x[0]) * sin(M_PI * x[1]); }; // analytic solution
		// create and solve problem
		FE_Solution FEM(n, polynomialDegree, dimension);
		solution = FEM.solve(f2, filenameNoExt);

		std::cout << "u vector: " << std::endl;
		for (int i=0; i<solution.size(); i++)
		{
			std::cout << "u" << i << " = " << solution[i] << std::endl;
		}
		
		int nGridPoints = 17;  // number of grid points to evaluate solution at
		// neat solution output
		// std::cout << "Solution:" << std::endl;
		// for (int i=0; i<nGridPoints; i++)
		// {
		// 	for (int j=0; j<nGridPoints; j++)
		// 	{
		// 		double x = (double)i / ((double)nGridPoints - 1.0);
		// 		double y = (double)j / ((double)nGridPoints - 1.0);
		// 		std::cout << "u(" << x << ", " << y << ") = " << FEM.evaluateSolution({x, y}) << std::endl;
		// 	}
		// }
		
		// solution for plotting
		std::cout << "Solution for plotting:" << std::endl;
		for (int i=0; i<nGridPoints; i++)
		{
			for (int j=0; j<nGridPoints; j++)
			{
				double x = (double)i / ((double)nGridPoints - 1.0);
				double y = (double)j / ((double)nGridPoints - 1.0);
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

	// // Define a constant function f(x,y)=1
	// auto f_constant = [](const std::vector<double>& x) -> double {
	// 	return 1.0;
	// };

	// // Define triangle vertices: (0,0), (1,0), (0,1)
	// Point2D P{0, 0.0, 0.0};
	// Point2D Q{1, 1.0, 0.0};
	// Point2D R{2, 0.0, 1.0};

	// std::vector<Point2D> triangle_nodes = {P, Q, R};

	// // Assume local_DoF indices for this element: {0, 1, 2}
	// std::vector<int> local_DoF = {0, 1, 2};

	// Element2D element(0, 1, local_DoF, triangle_nodes);

	// // Compute the load vector
	// std::vector<double> localLoad = element.getLocalLoad(f_constant);

	// std::cout << "Computed load vector:" << std::endl;
	// for (double v : localLoad) {
	// 	std::cout << v << std::endl;
	// }

	return 0;
}