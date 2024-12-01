#include "FE_Solution.hpp"
#include <iostream>

//TODO:
// - higher order polynomial spaces per element
// - proper wuadrature implementation
// - more general boundary conditions
// - fix boundary conditions
// - 2D
// - non-linear problems
// - consider using Eigen for the CSR matrix throughout
// - 

int main()
{
	const int n = 4; // number of elements
	auto f = [](double x) { return 1.0; }; // RHS function

	// boundary conditions
	double u0 = 0.0;
	double u1 = 0.0;
	bool boundary_u0 = true;
	bool boundary_u1 = true;

	std::vector<double> solution; // empty solution vector

	FE_Solution FEM(n);
	FEM.initialise();
	solution = FEM.solve(f, u0, u1, boundary_u0, boundary_u1);

	std::cout << "Solution: " << std::endl;
	for (int i=0; i<n+1; i++)
	{
		std::cout << "u(" << i * (1.0/n) << ") = " << solution[i] << std::endl;
	}

	return 0;
}