#ifndef FESOLUTIONHEADERDEF
#define FESOLUTIONHEADERDEF

#include <vector>
#include "CSR_Matrix.hpp"
#include "FE_Mesh.hpp"
#include "Solver.hpp"

class FE_Solution
{
	public:
		int n;
		CSRMatrix globalStiffness;
		std::vector<double> globalLoad;
		std::vector<double> solution;
		FE_Mesh mesh;

		FE_Solution(int n);
		FE_Mesh initialise();
		std::vector<double> solve(double (*f)(double), double u0=0.0, double u1=0.0, bool boundary_u0=true, bool boundary_u1=true);
};

#endif