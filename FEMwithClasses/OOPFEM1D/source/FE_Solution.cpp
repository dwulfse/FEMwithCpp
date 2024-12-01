#include "FE_Solution.hpp"

FE_Solution::FE_Solution(int n)
 : n(n), globalStiffness(), globalLoad(n+1), solution(n+1), mesh(n)
{
}

FE_Mesh FE_Solution::initialise()
{
	mesh.constructMesh();

	return mesh;
}

std::vector<double> FE_Solution::solve(double (*f)(double), double u0, double u1, bool boundary_u0, bool boundary_u1)
{
	mesh.allocateStiffness();
	globalStiffness = mesh.assembleStiffnessMatrix();
	globalLoad = mesh.assembleLoadVector(f);
	mesh.applyBoundaryConditions(u0, u1, boundary_u0, boundary_u1);

	Solver solver(n, globalStiffness, globalLoad, solution);
	// solution = solver.solveSystem();

	solver.setupEigen();
	solution = solver.solveEigen();

	return solution;
}