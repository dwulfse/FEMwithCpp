#include "FE_Solution.hpp"
#include "FE_Mesh1D.hpp"
#include "FE_Mesh2D.hpp"

// constructor
FE_Solution::FE_Solution(int n, int p, int d)
 : n(n), globalStiffness(), globalLoad(p*n+1), solution(p*n+1), p(p), d(d)
{
	if (d == 1)
	{
		mesh = std::make_unique<FE_Mesh1D>(n, p);
	}
	else if (d == 2)
	{
		mesh = std::make_unique<FE_Mesh2D>(n, p);
	}
}

// destructor
FE_Solution::~FE_Solution()
{
}

// controller function for whole program
std::vector<double> FE_Solution::solve(double (*f)(double), std::string fileNameNoExt, double u0, double u1, bool boundary_u0, bool boundary_u1)
{
	mesh->constructMesh(d == 2 ? fileNameNoExt : "");
	mesh->allocateStiffness();
	mesh->assembleStiffnessMatrix();
	mesh->assembleLoadVector(f);
	mesh->applyBoundaryConditions(u0, u1, boundary_u0, boundary_u1);

	globalStiffness = mesh->stiffness;
	globalLoad = mesh->load;

	Solver solver(n, p, globalStiffness, globalLoad, solution);
	// solution = solver.solveSystem();

	solver.setupEigen();
	solution = solver.solveEigen();

	return solution;
}

double FE_Solution::evaluateSolution(std::vector<double> x)
{
	return mesh->evaluateSolution(x, solution);
}

double FE_Solution::getL2Error(double (*f_analytic)(double))
{
	double L2error = 0.0;
	GaussQuadrature quad;
	quad.assembleQuadrature(p+1);
	double L2 = 0.0;

	for (int i=0; i<n; i++)
	{
		// for each element get L and R nodes
		double nodeL = mesh->getNode(i)[0];
		double nodeR = mesh->getNode(i+1)[0];
		for (int j=0; j<p+1; j++)
		{
			// gauss points on ref elem: zeta_i
			// gauss transformed onto actual domain: phi_i
			// for each quad point, find w_i * J * (u(x_i(zeta_i)) - u_h hat(zeta_i))^2,
			// J = h / 2
			double phi = 0.5 * (1.0 + quad.points[j]);
			double x = nodeL + phi * (nodeR - nodeL);
			double u_exact = f_analytic(x);
			double u_h = evaluateSolution({x});

			L2error += quad.weights[j] * 0.5 * (nodeR - nodeL) * pow(u_exact - u_h, 2);
		}
	}

	return pow(L2error, 0.5);
}

double FE_Solution::getL2Error(double (*f_analytic)(const std::vector<double>&))
{
	// FIXME: implement
	return 0.0;
}

void FE_Solution::sendSolutionToFile(int noGridPoints, double (*f_analytic)(double))
{
	std::ofstream file("solution.csv");
	file << "x,u,u_exact\n";

	for (int i=0; i<noGridPoints; i++)
	{
		double x = i / (noGridPoints - 1.0);
		file << x << "," << evaluateSolution({x}) << "," << f_analytic(x) << "\n";
	}

	file.close();
}

void FE_Solution::sendSolutionToFile(int noGridPoints, double (*f_analytic)(const std::vector<double>&))
{
	// FIXME: implement
}