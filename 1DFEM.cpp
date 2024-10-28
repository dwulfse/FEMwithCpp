#include <iostream>
#include <cmath>
#include <fstream>

// TODO:
// - Setup Git repo
// - Setup sparse matrix structure
// X Work out why it isnt accurate *
// X Make it more accurate
// X * Boundary conditions (!!!)
// X L2 norm using quadrature (!!!)
// - Generalised Quadrature of q points
// - Read Notes 3 - General Meshes
// - Understand Notes 3 - General Meshes

// Class Structure:
// - FE_Solution
// 		- PGFEM w/ polynomial spaces per element
// - FE_Mesh
//		- Element objects w/ methods
// - Solver?

void print_vec(double* vec, int N)
{
	for (int i=0; i<N - 1; i++)
	{
		std::cout << vec[i] << ",\t";
	}

	std::cout << vec[N-1] << std::endl;
}

void print_mat(double** mat, int N)
{
	for (int i=0; i<N; i++)
	{
		print_vec(mat[i], N);
	}
}

double** allocateMatrix(int N)
// Allocate memory for matrix initialised to 0.0
{
	double** A = new double*[N];
	for (int i=0; i<N; i++)
	{
		A[i] = new double[N];
		for (int j=0; j<N; j++)
		{
			A[i][j] = 0.0;
		}
	}

	return A;
}

void deallocateMatrix(double** A, int N)
{
	for (int i=0; i<N; i++)
	{
		delete[] A[i];
	}
	delete[] A;
}

double* allocateVector(int N)
// Allocate memory for vector initialised to 0.0
{
	return new double[N]();
}

void deallocateVector(double* b)
{
	delete[] b;
}

void assembleA(int N, double** A, double (*a)(double), double (*b)(double), double (*c)(double))
// Assemble matrix A for phi_i basis hat functions
// A = int_0^1 LHS * v dx
{
	double h = 1.0 / (N - 1);

	for (int i=0; i<N; i++)
	{
		double a_x = a(i * h);
		double b_x = b(i * h);
		double c_x = c(i * h);
		
		// Assemble A based on L(u) = a(x)u'' + b(x)u' + c(x)
		A[i][i] = 2.0 * a_x / h; // a contribution
		A[i][i] += c_x * h / 3.0; // c contribution
		if (i < N-1)
		{
			A[i][i+1] = -1.0 * a_x / h; // a contribution
			A[i][i+1] += b_x / 2.0; // b contribution
		}
		if (i > 0)
		{
			A[i][i-1] = -1.0 * a_x / h; // a contribution
			A[i][i-1] -= b_x / 2.0; // b contribution
		}
	}
}

// void assembleb(int N, double* b, double (*f)(double))
// {
// 	const double h = 1.0 / (N - 1);

// 	for (int i=0; i<N; i++)
// 	{
// 		b[i] = h * f(i * h);
// 	}

// 	print_vec(b, N);
// }

void assemblebGauss(int N, double* b, double (*f)(double))
// Assemble vector b using Gauss quadrature
// b = int_0^1 f v dx
{
	const double h = 1.0 / (N - 1);
	// Setup Gauss quadrature based on Legendre polynomials (???)
	double gaussQuad[2][2] = {{-1.0/sqrt(3.0), 1.0}, {1.0/sqrt(3.0), 1.0}};
	double xi, integral, gauss;

	for (int i=0; i<N; i++)
	{
		xi = i * h;
		integral = 0.0;

		for (int j=0; j<2; j++)
		{
			gauss = xi + h * (gaussQuad[j][0] + 1.0) / 2.0;
			integral += gaussQuad[j][1] * f(gauss);
		}

		// Take average of the two Gauss points
		b[i] = h * integral / 2.0;
	}
}

void applyBoundaryConditions(int N, double** A, double* b, double u0, double u1, bool boundary_u0, bool boundary_u1)
{
	if (boundary_u0)
	// for LHS boundary, set first diag term to 1, first row to 0 and RHS to u0
	{
		for (int i=0; i<N; i++)
		{
			A[0][i] = 0.0;
		}
		A[0][0] = 1.0;
		b[0] = u0;
	}

	if (boundary_u1)
	// for RHS boundary, set last diag term to 1, last row to 0 and RHS to u1
	{
		for (int i=0; i<N; i++)
		{
			A[N-1][i] = 0.0;
		}
		A[N-1][N-1] = 1.0;
		b[N-1] = u1;
	}
}

void factoriseQR(int N, double** A, double** Q, double** R)
// Gram-Schmidt QR factorisation
{
	for (int j=0; j<N; j++)
	{
		// Q[:, j] = A[:, j]
		for (int k=0; k<N; k++)
		{
			Q[k][j] = A[k][j];
		}

		for (int i=0; i<j; i++)
		{
			// R[i][j] = Q[:, i] . A[:, j]
			for (int k=0; k<N; k++)
			{
				R[i][j] += Q[k][i] * A[k][j];
			}
			// Q[:, j] -= R[i][j] * Q[:, i]
			for (int k=0; k<N; k++)
			{
				Q[k][j] -= R[i][j] * Q[k][i];
			}
		}

		// R[j][j] = ||Q[:, j]||
		for (int k=0; k<N; k++)
		{
			R[j][j] += Q[k][j] * Q[k][j];
		}
		R[j][j] = pow(R[j][j], 0.5);

		// Q[:, j] /= R[j][j]
		for (int k=0; k<N; k++)
		{
			Q[k][j] /= R[j][j];
		}
	}
}

void solveSystem(int N, double** A, double* b, double* u)
// Solve system using QR factorisation and backward substitution
{
	// Factorise A=QR
	double** Q = allocateMatrix(N);
	double** R = allocateMatrix(N);

	factoriseQR(N, A, Q, R);

	// Find RHS Qtb
	double* Qtb = allocateVector(N);
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			Qtb[i] += Q[j][i] * b[j];
		}
	}

	// Solve Ru = Qtb via backward substitution
	for (int i=N-1; i>= 0; i--)
	{
		u[i] = Qtb[i];
		for (int j=i+1; j<N; j++)
		{
			u[i] -= R[i][j] * u[j];
		}
		u[i] /= R[i][i];
	}

	// Clean up
	deallocateMatrix(Q, N);
	deallocateMatrix(R, N);
}

// double calculateError(int N, double* u, double (*u_analytic)(double))
// // Calculate L2 norm of error ||U - Uh||_(L^2)^2 (???)
// {
// 	const double h = 1.0 / (N - 1);
// 	double L2norm = 0.0;
// 	double u_exact, xi;

// 	for (int i=0; i<N; i++)
// 	{
// 		xi = i * h;
// 		u_exact = u_analytic(xi);
// 		L2norm += pow(u_exact - u[i], 2);
// 	}

// 	return L2norm * h;
// }

double calculateL2Norm(int N, double* u, double (*u_analytic)(double))
{
	// u_exact: u(x_i(zeta))
	// u_h: u_h hat(zeta)
	// phi_0: phi_1 hat(zeta)
	// phi_1: phi_2 hat(zeta)
	// xi_0: x_(i-1)
	// xi_1: x_i
	// gaussQuad: {{zeta_1, w_1}, {zeta_2, w_2}}
	// gauss: x_i(zeta)
	double h = 1.0 / (N - 1);
	double L2norm = 0.0;
	double u_exact, u_h, phi_0, phi_1, xi_0, xi_1, gauss;
	double gaussQuad[2][2] = {{-1.0/sqrt(3.0), 1.0}, {1.0/sqrt(3.0), 1.0}};

	for (int i=0; i<N; i++)
	{
		// for each element get x_(i-1) and x_i nodes
		xi_0 = (i - 1) * h;
		xi_1 = i * h;
		for (int j=0; j<2; j++)
		{
			// for each gauss point (zeta) in element find
			// (u(x_i(zeta_i)) - u_h hat(zeta_i))^2 * 0.5 * (x_i - x_(i-1)) * w_i
			// where 0.5 * (x_i - x_(i-1)) =: dx/dzeta
			phi_0 = 0.5 * (1.0 - gaussQuad[j][0]);
			phi_1 = 0.5 * (1.0 + gaussQuad[j][0]);
			gauss = xi_0 * phi_0 + xi_1 * phi_1;

			u_exact = u_analytic(gauss);
			u_h = u[i-1] * phi_0 + u[i] * phi_1;

			L2norm += gaussQuad[j][1] * 0.5 * (xi_1 - xi_0) * pow(u_exact - u_h, 2);
		}
	}

	return pow(L2norm, 0.5);
}

void calculateResidualsToFile(int N, double* u, double (*u_analytic)(double))
// Send residual calculations to file for plotting
{
	const double h = 1.0 / (N - 1);
	double u_exact, xi, residual;

	std::ofstream file("residuals.csv");
	file << "xi,u,u_exact,residual\n";

	for (int i=0; i<N; i++)
	{
		xi = i * h;
		u_exact = u_analytic(xi);
		residual = u_exact - u[i];
		file << xi << "," << u[i] << "," << u_exact << "," << residual << "\n";
	}

	file.close();
}

double solve1D(int N, double (*f)(double), double (*u_analytic)(double), double (*a_func)(double), double (*b_func)(double), double (*c_func)(double), double u0=0.0, double u1=0.0, bool boundary_u0=true, bool boundary_u1=true)
// Solve 1D FEM problem
// int_0^1 LHS v dx = int_0^1 f v dx
{
	// Initialise
	const double h = 1.0 / (N - 1);
	double** A = allocateMatrix(N);
	double* b = allocateVector(N);
	double* u = allocateVector(N);

	// Assemble  and solve system
	assembleA(N, A, a_func, b_func, c_func);
	assemblebGauss(N, b, f);
	applyBoundaryConditions(N, A, b, u0, u1, boundary_u0, boundary_u1);

	solveSystem(N, A, b, u);

	// std::cout << "Solution: ==============================" << std::endl;
	// for (int i=0; i<N; i++)
	// {
	// 	std::cout << "u(" << i * h << ") = " << u[i] << std::endl;
	// }

	double error = calculateL2Norm(N, u, u_analytic);
	std::cout << "Error: " << error << std::endl;

	calculateResidualsToFile(N, u, u_analytic);

	// Clean up
	deallocateMatrix(A, N);
	deallocateVector(b);
	deallocateVector(u);

	return error;
}

int main()
{
	// N: number of nodes
	// f: RHS function to solve
	// u_analytic: analytic solution u(x)
	// a, b, c: coefficients for LHS of PDE
	// u0, u1: boundary conditions
	// boundary_u0, boundary_u1: are boundary conditions applied?
	const int N = 8;
	auto f = [](double x) { return 1.0; };
	auto u_analytic = [](double x) { return 0.5 * (x*x - x); };
	auto a = [](double x) { return -1.0; };
	auto b = [](double x) { return 0.0; };
	auto c = [](double x) { return 0.0; };
	double u0 = 0.0;
	double u1 = 0.0;
	bool boundary_u0 = true;
	bool boundary_u1 = true;

	// solve1D(N, f, u_analytic, a, b, c);

	// calculate error against h to file for plotting
	std::ofstream file("errors.csv");
	file << "N,error\n";

	for (int i=3; i<=100; i++)
	{
		file << 1.0 /(i - 1) << "," << solve1D(i, f, u_analytic, a, b, c) << "\n";
	}

	file.close();

	return 0;
}