#include <iostream>

int print_vec(double* vec, int N)
{
	for (int i = 0; i < N; i++)
	{
		std::cout << vec[i] << ",\t";
	}

	std::cout << std::endl;

	return 0;
}

int print_mat(double** mat, int N)
{
	for (int i = 0; i < N; i++)
	{
		print_vec(mat[i], N);
	}

	return 0;
}

double** assembleFEMMatrix(int N)
// Petrov-Galerkin FEM matrix for 
// b(u, v) = \int_0^1 u' v dx
// N+1 nodes
// Uh = cts. piecewise linears { phi_j } hats
// Vh = discts. piecewise consts { psi_i } indicators
{
	double** A = new double*[N];

	for (int i = 0; i < N; i++)
	{
		A[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 0.0;
		}

		A[i][i] = 1.0;
		if (i < N-1)
		{
			A[i][i+1] = -1.0;
		}
	}

	return A;
}

double* assembleFEMVector1(int N, double (*func)(double))
// l(v) = \int_0^1 f v dx
// N+1 nodes
// Vh = discont. piecewise consts { psi_i } indicators
{
	double h = 1.0/N;
	double* b = new double[N];

	for (int i = 0; i < N; i++)
	{
		b[i] = h * func(h * i + h/2);
	}

	return b;
}

double** assembleFEMMatrix2(int N)
{
	double** A = new double*[N];

	for (int i = 0; i < N; i++)
	{
		A[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 0.0;
		}

		A[i][i] = 1.0;
		if (i < N-1)
		{
			A[i][i+1] = -1.0;
		}
	}

	return A;
}

double* assembleFEMVector2(int N, double (*func)(double))
{
	double h = 1/N;
	double* b = new double[N];

	b[0] = h/2 * func(h/2);
	
	for (int i = 1; i < N; i++)
	{
		b[i] = h * func(h * i + h/2);	
	}

	return b;
}

double** assembleFEMMatrix3(int N)
// Second order Galerkin FEM for
// b(u, v) = \int_0^1 u' v' dx
// N+1 nodes
// Vh = cts. piecewise linears { phi_j } hats, 0 @ x=0
{
	double h = 1.0/N;
	double** A = new double*[N];

	for (int i = 0; i < N; i++)
	{
		A[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 0.0;
		}

		A[i][i] = 2.0/h;
		if (i < N-1)
		{
			A[i][i+1] = -1.0/h;
		}

		if (i > 0)
		{
			A[i][i-1] = -1.0/h;
		}
	}

	A[N-1][N-1] = 1.0/h;

	return A;
}

double* assembleFEMVector3(int N, double (*func)(double), double g1)
// l(v) = \int_0^1 f v dx + g1 v(1)
// N+1 nodes
// Vh = cts. piecewise linears { phi_j } hats, 0 @ x=0
{
	double h = 1.0/N;
	double* b = new double[N];

	for (int i = 0; i < N-1; i++)
	{
		b[i] = h/2 * (func(h * i + h/2) + func(h * i + 3 * h/2));
	}

	b[N-1] = h/2 * func(h * (N-1) + h/2) + g1;

	return b;
}

double** assembleFEMMatrix4(int N, double alpha, double beta)
{
	double h = 1.0/N;
	double** A = new double*[N];
	
	for (int i = 0; i < N; i++)
	{
		A[i] = new double[N];
		for (int j = 0; j < N; j++)
		{
			A[i][j] = 0.0;
		}

		A[i][i] = 2 * alpha/h;
		if (i < N-1)
		{
			A[i][i+1] = -alpha/h + beta/2;
		}
		if (i > 0)
		{
			A[i][i-1] = -alpha/h - beta/2;
		}
	}

	A[N-1][N-1] = alpha/h + beta/2;

	return A;
}

double* assembleVector_u0(int N, double (*u0_func)(double))
{
	double h = 1.0/N;
	double* b = new double[N-1];
	double* u0 = new double[N-1];

	for (int i = 0; i < N; i++)
	{
		b[i] = h/2 * (u0_func(h * i + h/2) + u0_func(h * i + 3 * h/2));
	}

	// Assemble matrix M
	// Solve M u0 = b

	return u0;
}

double** assembleMatrix_M(int N)
{
	double h = 1.0/N;
	double** M = new double*[N-1];

	for (int i = 0; i < N-1; i++)
	{
		M[i] = new double[N-1];
		for (int j = 0; j < N-1; j++)
		{
			M[i][j] = 0.0;
		}

		M[i][i] = 2.0*h/3;
		if (i < N-2)
		{
			M[i][i+1] = h/6;
		}
		if (i > 0)
		{
			M[i][i-1] = h/6;
		}
	}

	return M;
}

double** assembleMatrix_K(int N)
{
	double h = 1.0/N;
	double** K = new double*[N-1];

	for (int i = 0; i < N-1; i++)
	{
		K[i] = new double[N-1];
		for (int j = 0; j < N; j++)
		{
			K[i][j] = 0.0;
		}

		K[i][i] = 2.0/h;
		if (i < N-2)
		{
			K[i][i+1] = -1.0/h;
		}
		if (i > 0)
		{
			K[i][i-1] = -1.0/h;
		}
	}
	return 0;
}

double* heatEquationFEM(double tau, double alpha, double Ntime, int N, double (*u0_func)(double))
{
	double** u_array = new double*[Ntime];
	u_array[0] = assembleVector_u0(N, u0_func);
	for (int i = 1; i < Ntime; i++)
	{
		u_array[i] = new double[N-1];
	}

	double** M = assembleMatrix_M(N);
	double** K = assembleMatrix_K(N);

	for (int i = 1; i < Ntime; i++)
	{
		// iterative solve of M u_i = M - tau * alpha * K @ u_{i-1}
	}

	return 0;
}

int main() {
	int N = 4;
	auto f = [](double x) { return x; };

	// Test 1 ========================================
	double** A = assembleFEMMatrix(N);
	double* b = assembleFEMVector1(N, f);
	std::cout << "Matrix A:" << std::endl;
	print_mat(A, N);
	std::cout << "Vector b:" << std::endl;
	print_vec(b, N);
	std::cout << std::endl;

	// Test 2 ========================================
	double** A2 = assembleFEMMatrix2(N);
	double* b2 = assembleFEMVector2(N, f);
	std::cout << "Matrix A2:" << std::endl;
	print_mat(A2, N);
	std::cout << "Vector b2:" << std::endl;
	print_vec(b2, N);
	std::cout << std::endl;

	// Test 3 ========================================
	double g1 = 1.0;
	double** A3 = assembleFEMMatrix3(N);
	double* b3 = assembleFEMVector3(N, f, g1);
	std::cout << "Matrix A3:" << std::endl;
	print_mat(A3, N);
	std::cout << "Vector b3:" << std::endl;
	print_vec(b3, N);
	std::cout << std::endl;

	// Test 4 ========================================
	double alpha = 1.0;
	double beta = 1.0;
	double** A4 = assembleFEMMatrix4(N, alpha, beta);
	std::cout << "Matrix A4:" << std::endl;
	print_mat(A4, N);
	std::cout << std::endl;

	// Test 5 ========================================
	double tau = 0.1;
	int Ntime = 10;
	auto u0_func = [](double x) { return x; };
	double* u0_vec = assembleVector_u0(N, u0_func);
	std::cout << "Vector u0:" << std::endl;
	print_vec(u0_vec, N);
	double** M = assembleMatrix_M(N);
	std::cout << "Matrix M:" << std::endl;
	print_mat(M, N);
	double** K = assembleMatrix_K(N);
	std::cout << "Matrix K:" << std::endl;
	print_mat(K, N);
	double* u = heatEquationFEM(tau, alpha, Ntime, N, u0_func);
	std::cout << "Vector u:" << std::endl;
	print_vec(u, N);

	return 0;
}