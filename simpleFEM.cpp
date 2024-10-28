#include <iostream>
#include <cmath>

// TODO:
// - Setup Git repo
// - Plot residuals on u ( in python? )
// - calculate integral using generalised quadrature
// - Setup sparse matrix structure

void print_vec(double* vec, int N)
{
	for (int i=0; i<N; i++)
	{
		std::cout << vec[i] << ",\t";
	}

	std::cout << std::endl;
}

void print_mat(double** mat, int N)
{
	for (int i=0; i<N; i++)
	{
		print_vec(mat[i], N);
	}
}

double** allocateMatrix(int N)
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
{
	return new double[N]();
}

void deallocateVector(double* b)
{
	delete[] b;
}

void assembleSystem(int N, double** A, double* b, double (*func)(double))
{
	double h = 1.0/N;

	for (int i=0; i<N; i++)
	{
		A[i][i] = 2.0/h;
		if (i < N-1)
		{
			A[i][i+1] = -1.0/h;
		}
		if (i > 0)
		{
			A[i][i-1] = -1.0/h;
		}

		b[i] = h * func(i * h);
	}
}

void factoriseQR(int N, double** A, double** Q, double** R)
{
	for (int j=0; j<N; j++)
	{
		for (int k=0; k<N; k++)
		{
			Q[k][j] = A[k][j];
		}

		for (int i=0; i<j; i++)
		{
			for (int k=0; k<N; k++)
			{
				R[i][j] += Q[k][i] * A[k][j];
			}
			for (int k=0; k<N; k++)
			{
				Q[k][j] -= R[i][j] * Q[k][i];
			}
		}

		for (int k=0; k<N; k++)
		{
			R[j][j] += Q[k][j] * Q[k][j];
		}
		R[j][j] = pow(R[j][j], 0.5);

		for (int k=0; k<N; k++)
		{
			Q[k][j] /= R[j][j];
		}
	}
}

void factoriseGSQR(int N, double** A, double** Q, double** R)
{
	double* temp = allocateVector(N);

	for (int j=0; j<N; j++)
	{
		for (int k=0; k<N; k++)
		{
			temp[k] = A[k][j];
		}

		for (int i=0; i<j; i++)
		{
			for (int k=0; k<N; k++)
			{
				R[i][j] += Q[k][i] * temp[k];
			}
			for (int k=0; k<N; k++)
			{
				temp[k] -= R[i][j] * Q[k][i];
			}
		}

		for (int k=0; k<N; k++)
		{
			R[j][j] += temp[k] * temp[k];
		}
		R[j][j] = pow(R[j][j], 0.5);

		for (int k=0; k<N; k++)
		{
			Q[k][j] = temp[k] / R[j][j];
		}
	}

	deallocateVector(temp);
}

void factoriseQR_3(int N, double** A, double** Q, double** R) {
	for (int j = 0; j < N; ++j) {
		// Compute R[j][j]
		double norm = 0.0;
		for (int i = 0; i < N; ++i) {
			norm += A[i][j] * A[i][j];
		}
		R[j][j] = std::sqrt(norm);
		
		// Compute Q[:, j]
		for (int i = 0; i < N; ++i) {
			Q[i][j] = A[i][j] / R[j][j];
		}

		// Update remaining columns of A
		for (int k = j + 1; k < N; ++k) {
			R[j][k] = 0.0;
			for (int i = 0; i < N; ++i) {
				R[j][k] += Q[i][j] * A[i][k];
			}
			for (int i = 0; i < N; ++i) {
				A[i][k] = A[i][k] - Q[i][j] * R[j][k];
			}
		}
	}
}

void solveSystem(int N, double** Q, double** R, double* b, double* u)
{
	double* Qtb = allocateVector(N);
	for (int i=0; i<N; i++)
	{
		for (int j=0; j<N; j++)
		{
			Qtb[i] += Q[j][i] * b[j];
		}
	}

	for (int i=N-1; i>= 0; i--)
	{
		u[i] = Qtb[i];
		for (int j=i+1; j<N; j++)
		{
			u[i] -= R[i][j] * u[j];
		}
		u[i] /= R[i][i];
	}

	std::cout << "Qtb:" << std::endl;
	print_vec(Qtb, N);

	deallocateVector(Qtb);
}

void findExactU(int N, double* u, double (*f_exact)(double))
{
	double h = 1.0/N;
	for (int i=0; i<N; i++)
	{
		u[i] = f_exact(i * h);
	}
}

void testLinSys()
{
	// Test linear system solver

	double** A2 = allocateMatrix(3);
	double* b2 = allocateVector(3);

	std::cout << "A:" << std::endl;
	std::cin >> A2[0][0] >> A2[0][1] >> A2[0][2] >> A2[1][0] >> A2[1][1] >> A2[1][2] >> A2[2][0] >> A2[2][1] >> A2[2][2];

	std::cout << "b:" << std::endl;
	std::cin >> b2[0] >> b2[1] >> b2[2];

	double** Q2 = allocateMatrix(3);
	double** R2 = allocateMatrix(3);
	double* u2 = allocateVector(3);

	factoriseQR(3, A2, Q2, R2);

	std::cout << "Q:" << std::endl;
	print_mat(Q2, 3);
	std::cout << "R:" << std::endl;
	print_mat(R2, 3);

	double** A2_check = allocateMatrix(3);
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			for (int k=0; k<3; k++)
			{
				A2_check[i][j] += Q2[i][k] * R2[k][j];
			}
		}
	}

	std::cout << "A_check:" << std::endl;
	print_mat(A2_check, 3);

	solveSystem(3, Q2, R2, b2, u2);

	std::cout << "u:" << std::endl;
	print_vec(u2, 3);

	deallocateMatrix(A2, 3);
	deallocateVector(b2);
	deallocateVector(u2);
	deallocateMatrix(Q2, 3);
	deallocateMatrix(R2, 3);
}

int main()
{
	// FEM parameters
	const int N = 8; 	// Number of elements (N+1 nodes)
	auto f = [](double x) { return sin(atan(1.0) * 4.0 * x); }; 	// Function to solve
	auto f_exact = [](double x) { return pow(atan(1.0) * 4.0, -1) * sin(atan(1.0) * 4.0 * x); };

	// Initialise A, b, and u
	double** A = allocateMatrix(N);
	double* b = allocateVector(N);
	double* u = allocateVector(N);

	assembleSystem(N, A, b, f);

	std::cout << "FEM System: ==============================" << std::endl;
	std::cout << "A:" << std::endl;
	print_mat(A, N);
	std::cout << std::endl << "b:" << std::endl;
	print_vec(b, N);
	std::cout << std::endl << std::endl;

	// Solve system
	double** Q = allocateMatrix(N);
	double** R = allocateMatrix(N);

	factoriseQR_3(N, A, Q, R);
	solveSystem(N, Q, R, b, u);

	std::cout << "FEM Solution: ============================" << std::endl;
	std::cout << "Q:" << std::endl;
	print_mat(Q, N);
	std::cout << std::endl << "R:" << std::endl;
	print_mat(R, N);
	std::cout << std::endl << "u:" << std::endl;
	print_vec(u, N);
	std::cout << std::endl;

	double* u_exact = allocateVector(N);
	findExactU(N, u_exact, f_exact);

	std::cout << "u_exact:" << std::endl;
	print_vec(u_exact, N);

	double* error = allocateVector(N);
	for (int i=0; i<N; i++)
	{
		error[i] = u[i] - u_exact[i];
	}

	double total_error = 0.0;
	for (int i=0; i<N; i++)
	{
		total_error += error[i] * error[i];
	}
	total_error = pow(total_error, 0.5);

	std::cout << std::endl;
	std::cout << "Error: " << total_error << std::endl << std::endl;
	std::cout << "Error vector:" << std::endl;
	print_vec(error, N);

	deallocateMatrix(A, N);
	deallocateVector(b);
	deallocateVector(u);
	deallocateMatrix(Q, N);
	deallocateMatrix(R, N);
	deallocateVector(u_exact);

	std::cout << std::endl << std::endl;
	testLinSys();

	return 0;
}