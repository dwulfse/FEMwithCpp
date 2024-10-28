#include <iostream>
#include <cmath>

// Function to perform QR factorization using Gram-Schmidt
void QRFactorization(int N, double** A, double** Q, double** R) {
    // Initialize R to zero
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            R[i][j] = 0.0;
        }
    }

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

// Function to solve the system QRu = b using back substitution
void SolveQR(int N, double** Q, double** R, double* b, double* u) {
    // First solve Q^T * y = b for y (since Q is orthogonal, Q^T = inverse(Q))
    double* y = new double[N];
    
    // Compute Q^T * b (this is the dot product of columns of Q with vector b)
    for (int i = 0; i < N; ++i) {
        y[i] = 0.0;
        for (int j = 0; j < N; ++j) {
            y[i] += Q[j][i] * b[j];  // Q^T is applied by switching indices
        }
    }

    // Debugging: Print intermediate vector y
    std::cout << "Intermediate vector y (Q^T * b):\n";
    for (int i = 0; i < N; ++i) {
        std::cout << y[i] << " ";
    }
    std::cout << std::endl;

		// Now solve R * u = y using back substitution
		for (int i = N - 1; i >= 0; --i) {
				u[i] = y[i];
				for (int j = i + 1; j < N; ++j) {
						u[i] -= R[i][j] * u[j];
				}
				u[i] /= R[i][i];
		}


    // Debugging: Print solution vector u after back substitution
    std::cout << "Solution vector u:\n";
    for (int i = 0; i < N; ++i) {
        std::cout << u[i] << " ";
    }
    std::cout << std::endl;

    delete[] y;
}


// Utility function to print a matrix
void printMatrix(int N, double** matrix) {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

// Utility function to print a vector
void printVector(int N, double* vec) {
    for (int i = 0; i < N; ++i) {
        std::cout << vec[i] << " ";
    }
    std::cout << std::endl;
}

int main() {
    int N;

    // Input matrix size
    std::cout << "Enter matrix size N: ";
    std::cin >> N;

    // Allocate memory for matrices A, Q, R and vectors b, u
    double** A = new double*[N];
    double** Q = new double*[N];
    double** R = new double*[N];
    double* b = new double[N];
    double* u = new double[N];

    for (int i = 0; i < N; ++i) {
        A[i] = new double[N];
        Q[i] = new double[N];
        R[i] = new double[N];
    }

    // Input matrix A
    std::cout << "Enter matrix A (row by row):\n";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            std::cin >> A[i][j];
        }
    }

    // Input vector b
    std::cout << "Enter vector b:\n";
    for (int i = 0; i < N; ++i) {
        std::cin >> b[i];
    }

    // Perform QR factorization
    QRFactorization(N, A, Q, R);

    // Solve QRu = b for u
    SolveQR(N, Q, R, b, u);

    // Output results
    std::cout << "\nMatrix Q:\n";
    printMatrix(N, Q);

    std::cout << "\nMatrix R:\n";
    printMatrix(N, R);

    std::cout << "\nSolution vector u:\n";
    printVector(N, u);

    // Free allocated memory
    for (int i = 0; i < N; ++i) {
        delete[] A[i];
        delete[] Q[i];
        delete[] R[i];
    }
    delete[] A;
    delete[] Q;
    delete[] R;
    delete[] b;
    delete[] u;

    return 0;
}
