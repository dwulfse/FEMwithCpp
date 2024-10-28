#include <iostream>

struct CSRMatrix
{
	double* vals;
	double* col_indicies;
	double* row_indicies;

	CSRMatrix(double* v, double* c, double* r)
	: vals(v), col_indicies(c), row_indicies(r)
	{

	}
};

int findNNZ(double** A, int n, int m)
{
	int NNZ = 0;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<m; j++)
		{
			if (A[i][j] != 0)
			{
				NNZ++;
			}
		}
	}

	return NNZ;
}

CSRMatrix convertToCSR(double** A, int n, int m, int NNZ)
{
	double* vals = new double[NNZ];
	double* col_indicies = new double[NNZ];
	double* row_indicies = new double[n+1];

	int val_ptr = 0;
	int col_ptr = 0;
	int row_ptr = 0;
	bool row_found = false;

	for (int i=0; i<n; i++)
	{
		row_found = false;
		for (int j=0; j<m; j++)
		{
			if (A[i][j] != 0)
			{
				vals[val_ptr++] = A[i][j];
				col_indicies[col_ptr++] = j;
				if (!row_found)
				{
					row_indicies[row_ptr++] = val_ptr;
					row_found = true;
				}
			}
		}
	}

	return CSRMatrix(vals, col_indicies, row_indicies);
}

int main()
{
	const int n = 8;
	double** test_mat = new double*[n];
	for (int i=0; i<n; i++)
	{
		test_mat[i] = new double[n];
		for (int j=0; j<n; j++)
		{
			test_mat[i][j] = 0.0;
			if ((i + j) % 5 == 3)
			{
				test_mat[i][j] = 4;
			}
			if ((i * j + 1) % 7 == 0)
			{
				test_mat[i][j] = 2;
			}
			if ((i - j) % 6 == 2)
			{
				test_mat[i][j] = 7;
			}
		}
	}

	std::cout << "Matrix: " << std::endl;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
		{
			std::cout << test_mat[i][j] << " ";
		}
		std::cout << std::endl;
	}

	int NNZ = findNNZ(test_mat, n, n);
	CSRMatrix test_csr = convertToCSR(test_mat, n, n, NNZ);

	std::cout << "CSR:\nvals\tcol\trow" << std::endl;
	for (int i=0; i<NNZ; i++)
	{
		std::cout << test_csr.vals[i] << "\t" << test_csr.col_indicies[i] << "\t";
		if (i < n)
		{
			std::cout << test_csr.row_indicies[i];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	return 0;
}