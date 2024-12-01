#include "CSR_Matrix.hpp"
#include <stdexcept>

CSRMatrix::CSRMatrix()
{
}

CSRMatrix::CSRMatrix(int nnz, int n)
 : entries(nnz), col_no(nnz), row_start(n+1)
{
}

CSRMatrix::CSRMatrix(std::vector<double> entries, std::vector<int> col_no, std::vector<int> row_start)
 : entries(entries), col_no(col_no), row_start(row_start)
{
} // change to copy constructor

double& CSRMatrix::operator()(int i, int j)
{
	if (i >= row_start.size() - 1 || i < 0 || j >= row_start.size() - 1 || j < 0)
	{
		throw std::out_of_range("Index out of range");
	}

	for (int k=row_start[i]; k<row_start[i+1]; k++)
	{
		if (col_no[k] == j)
		{
			return entries[k];
		}
	}

	throw std::out_of_range("Element not accessible");
}
