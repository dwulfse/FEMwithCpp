#ifndef CSRMATRIXHEADERDEF
#define CSRMATRIXHEADERDEF

#include <vector>

class CSRMatrix
{
	public:
		std::vector<double> entries;
		std::vector<int> col_no;
		std::vector<int> row_start;

		// constructors
		CSRMatrix();
		CSRMatrix(int nnz, int n);
		CSRMatrix(std::vector<double> entries, std::vector<int> col_no, std::vector<int> row_start);

		// destructor
		~CSRMatrix();

		// override () operator for read-only indexing
		double& operator()(int i, int j);
};

#endif