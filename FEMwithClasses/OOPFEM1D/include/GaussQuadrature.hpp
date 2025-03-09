#ifndef GAUSSQUADRATUREHEADERDEF
#define GAUSSQUADRATUREHEADERDEF

#include <vector>

class GaussQuadrature
{
	public:
		std::vector<double> points;
		std::vector<double> weights;

		// constructor
		GaussQuadrature();

		// destructor
		~GaussQuadrature();

		// n in constructor directy less clunky ?
		void assembleQuadrature(int n);
};

#endif