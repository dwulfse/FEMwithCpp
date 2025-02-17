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

		// ** change to accept n in constructor
		// ** actually find gauss points and weights, not just precomputed.
		void assembleQuadrature(int n);
};

#endif