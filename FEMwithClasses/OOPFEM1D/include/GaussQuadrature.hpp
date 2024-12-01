#ifndef GAUSSQUADRATUREHEADERDEF
#define GAUSSQUADRATUREHEADERDEF

#include <vector>

class GaussQuadrature
{
	public:
		std::vector<double> points;
		std::vector<double> weights;

		GaussQuadrature();
		void assembleQuadrature(int n);
};

#endif