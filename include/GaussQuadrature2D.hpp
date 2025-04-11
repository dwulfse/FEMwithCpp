#ifndef GAUSSQUADRATURE2DHEADERDEF
#define GAUSSQUADRATURE2DHEADERDEF

#include <vector>
#include "Point2D.hpp"

class GaussQuadrature2D
{
	public:
		std::vector<Point2D> points;
		std::vector<double> weights;

		// constructor
		GaussQuadrature2D();

		// destructor
		~GaussQuadrature2D();

		// Assemble triangle gauss quadrature by tensor product and map
		void assembleQuadrature(int n);
};

#endif