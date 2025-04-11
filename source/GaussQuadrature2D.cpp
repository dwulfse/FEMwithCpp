#include "GaussQuadrature2D.hpp"
#include "GaussQuadrature.hpp"

// constructor
GaussQuadrature2D::GaussQuadrature2D()
 : points(), weights()
{
}

// destructor
GaussQuadrature2D::~GaussQuadrature2D()
{
}

void GaussQuadrature2D::assembleQuadrature(int n)
{
	// create 1D gauss quadrature
	GaussQuadrature quad;
	quad.assembleQuadrature(n);

	// resize vectors
	points.resize(n * n);
	weights.resize(n * n);

	int index = 0;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
		{
			double u = 0.5 * (quad.points[i] + 1.0);
			double v = 0.5 * (quad.points[j] + 1.0);
			Point2D point;
			point.x = 2.0 * u * (1.0 - v) - 1.0;
			point.y = 2.0 * v - 1.0;
			points[index] = point;
			weights[index] = quad.weights[i] * quad.weights[j] * (1.0 - v);
			index ++;
		}
	}
}
