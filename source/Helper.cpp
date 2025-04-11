#include "Helper.hpp"
#include <cmath>

Point2D mapToPhysical(const Point2D& P1, const Point2D& P2, const Point2D& P3, const Point2D& xi)
{
	Point2D P;
	P.x = P1.x + 0.5 * ((xi.x + 1.0) * (P2.x - P1.x) + (xi.y + 1.0) * (P3.x - P1.x));
	P.y = P1.y + 0.5 * ((xi.x + 1.0) * (P2.y - P1.y) + (xi.y + 1.0) * (P3.y - P1.y));

	return P;
}

void computeAffineMatrix(const Point2D &P1, const Point2D &P2, const Point2D &P3, double A[2][2])
{
	A[0][0] = P2.x - P1.x;
	A[0][1] = P3.x - P1.x;
	A[1][0] = P2.y - P1.y;
	A[1][1] = P3.y - P1.y;
}

bool invertAffineMatrix(double A[2][2], double A_inv[2][2])
{
	double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
	if (fabs(det) < 1e-12)
	{
		return false;
	}

	A_inv[0][0] = A[1][1] / det;
	A_inv[0][1] = -A[0][1] / det;
	A_inv[1][0] = -A[1][0] / det;
	A_inv[1][1] = A[0][0] / det;

	return true;
}
