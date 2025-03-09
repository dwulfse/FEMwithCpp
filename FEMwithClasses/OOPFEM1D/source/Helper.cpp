#include "Helper.hpp"
#include <cmath>

Point2D mapToPhysical(const Point2D& P1, const Point2D& P2, const Point2D& P3, const Point2D& xi)
{
	Point2D P;
	P.x = P1.x + xi.x * (P2.x - P1.x) + xi.y * (P3.x - P1.x);
	P.y = P1.y + xi.x * (P2.y - P1.y) + xi.y * (P3.y - P1.y);
	return P;
}

double computeJacobian(const Point2D& P1, const Point2D& P2, const Point2D& P3)
{
	return fabs((P2.x - P1.x) * (P3.y - P1.y) - (P3.x - P1.x) * (P2.y - P1.y));
}