#ifndef HELPERHEADERDEF
#define HELPERHEADERDEF

#include "Point2D.hpp"

Point2D mapToPhysical(const Point2D& P1, const Point2D& P2, const Point2D& P3, const Point2D& xi);
void computeAffineMatrix(const Point2D &P1, const Point2D &P2, const Point2D &P3, double A[2][2]);
bool invertAffineMatrix(double A[2][2], double A_inv[2][2]);
double computeTriangleArea(Point2D P1, Point2D P2, Point2D P3);
double computeTriangleSignedArea(Point2D P1, Point2D P2, Point2D P3);

#endif