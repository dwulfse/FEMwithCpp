#ifndef HELPERHEADERDEF
#define HELPERHEADERDEF

#include "Point2D.hpp"

Point2D mapToPhysical(const Point2D& P1, const Point2D& P2, const Point2D& P3, const Point2D& xi);
double computeJacobian(const Point2D& P1, const Point2D& P2, const Point2D& P3);

#endif