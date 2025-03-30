#ifndef POINTHEADERDEF
#define POINTHEADERDEF

struct Point2D
{
	int id;
	double x, y;

	Point2D() : id(0), x(0.0), y(0.0) {};
	Point2D(int id, double x, double y) : id(id), x(x), y(y) {};
};

#endif