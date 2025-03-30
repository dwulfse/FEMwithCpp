#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include <vector>
#include "PolynomialSpace.hpp"
// #include "Point2D.hpp"

// double computeTriangleArea(const Point2D &P1, const Point2D &P2, const Point2D &P3);
// void computeAffineMatrix(const Point2D &P1, const Point2D &P2, const Point2D &P3, double A[2][2]);
// bool invertAffineMatrix(double A[2][2], double A_inv[2][2]);

class Element
{
	public:
		int id; // element id
		int p; // polynomial degree
		std::vector<int> local_DoF; // DoF
		PolynomialSpace poly;	// polynomial space

		// constructors
		Element() {};
		Element(int p, std::vector<int>& local_DoF) : p(p), poly(p), local_DoF(local_DoF) {};
		Element(int id, int p, std::vector<int>& local_DoF) : id(id), p(p), poly(p), local_DoF(local_DoF) {};

		// destructor
		virtual ~Element() {};

		// pure virtual methods
		// FIXME: get methods actually do caculation so revisit
		virtual std::vector<std::vector<double>> getLocalStiffness() = 0;
		virtual std::vector<double> getLocalLoad(double (*f)(const std::vector<double>&)) = 0;
};

#endif