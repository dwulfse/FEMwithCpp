#ifndef ELEMENT1DHEADERDEF
#define ELEMENT1DHEADERDEF

#include <vector>
#include "Element.hpp"

// double computeTriangleArea(const Point2D &P1, const Point2D &P2, const Point2D &P3);
// void computeAffineMatrix(const Point2D &P1, const Point2D &P2, const Point2D &P3, double A[2][2]);
// bool invertAffineMatrix(double A[2][2], double A_inv[2][2]);

class Element1D : public Element
{
	public:
		std::vector<double> nodes;
		std::vector<Element1D> elements;

		// constructors
		Element1D();
		Element1D(int p, std::vector<int>& local_DoF, std::vector<double>& nodes);

		// destructor
		virtual ~Element1D();

		// methods
		// FIXME: get methods actually do caculation so revisit
		virtual std::vector<std::vector<double>> getReferenceStiffness() override;
		virtual std::vector<std::vector<double>> getLocalStiffness() override;
		virtual std::vector<double> getLocalLoad(double (*f)(const std::vector<double>&)) override;
};

#endif