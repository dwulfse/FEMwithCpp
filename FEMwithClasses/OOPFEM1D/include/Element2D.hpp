#ifndef ELEMENT2DHEADERDEF
#define ELEMENT2DHEADERDEF

#include <vector>
#include "Element.hpp"
#include "Point2D.hpp"

double computeTriangleArea(const Point2D &P1, const Point2D &P2, const Point2D &P3);
void computeAffineMatrix(const Point2D &P1, const Point2D &P2, const Point2D &P3, double A[2][2]);
bool invertAffineMatrix(double A[2][2], double A_inv[2][2]);

class Element2D : public Element
{
	public:
		std::vector<Point2D> nodes;

		// constructors
		Element2D();
		Element2D(int id, int p, std::vector<int>& local_DoF, std::vector<Point2D>& nodes);

		// destructor
		virtual ~Element2D();

		// methods
		// FIXME: get methods actually do caculation so revisit
		virtual std::vector<std::vector<double>> getReferenceStiffness() override;
		virtual std::vector<std::vector<double>> getLocalStiffness() override;
		virtual std::vector<double> getLocalLoad(double (*f)(const std::vector<double>&)) override;
};

#endif