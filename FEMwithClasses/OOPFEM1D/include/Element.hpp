#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include <vector>
#include "PolynomialSpace.hpp"
#include "GaussQuadrature.hpp"
#include "Point2D.hpp"

double computeTriangleArea(const Point2D &P1, const Point2D &P2, const Point2D &P3);
void computeAffineMatrix(const Point2D &P1, const Point2D &P2, const Point2D &P3, double A[2][2]);
bool invertAffineMatrix(double A[2][2], double A_inv[2][2]);

class Element
{
	public:
		int id; // element id
		int p; // polynomial degree
		std::vector<double> nodes; // nodes
		std::vector<int> local_DoF; // DoF
		PolynomialSpace poly;	// polynomial space

		// constructors
		Element();
		Element(std::vector<double> nodes, std::vector<int> local_DoF, int p);

		// destructor
		~Element();

		// methods
		// ** get methods actually do caculation so revisit
		std::vector<std::vector<double>> getReferenceStiffness();
		std::vector<std::vector<double>> getLocalStiffness();
		std::vector<std::vector<double>> getLocalStiffness2D();
		std::vector<double> getLocalLoad(double (*f)(double));
};

#endif