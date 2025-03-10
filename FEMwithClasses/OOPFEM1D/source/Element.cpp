#include "Element.hpp"
#include <cmath>

// TODO:
// - nodes are never used, remove and just use p

// helper functions
double computeTriangleArea(const Point2D &P1, const Point2D &P2, const Point2D &P3)
{
	return 0.5 * fabs((P2.x - P1.x) * (P3.y - P1.y) - (P3.x - P1.x) * (P2.y - P1.y));
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

// constructors
Element::Element()
{
}

Element::Element(std::vector<double> nodes, std::vector<int> local_DoF, int p)
 : nodes(nodes), local_DoF(local_DoF), poly(p), p(p)
{
}

// destructor
Element::~Element()
{
}

std::vector<std::vector<double>> Element::getReferenceStiffness()
{
	std::vector<std::vector<double>> refStiffness(p+1, std::vector<double>(p+1, 0.0));
	GaussQuadrature quad; // change constructor to accept p+1 directly
	quad.assembleQuadrature(p+1);
	
	double gauss;
	for (int k=0; k<quad.points.size(); k++) // can just be p+1 ??
	{
		gauss = quad.points.at(k);
		for (int i=0; i<p+1; i++)
		{
			for (int j=0; j<p+1; j++)
			{
				refStiffness.at(i).at(j) += quad.weights.at(k) * poly.evaluate_lobatto_deriv(i, gauss) * poly.evaluate_lobatto_deriv(j, gauss);
			}
		}
	}

	return refStiffness;
}

std::vector<std::vector<double>> Element::getLocalStiffness()
{
	std::vector<std::vector<double>> stiffness(p+1, std::vector<double>(p+1, 0.0));
	std::vector<std::vector<double>> refStiffness = getReferenceStiffness();
	double J = (nodes.back() - nodes.front()) / 2.0;

	for (int i=0; i<p+1; i++)
	{
		for (int j=0; j<p+1; j++)
		{
			stiffness.at(i).at(j) = refStiffness.at(i).at(j) / J;
		}
	}

	return stiffness;
}

std::vector<std::vector<double>> Element::getLocalStiffness2D()
{
	// std::vector<std::vector<double>> stiffness(3, std::vector<double>(3, 0.0));
	// double area = computeTriangleArea({nodes.at(0), 0, 0}, {nodes.at(1), 0, 0}, {nodes.at(2), 0, 0});
	// FIXME: implement
	return std::vector<std::vector<double>>();
}

std::vector<double> Element::getLocalLoad(double (*f)(double))
// "globalize" quadrature cos irts just nodes and weights, not specific per element
// TODO:
{
	std::vector<double> load(p+1, 0.0);
	GaussQuadrature quad;
	quad.assembleQuadrature(p+1);
	double J = (nodes.back() - nodes.front()) / 2.0;
	double gauss, physical;

	for (int k=0; k<quad.points.size(); k++)
	{
		gauss = quad.points.at(k);
		physical = (1.0 - gauss) / 2.0 * nodes.front() + (1.0 + gauss) / 2.0 * nodes.back();
		for (int i=0; i<p+1; i++)
		{
			load.at(i) += quad.weights.at(k) * J * f(physical) * poly.evaluate_lobatto(i, gauss);
		}
	}

	return load;
}