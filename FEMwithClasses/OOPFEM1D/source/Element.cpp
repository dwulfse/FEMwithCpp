#include "Element.hpp"

// TODO:
// - nodes are never used, remove and just use p

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