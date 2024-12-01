#include "PolynomialSpace.hpp"

PolynomialSpace::PolynomialSpace()
{
}

PolynomialSpace::PolynomialSpace(std::vector<double> nodes)
 : nodes(nodes)
{
}

double PolynomialSpace::evaluate(int i, double x)
{
	if (x >= nodes[i-1] && x <= nodes[i])
	{
		return (x - nodes[i-1]) / (nodes[i] - nodes[i-1]);
	}
	else if (x >= nodes[i] && x <= nodes[i+1])
	{
		return (nodes[i+1] - x) / (nodes[i+1] - nodes[i]);
	}
	else
	{
		return 0.0;
	}
}

double PolynomialSpace::evaluate_deriv(int i, double x)
{
	if (x >= nodes[i-1] && x <= nodes[i])
	{
		return 1 / (nodes[i] - nodes[i-1]);
	}
	else if (x >= nodes[i] && x <= nodes[i+1])
	{
		return -1 / (nodes[i+1] - nodes[i]);
	}
	else
	{
		return 0.0;
	}
}