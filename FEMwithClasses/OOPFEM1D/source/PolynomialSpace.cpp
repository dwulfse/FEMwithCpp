#include "PolynomialSpace.hpp"
#include "Polynomial.hpp"

#include <cmath>

// TODO:
// - work out where p is used
// - change to precompute and store (nodes and denominators??)

// constructors
PolynomialSpace::PolynomialSpace()
{
}

PolynomialSpace::PolynomialSpace(int p)
 : p(p)
{
}

// destructor
PolynomialSpace::~PolynomialSpace()
{
}

double PolynomialSpace::evaluate(int i, double x)
{
	if (i == 0)
	{
		return (1.0 - x) / 2.0;
	}
	else if (i == 1)
	{
		return (1.0 + x) / 2.0;
	}
	return 0.0;
}

double PolynomialSpace::evaluate_deriv(int i, double x)
{
	if (i == 0)
	{
		return 1.0 / 2.0;
	}
	else if (i == 1)
	{
		return -1.0 / 2.0;
	}
	return 0.0;
}

double PolynomialSpace::evaluate_lobatto(int i, double x)
{
	return lobatto(i, x);
}

double PolynomialSpace::evaluate_lobatto_deriv(int i, double x)
{
	if (i == 0)
	{
		return -0.5;
	}
	else if (i == 1)
	{
		return 0.5;
	}
	else
	{
		return std::sqrt(i - 1.5) * (legendre_deriv(i, x) - legendre_deriv(i-2, x));
	}
}