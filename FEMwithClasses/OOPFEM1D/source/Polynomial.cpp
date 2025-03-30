#include "Polynomial.hpp"
#include <cmath>

double legendre(int i, double x)
{
	if (i == 0)
	{
		return 1.0;
	}
	else if (i == 1)
	{
		return x;
	}
	else
	{
		return ((2.0 * i - 1.0) * x * legendre(i - 1, x) - (i - 1) * legendre(i - 2, x)) / i;
	}
}

double legendre_deriv(int i, double x)
{
	if (i == 0)
	{
		return 0.0;
	}
	else if (i == 1)
	{
		return 1.0;
	}
	else
	{
		if (abs(abs(x) - 1.0) < 1e-10)
		{
			return 0.5 * i * (i + 1.0) * (i % 2 == 0 ? 1.0 : -1.0);
		}
		return (i / (1.0 - x * x)) * (legendre(i - 1, x) - x * legendre(i, x)); // 1-x^2 or x^2-1??
	}
}

double lobatto(int i, double x)
{
	if (i == 0)
	{
		return (1.0 - x) / 2.0;
	}
	else if (i == 1)
	{
		return (1.0 + x) / 2.0;
	}
	else
	{
		return std::sqrt(i - 1.5) * (legendre(i, x) - legendre(i-2, x));
	}
}