#include "Polynomial.hpp"
#include <cmath>

double legendre(int i, double x)
{
	if (i == 0)
	{
		return 1.0;
	}
	if (i == 1)
	{
		return x;
	}

	double L_nm2 = 1.0;
	double L_nm1 = x;
	double L = 0.0;

	for (int n = 2; n <= i; n++)
	{
		L = ((2.0 * n - 1.0) * x * L_nm1 - (n - 1) * L_nm2) / n;
		L_nm2 = L_nm1;
		L_nm1 = L;
	}

	return L;
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
		if (abs(abs(x) - 1.0) < 1e-10) // if x on the boundary
		{
			return 0.5 * i * (i + 1.0) * (i % 2 == 0 ? 1.0 : -1.0);
		}
		return (i / (1.0 - x * x)) * (legendre(i - 1, x) - x * legendre(i, x));
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
		return (1.0 / sqrt(2.0 * (2.0 * i - 1.0))) * (legendre(i-2, x) - legendre(i, x));
	}
}