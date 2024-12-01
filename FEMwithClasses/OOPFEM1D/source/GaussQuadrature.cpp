#include "GaussQuadrature.hpp"
#include <stdexcept>
#include <cmath>

// TODO:
// - massive overhaul needed

GaussQuadrature::GaussQuadrature()
{
}

void GaussQuadrature::assembleQuadrature(int n)
{
	if (n == 1)
	{
		points = {0.0};
		weights = {2.0};
	}
	else if (n == 2)
	{
		points = {-1/sqrt(3), 1/sqrt(3)};
		weights = {1.0, 1.0};
	}
	else if (n == 3)
	{
		points = {-sqrt(3.0/5), 0.0, sqrt(3.0/5)};
		weights = {5.0/9, 8.0/9, 5.0/9};
	}
	else if (n == 4)
	{
		points = {-sqrt(3.0/7 + 2.0/7*sqrt(6.0/5)), -sqrt(3.0/7 - 2.0/7*sqrt(6.0/5)), sqrt(3.0/7 - 2.0/7*sqrt(6.0/5)), sqrt(3.0/7 + 2.0/7*sqrt(6.0/5))};
		weights = {(18 - sqrt(30))/36, (18 + sqrt(30))/36, (18 + sqrt(30))/36, (18 - sqrt(30))/36};
	}
	else
	{
		throw std::invalid_argument("Invalid number of quadrature points");
	}
}