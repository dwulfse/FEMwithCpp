#ifndef POLYNOMIALSPACEHEADERDEF
#define POLYNOMIALSPACEHEADERDEF

#include <vector>

class PolynomialSpace
{
	public:
		std::vector<double> nodes;

		PolynomialSpace();
		PolynomialSpace(std::vector<double> nodes);

		double evaluate(int i, double x);
		double evaluate_deriv(int i, double x);
};

#endif