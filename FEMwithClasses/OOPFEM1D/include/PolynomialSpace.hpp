#ifndef POLYNOMIALSPACEHEADERDEF
#define POLYNOMIALSPACEHEADERDEF

#include <vector>

class PolynomialSpace
{
	public:
		int p; // polynomial degree

		// constructors
		PolynomialSpace();
		PolynomialSpace(int p);

		// destructor
		~PolynomialSpace();

		double evaluate(int i, double x);
		double evaluate_deriv(int i, double x);
		double evaluate_lobatto(int i, double x);
		double evaluate_lobatto_deriv(int i, double x);
};

#endif