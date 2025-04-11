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

		// 1d basis
		double evaluate(int i, double x);
		double evaluate_deriv(int i, double x);
		double evaluate_lobatto(int i, double x);
		double evaluate_lobatto_deriv(int i, double x);
		// 2d basis
		void evaluate_affine(double xi1, double xi2, double lambda[3]);
		void evaluate_affine_grad(double xi1, double xi2, double grad[3][2]);
		double kernel(int k_2, double x);
		double kernel_deriv(int k_2, double x);
		double evaluate_edge(int k, double lambda_a, double lambda_b);
		double evaluate_bubble(int n1, int n2, const double lambda[3]);
		double basis_2D(int i, double x1, double x2);
		void basis_2D_grad(int i, double x1, double x2, double grad[2]);
};

#endif