#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include <vector>
#include "PolynomialSpace.hpp"
#include "GaussQuadrature.hpp"

class Element
{
	public:
		int p; // polynomial degree
		std::vector<double> nodes; // nodes
		std::vector<int> local_DoF; // DoF
		PolynomialSpace poly;	// polynomial space

		// constructors
		Element();
		Element(std::vector<double> nodes, std::vector<int> local_DoF, int p);

		// destructor
		~Element();

		// methods
		// ** get methods actually do caculation so revisit
		std::vector<std::vector<double>> getReferenceStiffness();
		std::vector<std::vector<double>> getLocalStiffness();
		std::vector<double> getLocalLoad(double (*f)(double));
};

#endif