#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include <vector>
#include "PolynomialSpace.hpp"
#include "GaussQuadrature.hpp"

class Element
{
	public:
		std::vector<double> nodes; // nodes
		std::vector<int> node_indicies; // node indicies
		PolynomialSpace poly;	// polynomial space

		Element();
		Element(std::vector<double> nodes, std::vector<int> node_indicies);

		std::vector<std::vector<double>> getLocalStiffness();
		std::vector<double> getLocalLoad(double (*f)(double));
};

#endif