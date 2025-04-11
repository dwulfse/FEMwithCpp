#ifndef ELEMENT1DHEADERDEF
#define ELEMENT1DHEADERDEF

#include <vector>
#include "Element.hpp"

class Element1D : public Element
{
	public:
		std::vector<double> nodes;
		std::vector<Element1D> elements;

		// constructors
		Element1D();
		Element1D(int p, std::vector<int>& local_DoF, std::vector<double>& nodes);

		// destructor
		virtual ~Element1D();
		
		std::vector<std::vector<double>> getReferenceStiffness();

		// methods
		virtual std::vector<std::vector<double>> getLocalStiffness() override;
		virtual std::vector<double> getLocalLoad(double (*f)(const std::vector<double>&)) override;
};

#endif