#include "Element.hpp"

Element::Element()
{
}

Element::Element(std::vector<double> nodes, std::vector<int> node_indicies)
 : nodes(nodes), node_indicies(node_indicies), poly(nodes)
{
}

std::vector<std::vector<double>> Element::getLocalStiffness()
{
	std::vector<std::vector<double>> stiffness(nodes.size(), std::vector<double>(nodes.size()));
	GaussQuadrature quad;
	quad.assembleQuadrature(2);
	double gauss;

	for (int i=0; i<2; i++)
	{
		for (int j=0; j<2; j++)
		{
			stiffness[i][j] = 0.0;
			for (int k=0; k<2; k++)
			{
				gauss = (quad.points[k] + 1) / 2; // transform to [0, 1]
				// stiffness[i][j] += quad.weights[k] * poly.evaluate_deriv(node_indicies[i], gauss) * poly.evaluate_deriv(node_indicies[j], gauss);
				stiffness[i][j] += quad.weights[k] * poly.evaluate_deriv(i, gauss) * poly.evaluate_deriv(j, gauss);
			}
		}
	}

	return stiffness;
}

std::vector<double> Element::getLocalLoad(double (*f)(double))
{
	std::vector<double> load(nodes.size());
	GaussQuadrature quad;
	quad.assembleQuadrature(2);
	double gauss;

	for (int i=0; i<2; i++)
	{
		load[i] = 0.0;
		for (int k=0; k<2; k++)
		{
			gauss = (quad.points[k] + 1) / 2; // transform to [0, 1]
			// load[i] += quad.weights[k] * f(gauss) * poly.evaluate(node_indicies[i], gauss);
			load[i] += quad.weights[k] * f(gauss) * poly.evaluate(i, gauss);
		}
	}

	return load;
}