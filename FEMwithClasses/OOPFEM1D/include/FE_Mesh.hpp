#ifndef FEMESHHEADERDEF
#define FEMESHHEADERDEF

#include <vector>
#include <string>
#include "Element.hpp"
#include "CSR_Matrix.hpp"

class FE_Mesh
{
	public:
		int n; // number of elements
		int p; // polynomial degree
		std::vector<Element> elements;
		CSRMatrix stiffness;
		std::vector<double> load;

		// constructors
		FE_Mesh(int n, int p) : n(n), p(p), stiffness(), load(p*n+1, 0.0), elements(n) {};
		FE_Mesh() {};

		// destructor
		virtual ~FE_Mesh() {};

		// pure virtual method
		virtual std::vector<double> getNode(int i) = 0;
		virtual void constructMesh(std::string filename = "") = 0;
		virtual double evaluateSolution(std::vector<double> x, std::vector<double> solution) = 0;

		// virtual methods
		virtual void allocateStiffness();
		virtual void assembleStiffnessMatrix();
		virtual void assembleLoadVector(double (*f)(double));
		virtual void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1);
};

#endif