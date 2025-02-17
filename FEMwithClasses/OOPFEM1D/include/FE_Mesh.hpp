#ifndef FEMESHHEADERDEF
#define FEMESHHEADERDEF

#include <vector>
#include "Element.hpp"
#include "CSR_Matrix.hpp"

class FE_Mesh
{
	public:
		int n; // number of elements
		int p; // polynomial degree
		std::vector<double> nodes;
		std::vector<Element> elements;
		CSRMatrix stiffness;
		std::vector<double> load;

		// constructors
		FE_Mesh(int n, int p);
		FE_Mesh();

		// destructor
		~FE_Mesh();

		// methods
		void constructMesh();
		void allocateStiffness();
		void assembleStiffnessMatrix();
		void assembleLoadVector(double (*f)(double));
		void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1);
};

#endif