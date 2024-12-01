#ifndef FEMESHHEADERDEF
#define FEMESHHEADERDEF

#include <vector>
#include "Element.hpp"
#include "CSR_Matrix.hpp"

class FE_Mesh
{
	public:
		int n; // number of elements
		std::vector<double> nodes; // nodes
		std::vector<Element> elements; // elements
		CSRMatrix stiffness; // stiffness matrix
		std::vector<double> load; // load vector

		FE_Mesh(int n);
		FE_Mesh();

		void constructMesh();
		void allocateStiffness();
		CSRMatrix assembleStiffnessMatrix();
		std::vector<double> assembleLoadVector(double (*f)(double));
		void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1);
};

#endif