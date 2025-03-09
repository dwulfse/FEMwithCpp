#ifndef FEMESHHEADERDEF
#define FEMESHHEADERDEF

#include <vector>
#include <string>
#include "Element.hpp"
#include "CSR_Matrix.hpp"
#include "Point2D.hpp"

template <typename NodeType>
class FE_Mesh
{
	public:
		int n; // number of elements
		int p; // polynomial degree
		std::vector<NodeType> nodes;
		std::vector<Element> elements;
		CSRMatrix stiffness;
		std::vector<double> load;

		// constructors
		FE_Mesh(int n, int p);
		FE_Mesh();

		// destructor
		~FE_Mesh();

		// methods
		void parseMesh(std::string filename);
		void constructMesh();
		void allocateStiffness();
		void assembleStiffnessMatrix();
		void assembleLoadVector(double (*f)(double)); // perhaps (*f)(const NodeType&) ?
		void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1);
};

#include "FE_Mesh.tpp"

#endif