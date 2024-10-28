#include <iostream>

// TODO:
// - Setup Git repo
// ~ Setup sparse matrix structure
// 		X CSR Matrix
// 		X Full to CSR
// 		- CSR to Full
// - Generalised Quadrature of q points *****
// X Read Notes 3 - General Meshes
// - Understand Notes 3 - General Meshes
// - Decide dynamic allocation vs static

// Class Structure:
// - FE_Mesh
//    - mesh coords
//		- Element objects w/ methods
// 				- mapping, jacobian
//     		- quadrature
//    - boundary conditions?
// - FE_Solution
// 		- PGFEM w/ polynomial spaces per element (CG/PG?)
// 		- solution vector
//		- polynomial spaces - element dependent (class)
// 		- DoF handler (routine)
// - Solver
//    - Matrix
//    - direct Linear solver
//    - nonlinear solver?

// Resources:
// - Meshes: Lecture Notes 3 - General Meshes
// - Quadrature: Coursework Question on Legendre Polynomials

class FE_Mesh
{
	public:
		int n; 	// number of nodes
		double* nodes; // coordinates of nodes
		Element* elements; // elements of mesh
		double** connectivity; // adjacency matrix (??)
		double** global_stiffness; // global stiffness matrix

		FE_Mesh(int n, double* nodes, Element* elements)
		 : n(n), nodes(nodes), elements(elements)
		{
		}

		void constructConnectivity()
		{

		}

		void assembleGlobalStiffness()
		{

		}

		void applyBoundaryConditions()
		{

		}

		~FE_Mesh()
		{
			delete[] nodes;
			delete[] elements;
			for (int i=0; i<n; i++)
			{
				delete[] connectivity[i];
			}
			delete[] connectivity;
		}
};

class Element
// TODO:
// - map from local to actual
// - jacobian
// - node class?
{
	public:
		int id; // element id
		double* coords; // coordinates of nodes
		double** stiffness; // local stiffness matrix
		double** load; // local load vector

		Element(int id, double* coords)
		 : id(id), coords(coords)
		{
		}

		void getLocalStiffness()
		{

		}

		void getLocalLoad()
		{

		}

		~Element()
		{
			delete[] coords;
		}
};

class PolynomialSpace
{
	public:
		int order; // order of polynomial space
		// function
		// derivative

		PolynomialSpace(int order)
		 : order(order)
		{
		}

		~PolynomialSpace()
		{
		}
};

class FE_Solution
{

};

class Solver
{
	public:
		int n; // size of system
		double** A; // stiffness
		double* b; // load
		double* u; // solution vector

		Solver(double** A, double* b, double* u)
		 : A(A), b(b), u(u)
		{
		}

		void factoriseQR()
		{

		}

		void solveSystem()
		{

		}

		~Solver()
		{
			delete[] b;
			delete[] u;
			for (int i=0; i<n; i++)
			{
				delete[] A[i];
			}
			delete[] A;
		}
};

int main()
{
	return 0;
}