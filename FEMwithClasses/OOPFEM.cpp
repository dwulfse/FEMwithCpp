#include <iostream>
#include <vector>

// Project TODO:
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
// ~ FE_Mesh
//    - mesh coords
//		- Element objects w/ methods
// 				- mapping, jacobian
//     		- quadrature
//    - boundary conditions?
// ~ FE_Solution
// 		- PGFEM w/ polynomial spaces per element (CG/PG?)
// 		- solution vector
//		- polynomial spaces - element dependent (class)
// 		- DoF handler (routine)
// ~ Solver
//    - Matrix
//    - direct Linear solver
//    - nonlinear solver?

// Resources:
// - Meshes: Lecture Notes 3 - General Meshes
// - Quadrature: Coursework Question on Legendre Polynomials

class FE_Mesh
{
	public:
		int n; 	// number of elements
		std::vector<Node> nodes; // coordinates of n+1 nodes
		std::vector<Element> elements; // elements of mesh
		std::vector<std::vector<double>> connectivity; // adjacency matrix (??)
		std::vector<std::vector<double>> global_stiffness; // global stiffness matrix

		FE_Mesh(int n, std::vector<Node> nodes, std::vector<Element> elements)
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
		int order; // order of polynomial space
		std::vector<double> coords; // coordinates of nodes
		std::vector<Node> nodes; // node ids
		PolynomialSpace space; // polynomial space
		double** stiffness; // local stiffness matrix
		double** load; // local load vector

		Element(int id, int order, std::vector<double> coords, std::vector<Node> nodes, PolynomialSpace space)
		 : id(id), order(order), coords(coords), nodes(nodes), space(space)
		{
		}

		void getLocalStiffness()
		{
			int n = nodes.size();
			stiffness = new double*[n];
			for (int i=0; i<n; i++)
			{
				stiffness[i] = new double[n];
				for (int j=0; j<n; j++)
				{
					stiffness[i][j] = space.evaluateDeriv(i, coords[i]) * space.evaluateDeriv(j, coords[j]) * 0.5 * (coords[j] - coords[i]);
				}
			}
		}

		void getLocalLoad()
		{

		}

		~Element()
		{
		}
};

class Node
{
	public:
		int id; // node id
		int dim; // number of coordinates ??
		std::vector<double> coords; // coordinates of node

		Node(int id, int dim, std::vector<double> coords)
		 : id(id), dim(dim), coords(coords)
		{
		}

		~Node()
		{
		}
};

class PolynomialSpace
{
	public:
		int order; // order of polynomial space
		std::vector<Node> nodes; // nodes

		PolynomialSpace(int order)
		 : order(order)
		{
		}

		double evaluateFunc(int i, double x)
		{
			double result = 1.0;
			for (int j=0; j<order+1; j++)
			{
				if (i != j)
				{
					result *= (x - nodes[j].coords[0]) / (nodes[i].coords[0] - nodes[j].coords[0]);
				}
			}
			return result;
		}

		double evaluateDeriv(int i, double x)
		{
			double result = 0.0;
			for (int j=0; j<order+1; j++)
			{
				if (i != j)
				{
					double term = 1.0 / (nodes[i].coords[0] - nodes[j].coords[0]);
					for (int k=0; k<order+1; k++)
					{
						if (k != i && k!= j)
						{
							term *= (x - nodes[i].coords[0]) / (nodes[i].coords[0] - nodes[k].coords[0]);
						}
					}
				}
			}
			return result;
		}

		~PolynomialSpace()
		{
		}
};

class GaussQuadrature
{
	public:
		int n; // number of points
		std::vector<double> points; // quadrature points
		std::vector<double> weights; // quadrature weights

		GaussQuadrature(int n, std::vector<double> points, std::vector<double> weights)
		 : n(n), points(points), weights(weights)
		{
		}

		~GaussQuadrature()
		{
		}
};

class FE_Solution
// DoF handler
// polynomial spaces (per element)
{
	public:
		int n; // number of elements
		// int dim; // dimension of problem
		// std::vector<Node> nodes; // nodes
		// std::vector<Element> elements; // elements
		// FE_Mesh mesh; // mesh
		std::vector<double> u; // solution vector
		// std::vector<PolynomialSpace> space; // polynomial space

		FE_Solution(int n)
		 : n(n)
		{
		}

		void initialise()
		{
			// create nodes

			// create elements
			// create mesh

		}
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
	int n = 8;
	int dim = 1;

	FE_Solution FEM(n);
	// FEM.initialise();

	return 0;
}