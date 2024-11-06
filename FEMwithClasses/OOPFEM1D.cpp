#include <iostream>
#include <vector>

// TODO:
// - full GaussQuadrature class
// - full local stiffness implementation
// - proper h calculation in Element

class GaussQuadrature
{
	public:
		std::vector<double> points;
		std::vector<double> weights;

		GaussQuadrature()
		{
		}

		void assembleQuadrature(int n)
		{
			// temporary only 2 gauss points
			points.push_back(-1.0/sqrt(3.0));
			points.push_back(1.0/sqrt(3.0));
			weights.push_back(1.0);
			weights.push_back(1.0);
		}
};

class PolynomialSpace
{
	public:
		std::vector<double> nodes; // nodes

		PolynomialSpace(std::vector<double> nodes)
		 : nodes(nodes)
		{
		}

		double evaluate(int i, double x)
		{
			if (x >= nodes[i-1] && x <= nodes[i])
			{
				return (x - nodes[i-1]) / (nodes[i] - nodes[i-1]);
			}
			else if (x >= nodes[i] && x <= nodes[i+1])
			{
				return (nodes[i+1] - x) / (nodes[i+1] - nodes[i]);
			}
			else
			{
				return 0.0;
			}
		}

		double evaluate_deriv(int i, double x)
		{
			if (x >= nodes[i-1] && x <= nodes[i])
			{
				return 1 / (nodes[i] - nodes[i-1]);
			}
			else if (x >= nodes[i] && x <= nodes[i+1])
			{
				return -1 / (nodes[i+1] - nodes[i]);
			}
			else
			{
				return 0.0;
			}
		}
};

class Element
{
	public:
		std::vector<double> nodes; // nodes
		std::vector<int> node_indicies; // node indicies

		Element()
		{
		}

		Element(std::vector<double> nodes, std::vector<int> node_indicies)
		 : nodes(nodes), node_indicies(node_indicies)
		{
		}

		std::vector<std::vector<double>> getLocalStiffness()
		{
			const double h = 1.0/8.0;
			return {{1.0 / h, -1.0 / h}, {-1.0 / h, 1.0 / h}};
		}

		void getLocalStiffnessGOLDPLATED(double** A, PolynomialSpace phi)
		{
			GaussQuadrature quad;
			quad.assembleQuadrature(2);
			double gauss;

			for (int i=0; i<2; i++)
			{
				for (int j=0; j<2; j++)
				{
					A[i][j] = 0.0;
					for (int k=0; k<2; k++)
					{
						gauss = (quad.points[k] + 1) / 2; // transform to [0, 1]
						A[i][j] += quad.weights[k] * phi.evaluate_deriv(node_indicies[i], gauss) * phi.evaluate_deriv(node_indicies[j], gauss);
					}
				}
			}
		}

		std::vector<double> getLocalLoad(double (*f)(double))
		{
			const double h = nodes[1] - nodes[0];
			return {0.5 * f(0.5) * h, 0.5 * f(0.5) * h};
		}

		void getLocalLoadGOLDPLATED(double* b, PolynomialSpace phi, double (*f)(double))
		{
			GaussQuadrature quad;
			quad.assembleQuadrature(2);
			double gauss;

			for (int i=0; i<2; i++)
			{
				b[i] = 0.0;
				for (int k=0; k<2; k++)
				{
					gauss = (quad.points[k] + 1) / 2; // transform to [0, 1]
					b[i] += quad.weights[k] * f(gauss) * phi.evaluate(node_indicies[i], gauss);
				}
			}
		}
};

class FE_Mesh
{
	public:
		int n; // number of elements
		std::vector<double> nodes; // nodes
		std::vector<Element> elements; // elements
		std::vector<std::vector<double>> stiffness; // stiffness matrix
		std::vector<double> load; // load vector

		FE_Mesh(int n)
		 : n(n), nodes(n+1), elements(n), stiffness(n+1, std::vector<double>(n+1)), load(n+1)
		{
		}

		void constructMesh()
		{
			// nodes.reserve(n+1);
			// elements.reserve(n);
			const double h = 1.0/n;
			nodes[0] = 0.0;
			for (int i=1; i<n+1; i++)
			{
				nodes[i] = i * h;
				Element elem;
				elem.nodes = {nodes[i-1], nodes[i]};
				elem.node_indicies = {i-1, i};
				elements[i-1] = elem;
			}
		}

		std::vector<std::vector<double>> assembleStiffnessMatrix()
		{
			//PolynomialSpace phi(nodes);

			for (int k=0; k<n; k++)
			{
				std::vector<std::vector<double>> local_stiffness = elements[k].getLocalStiffness();
				std::vector<int> globalDoF = elements[k].node_indicies;
				// handle DoF
				for (int i=0; i<2; i++)
				{
					for (int j=0; j<2; j++)
					{
						stiffness[globalDoF[i]][globalDoF[j]] += local_stiffness[i][j];
					}
				}
			}

			return stiffness;
		}

		std::vector<double> assembleLoadVector(double (*f)(double))
		{
			//PolynomialSpace phi(nodes);

			for (int k=0; k<n; k++)
			{
				std::vector<double> local_load = elements[k].getLocalLoad(f);
				std::vector<int> globalDoF = elements[k].node_indicies;
				// handle DoF
				for (int i=0; i<2; i++)
				{
					load[globalDoF[i]] += local_load[i];
				}
			}

			return load;
		}
};

class Solver
{
	public:
		int n; // size of system
		std::vector<std::vector<double>> A; // stiffness
		std::vector<std::vector<double>> Q; // orthogonal matrix
		std::vector<std::vector<double>> R; // upper triangular matrix
		std::vector<double> b; // load
		std::vector<double> u; // solution vector

		Solver(int n, std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double> u)
		 : n(n), A(A), Q(n+1, std::vector<double>(n+1)), R(n+1, std::vector<double>(n+1)),  b(b), u(u)
		{
		}

		void factoriseQR()
		{
			for (int j=0; j<n; j++)
			{
				// Q[:, j] = A[:, j]
				for (int k=0; k<n; k++)
				{
					Q[k][j] = A[k][j];
				}

				for (int i=0; i<j; i++)
				{
					// R[i][j] = Q[:, i] . A[:, j]
					R[i][j] = 0.0;
					for (int k=0; k<n; k++)
					{
						R[i][j] += Q[k][i] * A[k][j];
					}
					// Q[:, j] -= R[i][j] * Q[:, i]
					for (int k=0; k<n; k++)
					{
						Q[k][j] -= R[i][j] * Q[k][i];
					}
				}

				// R[j][j] = ||Q[:, j]||
				R[j][j] = 0.0;
				for (int k=0; k<n; k++)
				{
					R[j][j] += Q[k][j] * Q[k][j];
				}
				R[j][j] = pow(R[j][j], 0.5);

				// Q[:, j] /= R[j][j]
				for (int k=0; k<n; k++)
				{
					Q[k][j] /= R[j][j];
				}
			}
		}

		std::vector<double> solveSystem()
		{

			factoriseQR();

			// Find RHS Qtb
			std::vector<double> Qtb(n, 0.0);
			for (int i=0; i<n; i++)
			{
				for (int j=0; j<n; j++)
				{
					Qtb[i] += Q[j][i] * b[j];
				}
			}

			// Solve Ru = Qtb via backward substitution
			for (int i=n-1; i>= 0; i--)
			{
				u[i] = Qtb[i];
				for (int j=i+1; j<n; j++)
				{
					u[i] -= R[i][j] * u[j];
				}
				u[i] /= R[i][i];
			}

			return u;
		}
};

class FE_Solution
{
	public:
		int n; // number of nodes
		std::vector<std::vector<double>> globalStiffness; // stiffness matrix
		std::vector<double> globalLoad; // RHS vector
		std::vector<double> solution; // solution vector

		FE_Solution(int n)
		 : n(n), globalStiffness(n+1, std::vector<double>(n+1)), globalLoad(n+1), solution(n+1)
		{
		}

		FE_Mesh initialise()
		{
			// create nodes
			// create elements
			// create mesh
			auto f = [](double x) { return 1.0; };
			FE_Mesh mesh(n);
			mesh.constructMesh();

			return mesh;
		}

		std::vector<double> solve(FE_Mesh mesh, double (*f)(double))
		{
			globalStiffness = mesh.assembleStiffnessMatrix();
			globalLoad = mesh.assembleLoadVector(f);
			Solver solver(n, globalStiffness, globalLoad, solution);
			solution = solver.solveSystem();
			return solution;
		}
};

int main()
{
	int n = 8;
	auto f = [](double x) { return 1.0; };
	FE_Solution FEM(n);
	FE_Mesh mesh = FEM.initialise();
	std::vector<double> solution = FEM.solve(mesh, f);
	std::cout << "solution: " << std::endl;
	for (int i=0; i<n; i++)
	{
		std::cout << "u(" << i * (1.0/n) << ") = " << solution[i] << std::endl;
	}

	return 0;
}