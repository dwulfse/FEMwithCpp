#include "../Eigen/Sparse"
#include "../Eigen/SparseLU"
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

// TODO:
// - full GaussQuadrature class
// - full local stiffness implementation
// - proper h calculation in Element

class CSRMatrix
{
	public:
		std::vector<double> entries;
		std::vector<int> col_no;
		std::vector<int> row_start;

		CSRMatrix()
		{
		}

		CSRMatrix(int nnz, int n)
		 : entries(nnz), col_no(nnz), row_start(n+1)
		{
		}

		CSRMatrix(std::vector<double> entries, std::vector<int> col_no, std::vector<int> row_start)
		 : entries(entries), col_no(col_no), row_start(row_start)
		{
		}

		double& operator()(int i, int j)
		{
			if (i >= row_start.size() - 1 || i < 0 || j >= row_start.size() - 1 || j < 0)
			{
				throw std::out_of_range("Index out of range");
			}

			for (int k=row_start[i]; k<row_start[i+1]; k++)
			{
				if (col_no[k] == j)
				{
					return entries[k];
				}
			}

			static double zero = 0.0;
			return zero;
			// throw std::out_of_range("Element not accessible");
		}
};

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

		PolynomialSpace()
		{
		}

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
		PolynomialSpace poly;	// polynomial space

		Element()
		{
		}

		Element(std::vector<double> nodes, std::vector<int> node_indicies)
		 : nodes(nodes), node_indicies(node_indicies), poly(nodes)
		{
		}

		std::vector<std::vector<double>> getLocalStiffness()
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

		std::vector<double> getLocalLoad(double (*f)(double))
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
};

class FE_Mesh
{
	public:
		int n; // number of elements
		std::vector<double> nodes; // nodes
		std::vector<Element> elements; // elements
		CSRMatrix stiffness; // stiffness matrix
		std::vector<double> load; // load vector

		FE_Mesh(int n)
		 : n(n), nodes(n+1), elements(n), stiffness(), load(n+1)
		{
		}

		void constructMesh()
		{
			const double h = 1.0/n;

			nodes[0] = 0.0;
			for (int i=1; i<n+1; i++)
			{
				nodes[i] = i * h;
				Element elem({nodes[i-1], nodes[i]}, {i-1, i});
				// elem.nodes = {nodes[i-1], nodes[i]};
				// elem.node_indicies = {i-1, i};
				elements[i-1] = elem;
			}
		}

		void allocateStiffness()
		{
			std::vector<int> nnz_per_row(n+1);
			std::vector<int> DoF(elements[0].node_indicies.size());

			for (int i=0; i<n; i++) // for each element
			{
				DoF = elements[i].node_indicies;
				for (int j=0; j<DoF.size(); j++) // for each connection
				{
					nnz_per_row[DoF[j]] += DoF.size();
				}
			}
			// combine into one loop?
			stiffness.row_start.resize(n+2);
			for (int i=1; i<n+2; i++) // for each node (row)
			{
				stiffness.row_start[i] = stiffness.row_start[i-1] + nnz_per_row[i-1];
			}

			int nnz = stiffness.row_start[n+1];
			stiffness.entries.resize(nnz);
			stiffness.col_no.resize(nnz);

			std::vector<int> row_start_copy = stiffness.row_start;

			for (int k=0; k<n; k++)
			{
				DoF = elements[k].node_indicies;
				for (int i=0; i<DoF.size(); i++)
				{
					for (int j=0; j<DoF.size(); j++)
					{
						stiffness.col_no[row_start_copy[DoF[i]]++] = DoF[j];
					}
				}
			}
		}

		CSRMatrix assembleStiffnessMatrix()
		{
			std::vector<std::vector<double>> local_stiffness(2, std::vector<double>(2));
			std::vector<int> globalDoF;

			for (int k=0; k<n; k++)
			{
				local_stiffness = elements[k].getLocalStiffness();
				globalDoF = elements[k].node_indicies;
				// handle DoF
				for (int i=0; i<globalDoF.size(); i++)
				{
					for (int j=0; j<globalDoF.size(); j++)
					{
						stiffness(globalDoF[i], globalDoF[j]) += local_stiffness[i][j];
					}
				}
			}

			return stiffness;
		}

		std::vector<double> assembleLoadVector(double (*f)(double))
		{
			std::vector<double> local_load(2);
			std::vector<int> globalDoF;

			for (int k=0; k<n; k++)
			{
				local_load = elements[k].getLocalLoad(f);
				globalDoF = elements[k].node_indicies;
				// handle DoF
				for (int i=0; i<globalDoF.size(); i++)
				{
					load[globalDoF[i]] += local_load[i];
				}
			}

			return load;
		}

		void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1)
		{
			if (boundary_u0)
			{
				for (int i=stiffness.row_start[0]; i<stiffness.row_start[1]; i++)
				{
					stiffness.entries[i] = 0.0;
				}
				stiffness(0, 0) = 1.0;
				load[0] = u0;
			}

			if (boundary_u1)
			{
				for (int i=stiffness.row_start[n]; i<stiffness.row_start[n+1]; i++)
				{
					stiffness.entries[i] = 0.0;
				}
				stiffness(n, n) = 1.0;
				load[n] = u1;
			}
		}
};

class Solver
{
	public:
		int n; // size of system
		CSRMatrix A; // stiffness
		std::vector<std::vector<double>> Q; // orthogonal matrix
		std::vector<std::vector<double>> R; // upper triangular matrix
		std::vector<double> b; // load
		std::vector<double> u; // solution vector
		Eigen::SparseMatrix<double> EigenA;
		Eigen::VectorXd Eigenb;

		Solver(int n, CSRMatrix A, std::vector<double> b, std::vector<double> u)
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
					Q[k][j] = A(k,j);
				}

				for (int i=0; i<j; i++)
				{
					// R[i][j] = Q[:, i] . A[:, j]
					R[i][j] = 0.0;
					for (int k=0; k<n; k++)
					{
						R[i][j] += Q[k][i] * A(k,j);
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
			std::vector<double> Qtb(n);
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

		void setupEigen()
		{
			EigenA.resize(n, n);
			std::vector<Eigen::Triplet<double>> values;

			for (int i=0; i<n; i++)
			{
				for (int j=A.row_start[i]; j<A.row_start[i+1]; j++)
				{
					if (A.entries[j] != 0.0)
					{
						values.emplace_back(i, A.col_no[j], A.entries[j]);
					}
				}
			}

			EigenA.setFromTriplets(values.begin(), values.end());

			Eigenb.resize(n);
			for (int i=0; i<n; i++)
			{
				Eigenb[i] = b[i];
			}
		}

		std::vector<double> solveEigen()
		{
			Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

			solver.analyzePattern(EigenA);
			solver.factorize(EigenA);

			std::cout << "EigenA: \n" << EigenA << std::endl;
			std::cout << "A.entries: \n" << std::endl;
			for (int i=0; i<n; i++)
			{
				std::cout << "( " << A.entries[i] << ", " << A.col_no[i] << " )" << std::endl;
			}

			if (solver.info() != Eigen::Success)
			{
				throw std::runtime_error("Eigen factorisation failed");
			}

			// solver.compute(EigenA);
			Eigen::VectorXd solution = solver.solve(Eigenb); // PROBLEM IS HERE atm

			//convert Eigen::VectorXd to std::vector<double> u
			for (int i=0; i<n; i++)
			{
				u[i] = solution[i];
			}

			return u;
		}
};

class FE_Solution
{
	public:
		int n; // number of nodes
		CSRMatrix globalStiffness; // stiffness matrix
		std::vector<double> globalLoad; // RHS vector
		std::vector<double> solution; // solution vector

		FE_Solution(int n)
		 : n(n), globalStiffness(), globalLoad(n+1), solution(n+1)
		{
		}

		FE_Mesh initialise()
		{
			FE_Mesh mesh(n);
			mesh.constructMesh();

			return mesh;
		}

		std::vector<double> solve(FE_Mesh mesh, double (*f)(double), double u0=0.0, double u1=0.0, bool boundary_u0=true, bool boundary_u1=true)
		{
			mesh.allocateStiffness();
			globalStiffness = mesh.assembleStiffnessMatrix();
			globalLoad = mesh.assembleLoadVector(f);
			mesh.applyBoundaryConditions(u0, u1, boundary_u0, boundary_u1);

			Solver solver(n, globalStiffness, globalLoad, solution);
			// solution = solver.solveSystem();

			solver.setupEigen();
			solution = solver.solveEigen();

			return solution;
		}
};

int main()
{
	const int n = 8;
	auto f = [](double x) { return 1.0; };
	double u0 = 0.0;
	double u1 = 0.0;
	bool boundary_u0 = true;
	bool boundary_u1 = true;

	FE_Solution FEM(n);
	FE_Mesh mesh = FEM.initialise();
	std::vector<double> solution = FEM.solve(mesh, f, u0, u1, boundary_u0, boundary_u1);

	std::cout << "solution: " << std::endl;
	for (int i=0; i<n+1; i++)
	{
		std::cout << "u(" << i * (1.0/n) << ") = " << solution[i] << std::endl;
	}

	// CSRMatrix A(5, 13);
	// A.entries = {2.0, 2.0, 3.0, 4.0, 2.0, 5.0, 5.0, 8.0, 17.0, 10.0, 16.0, 14.0};
	// A.col_no = {0, 3, 0, 1, 2, 3, 0, 3, 4, 2, 3 ,4};
	// A.row_start = {0, 2, 6, 9, 11, 12};

	// std::cout << A(0, 0) << std::endl;

	return 0;
}