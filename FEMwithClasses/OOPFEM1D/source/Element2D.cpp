#include "Element2D.hpp"
#include "GaussQuadrature2D.hpp"
#include "Helper.hpp"
#include <cmath>

// helper functions
double computeTriangleArea(const Point2D &P1, const Point2D &P2, const Point2D &P3)
{
	return 0.5 * fabs((P2.x - P1.x) * (P3.y - P1.y) - (P3.x - P1.x) * (P2.y - P1.y));
}

void computeAffineMatrix(const Point2D &P1, const Point2D &P2, const Point2D &P3, double A[2][2])
{
	A[0][0] = P2.x - P1.x;
	A[0][1] = P3.x - P1.x;
	A[1][0] = P2.y - P1.y;
	A[1][1] = P3.y - P1.y;
}

bool invertAffineMatrix(double A[2][2], double A_inv[2][2])
{
	double det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
	if (fabs(det) < 1e-12)
	{
		return false;
	}

	A_inv[0][0] = A[1][1] / det;
	A_inv[0][1] = -A[0][1] / det;
	A_inv[1][0] = -A[1][0] / det;
	A_inv[1][1] = A[0][0] / det;

	return true;
}

// constructors
Element2D::Element2D()
{
}

Element2D::Element2D(int id, int p, std::vector<int>& local_DoF, std::vector<Point2D>& nodes)
 : Element(id, p, local_DoF), nodes(nodes)
{
}

// destructor
Element2D::~Element2D()
{
}

std::vector<std::vector<double>> Element2D::getReferenceStiffness()
{
	// reference stiffness for triangle (-1, -1), (1, -1), (-1, 1)
	// |T| = 2
	// shape functions N_1 = -0.5x -0.5y, N_2 = 0.5x + 0.5, N_3 = 0.5y + 0.5
	// grad N_1 = (-0.5, -0.5), grad N_2 = (0.5, 0), grad N_3 = (0, 0.5)
	// K_ref_{ij} = |T| * (grad N_i . grad N_j)

	double gradRef[3][2] = { 
		{-0.5, -0.5},
		{ 0.5,  0.0},
		{ 0.0,  0.5}
	};

	std::vector<std::vector<double>> refStiffness(3, std::vector<double>(3, 0.0));

	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			refStiffness[i][j] = 2.0 * (gradRef[i][0] * gradRef[j][0] + gradRef[i][1] * gradRef[j][1]);
		}
	}

	return refStiffness;
}

std::vector<std::vector<double>> Element2D::getLocalStiffness()
{
	// find Jacobian of transforation to physsical triangle
	double A[2][2];
	computeAffineMatrix(nodes[0], nodes[1], nodes[2], A);
	double J = 0.5 * fabs(A[0][0] * A[1][1] - A[0][1] * A[1][0]);
	
	// get reference stiffness
	std::vector<std::vector<double>> refStiffness = getReferenceStiffness();

	// transform reference stiffness to physical triangle
	std::vector<std::vector<double>> stiffness(3, std::vector<double>(3, 0.0));
	for (int i=0; i<3; i++)
	{
		for (int j=0; j<3; j++)
		{
			stiffness[i][j] = refStiffness[i][j] * J;
		}
	}
	return stiffness;
}

std::vector<double> Element2D::getLocalLoad(double (*f)(const std::vector<double>&))
{
	std::vector<double> load(3, 0.0);

	GaussQuadrature2D quad;
	quad.assembleQuadrature(3);

	double J = 0.5 * computeTriangleArea(nodes[0], nodes[1], nodes[2]);

	Point2D P = nodes[0], Q = nodes[1], R = nodes[2];

	for (int k=0; k<3; k++)
	{
		Point2D point = mapToPhysical(P, Q, R, quad.points[k]);
		double fx = f({point.x, point.y});	// FIXME: and y

		load[0] += quad.weights[k] * fx * (-0.5 * quad.points[k].x - 0.5 * quad.points[k].y);
		load[1] += quad.weights[k] * fx * (0.5 * quad.points[k].x + 0.5);
		load[2] += quad.weights[k] * fx * (0.5 * quad.points[k].y + 0.5);
	}
	for (int i=0; i<3; i++)
	{
		load[i] *= J;
	}

	return load;
}