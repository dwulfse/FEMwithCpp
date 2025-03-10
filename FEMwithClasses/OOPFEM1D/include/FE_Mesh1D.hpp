#ifndef FEMESH1DHEADERDEF
#define FEMESH1DHEADERDEF

#include "FE_Mesh.hpp"

class FE_Mesh1D : public FE_Mesh
{
	public:
		std::vector<double> nodes;

		// constructors
		FE_Mesh1D(int n, int p);
		FE_Mesh1D();

		// destructor
		virtual ~FE_Mesh1D();

		// methods
		virtual std::vector<double> getNode(int i) override;
		virtual void constructMesh(std::string filename = "") override;
		virtual double evaluateSolution(std::vector<double> x, std::vector<double> solution) override;
		// virtual void allocateStiffness() override;
		// virtual void assembleStiffnessMatrix() override;
		// virtual void assembleLoadVector(double (*f)(double)) override;
		// virtual void applyBoundaryConditions(double u0, double u1, bool boundary_u0, bool boundary_u1) override;
};

#endif