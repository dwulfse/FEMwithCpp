struct Quadrature
{
	double* points;
	double* weights;

	Quadrature(double* p, double* w)
	: points(p), weights(w)
	{

	}
};

Quadrature generateGaussQuad(int n)
{
	
}