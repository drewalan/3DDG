/*
kernelFuncs.cpp
functions that implement the conservation law through the flux
Drew Murray
8/5/24
*/

#include "kernelFuncs.h"
#include "Matrix.h"
#include "Element.h"

//TODO move nonlinear if in mass to operator
Matrix momentum(Matrix u, Matrix dRho, float c, float Rho0, float BOverA)
{
	return (dRho + dRho % dRho * BOverA * .5F * (1.0F / Rho0)) * (c * c / Rho0);//% is overloaded to be elementwise-multiplication
}
float momentum(float u, float dRho, float c, float Rho0, float BOverA)
{
	return (dRho + dRho * dRho * BOverA * .5F * (1.0F / Rho0)) * (c * c / Rho0);
}
Matrix mass(Matrix u, Matrix dRho, float c, float Rho0, float BOverA)
{
	if (BOverA == 0.0F) return u * Rho0;
	else return u % (dRho.addScalar(Rho0));//% is overloaded to be elementwise-multiplication
}
float mass(float u, float dRho, float c, float Rho0, float BOverA)
{
	if (BOverA == 0.0F) return u * Rho0;
	else return u * (dRho + Rho0);
}

