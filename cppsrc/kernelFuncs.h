/*
kernelFuncs.h
functions that implement the conservation law through the flux
Drew Murray
8/5/24
*/

#ifndef KERNELFUNC
#define KERNELFUNC

#include "Matrix.h"
#include "Element.h"

Matrix momentum(Matrix u, Matrix dRho, float c, float Rho0, float BOverA);
float momentum(float u, float dRho, float c, float Rho0, float BOverA);

Matrix mass(Matrix u, Matrix dRho, float c, float Rho0, float BOverA);
float mass(float u, float dRho, float c, float Rho0, float BOverA);

#endif
