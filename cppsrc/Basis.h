/*
Basis.h
creates coefficients from polynomial basis functions
Drew Murray
7/2/16
*/

//circular build protection
#ifndef BASIS_H
#define BASIS_H

//libraries
#include <vector>
using std::vector;
#include <string>
using std::string;
#include "Matrix.h"

class Basis
{
	public:
		vector<double>* getZeros(unsigned short int order){return &mZeros[order];}
		vector<double>* getWeights(unsigned short int order){ return &mWeights[order]; }
		vector<vector<double>>* getDPsi(unsigned short int order){ return &mDPsi[order]; }
		vector<vector<double>>* getPsi(unsigned short int order){ return &mPsi[order]; }		
		Basis(string filename, string filename2);
		void legendre_poly(double n, double x, vector<double>&  p0, vector<double>&  p1, vector<double>&  p2, vector<double>&  p00);
		void legendre_gauss_lobatto(int ngl, unsigned short int order);
		void legendre_gauss(int ngl, unsigned short int order);
		void legendre_basis(int ngl, unsigned short int order);


	private:
		vector<vector<vector<double>>> mPsi;
		vector<vector<vector<double>>> mDPsi;
		vector<vector<double>> mZeros;
		vector<vector<double>> mWeights;
		Basis();//disable default constructor
		Basis(const Basis& b);//disable default constructor
};
#endif
