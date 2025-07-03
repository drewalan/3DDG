/*
Basis.cpp
creates coefficients from polynomial basis functions (imported as hex doubles from a MATLAB script)
Drew Murray
7/2/16
*/

//class includes
#include "Basis.h"
#include "Matrix.h"


//libraries
#include <vector>
using std::vector;
#include <string>
using std::string;
using std::to_string;
using std::stoi;
#include <fstream>
using std::ifstream;
#include <iostream>
using std::cerr;
using std::cout;
#include <sstream>
using std::stringstream;
using std::istringstream;
#include <bitset>
using std::bitset;
#include <limits>
using std::numeric_limits;
#include <cmath>
using std::pow;
#include <cstring>
#include <mpi.h>


//taken from http://stackoverflow.com/questions/8616573/convert-64-bit-binary-string-representation-of-a-double-number-back-to-double-nu
double bitstring_to_double2(const char* p)
{
    unsigned long long x = 0;
    for (; *p; ++p)
    {
        x = (x << 1) + (*p - '0');
    }
    double d;
    memcpy(&d, &x, 8);
    return d;
}

//shockingly, an existing implementation doesn't seem to exist in C or C++
double hex2double(const string& hexstr)
{
	//convert to binary
	stringstream ss;
	ss << std::hex << hexstr;
	unsigned long long n;
	ss >> n;
	bitset<64> b(n);
	bitset<1> s(b.to_string().substr(0,1));
	bitset<11> e(b.to_string().substr(1,11));
	bitset<52> m(b.to_string().substr(12,52));

	return bitstring_to_double2(b.to_string().c_str());
}


Basis::Basis()//create basis from Jim Kelly's routines
{
	cout << " - Creating Basis from routines\n";
	//As of yet unimplemented, instead use Basis::Basis(string filename, string filename2)
}

Basis::Basis(string filename, string filename2)//create basis from saved coeffs: zeros and weights, and second is dpsi
{
	int rank = 0;	
	MPI_Comm_rank (MPI_COMM_WORLD, &rank); // get current MPI-process ID. O, 1, ...
	if(rank==0)cout << " - Creating Basis from coeffs files\n";
	///////////////////////////////////////////////////////////////////////////
	//get zeros and weights
	///////////////////////////////////////////////////////////////////////////
	ifstream file;
	string line;
	string tempString;

	//Lines in the read-in file have a patterned order
	/*
		order
		zero, zero, ...
		weight, weight, ...

		order
		zero, zero, zero, ...
		weight, weight, weight, ...

		...
	*/

	unsigned char itr = 0;//which zero within an order
	unsigned char modItr = 0;//order, zeros, weights, ...
	unsigned char order = 0;

	//check if files exists
	std::ifstream inStream;
	std::ifstream inStream2;
	inStream.exceptions(std::ifstream::failbit);
	inStream2.exceptions(std::ifstream::failbit);
	try {
		inStream.open(filename);
	}
	catch (const std::exception& e) {
		std::ostringstream msg;
		const std::exception temp = e;
		msg << "Opening file '" << filename
			<< "' failed, it either doesn't exist or is not accessible.\n";
		throw std::runtime_error(msg.str());
	}
	try {
		inStream2.open(filename2);
	}
	catch (const std::exception& e) {
		std::ostringstream msg;
		const std::exception temp = e;
		msg << "Opening file '" << filename2
			<< "' failed, it either doesn't exist or is not accessible.\n";
		throw std::runtime_error(msg.str());
	}


	file.open(filename);
	if (file.is_open())
	{
		while (getline(file, line))
		{
			switch (modItr)
			{
			case 0://order
			{
				order = stoi(line.substr(2, 100));
				break;
			}
			case 1://zeros
			{
				vector<double> temp;
				istringstream iss(line);
				for (itr = 0;itr < order;itr++)
				{
					iss >> tempString;
					temp.push_back(hex2double(tempString));
				}
				mZeros.push_back(temp);
				break;
			}
			case 2://weights
			{
				vector<double> temp;
				istringstream iss(line);
				for (itr = 0;itr < order;itr++)
				{
					iss >> tempString;
					temp.push_back(hex2double(tempString));
				}
				mWeights.push_back(temp);
				break;
			}
			}

			modItr = (modItr + 1) % 3;
		}
		file.close();
	}
	else cerr << "Unable to open file " << filename << "\n";


	///////////////////////////////////////////////////////////////////////////
	//get dpsi
	///////////////////////////////////////////////////////////////////////////
	ifstream file2;
	line = "";
	tempString = "";

	//Lines in the read-in file have a patterned order
	/*
		order
		dpsi, dpsi, ...
		dpsi, dpsi, ...

		order
		dpsi, dpsi, dpsi, ...
		dpsi, dpsi, dpsi, ...
		dpsi, dpsi, dpsi, ...

		...
	*/

	unsigned int colNum = 0;
	order = 0;
	bool first = true;
	vector<vector<double>> tempGrid;//DPsi is a 2d grid

	file2.open(filename2);
	if (file2.is_open())
	{
		while (getline(file2, line))
		{
			if (line.substr(0, 1) == "p")//next order
			{

				if (!first)
				{//don't push empty grid on first loop
					//mPsi.push_back(tempGrid);
					mDPsi.push_back(tempGrid);
					mPsi.push_back(tempGrid);
				}
				first = false;

				order = stoi(line.substr(2, 100));
				tempGrid.clear();
			}
			else
			{
				vector<double> tempLine;
				istringstream iss(line);
				for (colNum = 0;colNum < order;colNum++)
				{
					iss >> tempString;
					tempLine.push_back(hex2double(tempString));
				}
				tempGrid.push_back(tempLine);
			}
		}
		mPsi.push_back(tempGrid);//don't forget to push final grid after last loop
		mDPsi.push_back(tempGrid);//don't forget to push final grid after last loop
		file2.close();
	}
	else cerr << "Unable to open file " << filename2 << "\n";
}

void Basis::legendre_poly(double n, double x, vector<double>&  p0, vector<double>&  p1, vector<double>&  p2, vector<double>&  p00){
	p1[0] = 0;
	p1[1] = 0;
	p1[2] = 0;
	p0[0] = 1;
	p0[1] = 0;
	p0[2] = 0;

	double a, b;
	//Construct Nth Order Legendre Polynomial
	for (int j = 1; j <= n; j++){
		p2[0] = p1[0];
		p2[1] = p1[1];
		p2[2] = p1[2];
		p1[0] = p0[0];
		p1[1] = p0[1];
		p1[2] = p0[2];

		a = (2.0*j - 1.0) / (1.0*j);
		b = (1.0*j - 1.0) / (1.0*j);
		p0[0] = a*x*p1[0] - b*p2[0];
		p0[1] = a*(p1[0] + x*p1[1]) - b*p2[1];
		p0[2] = a*(2.0*p1[1] + x*p1[2]) - b*p2[2];

		a = (2.0*j + 1.0) / (1.0*j + 1.0);
		b = (1.0*j) / (1.0*j + 1.0);
		p00[0] = a*x*p0[0] - b*p1[0];
		p00[1] = a*(p0[0] + x*p0[1]) - b*p1[1];
		p00[2] = a*(2.0*p0[1] + x*p0[2]) - b*p1[2];
	}
}

void Basis::legendre_gauss_lobatto(int ngl, unsigned short int order){
	int kmax = 20;
	double pi;
	double x, dx;
	pi = 4.0 * atan(1.0);
	int n = ngl - 1;
	int nh = (n + 1) / 2;
	double tmp[] = { 0, 0, 0 };
	std::vector<double> mP0(tmp, tmp + 3);
	std::vector<double> mP1(tmp, tmp + 3);
	std::vector<double> mP2(tmp, tmp + 3);
	std::vector<double> mP00(tmp, tmp + 3);

	for (int i = 0; i <= order; i++){
		mZeros[order][i] = 0;
		mWeights[order][i] = 0;
	}

	// First Find Half of the Roots
	for (int i = 1; i <= nh; i++){
		x = cos((2.0*i - 1.0) / (2.0*n + 1.0)*pi);
		for (int k = 1; k <= kmax; k++){
			// Construct Legendre Polynomial and Derivatives
			Basis::legendre_poly(n, x, mP0, mP1, mP2, mP00);

			// Get next Newton Iterative
			dx = -(1.0 - x*x)*mP0[1] / (-2.*x*mP0[1] + (1.0 - x*x)*mP0[2]);
			x = x + dx;
			//cout << mP0[1] << " " << mP0[2] << " " << x << " " << dx<<" " <<fabs(dx) << " after\n";
			if (fabs(dx) < 1.0e-20) break;
		}

		mZeros[order][n + 1 - i] = x;//mXgl[n + 1 - i] = x;
		mWeights[order][n + 1 - i] = 2.0 / (1.0*(n*(n + 1))*mP0[0] * mP0[0]);// mWgl[n + 1 - i] = 2.0 / (1.0*(n*(n + 1))*mP0[0] * mP0[0]);
	}


	//Check for Zero
	if (n + 1 != 2 * nh){
		x = 0;
		Basis::legendre_poly(n, x, mP0, mP1, mP2, mP00);
		mZeros[order][nh] = x;//mXgl[nh] = x;
		mWeights[order][nh] = 2.0 / (1.0*(n*(n + 1))*mP0[0] * mP0[0]);//mWgl[nh] = 2.0 / (1.0*(n*(n + 1))*mP0[0] * mP0[0]);
	}


	//Find Remainder of Roots via Symmetry
	for (int i = 1; i <= nh; i++){
		mZeros[order][i - 1] = -mZeros[order][n + 1 - i];//mXgl[i-1] = -mXgl[n + 1 - i];
		mWeights[order][i - 1] = +mWeights[order][n + 1 - i];//mWgl[i - 1] = +mWgl[n + 1 - i];
	}

}

void Basis::legendre_gauss(int ngl, unsigned short int order){
	int kmax = 20;
	double pi = 4.0*atan(1.0);
	int n = ngl - 1;
	int nh = (n + 1) / 2;
	double x, dx;
	double tmp[] = { 0, 0, 0 };
	std::vector<double> mP0(tmp, tmp + 3);
	std::vector<double> mP1(tmp, tmp + 3);
	std::vector<double> mP2(tmp, tmp + 3);
	std::vector<double> mP00(tmp, tmp + 3);
	for (int i = 1; i <= nh; i++){
		x = cos((2.0*i - 1.0) / (2.0*n + 1.0)*pi);
		for (int k = 1; k <= kmax; k++){
			// Construct Legendre Polynomial and Derivatives
			Basis::legendre_poly(n, x, mP0, mP1, mP2, mP00);
			// Get next Newton Iterative
			dx = -mP00[0] / mP00[1];
			x = x + dx;
			if (abs(dx) < 1.0e-20) break;
		}
		mZeros[order][n + 1 - i] = x;// mXgl[n + 1 - i] = x;
		mWeights[order][n + 1 - i] = 2.0 / ((1.0 - x*x)*mP00[1] * mP00[1]);//mWgl[n + 1 - i] = 2.0 / ((1.0 - x*x)*mP00[1] * mP00[1]);
	}

	//Check for Zero
	if (n + 1 != 2 * nh){
		x = 0;
		Basis::legendre_poly(n, x, mP0, mP1, mP2, mP00);
		mZeros[order][nh] = x;//mXgl[nh] = x;
		mWeights[order][nh] = 2.0 / ((1.0 - x*x)*mP00[1] * mP00[1]);//mWgl[nh] = 2.0 / ((1.0 - x*x)*mP00[1] * mP00[1]);
	}
	//Find Remainder of Roots via Symmetry
	for (int i = 1; i <= nh; i++){
		mZeros[order][i - 1] = -mZeros[order][n + 1 - i];//mXgl[i-1] = -mXgl[n + 1 - i];
		mWeights[order][i - 1] = +mWeights[order][n + 1 - i];//mWgl[i-1] = +mWgl[n + 1 - i];
	}
}

void Basis::legendre_basis(int ngl, unsigned short int order){
	int n = ngl - 1;
	double xj, ksi;
	double tmp[] = { 0, 0, 0 };
	std::vector<double> mP0j(tmp, tmp + 3);
	std::vector<double> mP1j(tmp, tmp + 3);
	std::vector<double> mP2j(tmp, tmp + 3);
	std::vector<double> mP00(tmp, tmp + 3);
	std::vector<double> mP0i(tmp, tmp + 3);
	std::vector<double> mP1i(tmp, tmp + 3);
	std::vector<double> mP2i(tmp, tmp + 3);

	for (int i = 0; i <= order; i++){
		for (int j = 0; j <= order; j++){
			mDPsi[order][i][j] = 0;
			mPsi[order][i][j] = 0;
		}
	}

	for (int j = 1; j <= ngl; j++){
		xj = mZeros[order][j - 1];// mXgl[j - 1];
		Basis::legendre_poly(n, xj, mP0j, mP1j, mP2j, mP00);
		for (int i = 1; i <= ngl; i++){
			ksi = mZeros[order][i - 1];//mXgl[i - 1];
			Basis::legendre_poly(n, ksi, mP0i, mP1i, mP2i, mP00);

			if (i == j) mPsi[order][i - 1][j - 1] = 1;
			else mPsi[order][i - 1][j - 1] = 0;

			if (i == j) {
				if (i != 1 && i != ngl) mDPsi[order][i - 1][j - 1] = 0;
				else if (i == 1) mDPsi[order][i - 1][j - 1] = -1.0*(n*(n + 1)) / 4.0;
				else if (i == ngl) mDPsi[order][i - 1][j - 1] = +1.0*(n*(n + 1)) / 4.0;
			}
			else if (i != j) {
				mDPsi[order][i - 1][j - 1] = mP0j[0] / (mP0i[0] * (xj - ksi));
			}
		}
	}
}
