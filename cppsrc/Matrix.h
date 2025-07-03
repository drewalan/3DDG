/*
Matrix.h
minimally library dependent matrix operations
adapted from
https://github.com/aloispichler/Matrix-Class
© 2020, Alois Pichler

No corresponding .cpp to prevent needing to modify these functions at all.
Not a problem as this file does not #include anything other than standard c libraries
*/


#ifndef matrix_H		
#define matrix_H


#include <cstring>		// enable memcpy
#include <cassert>		// enable assert
#include <iomanip>		// std::cout, std::endl
#include <iostream>
#include <cmath>

//suppress warnings from Matrix.h only https://stackoverflow.com/questions/6321839/how-to-disable-warnings-for-particular-include-files
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wreorder"

using std::isinf;

class Matrix;			// forward declaration (see below)


template<typename T = double>				// the default type is double
class Vector			// declaration of the class Vector
{
public:
	Vector<T>(unsigned _dim = 0)				// default constructor
	{
		vLen = _dim; vData = new T[vLen];
	}
	Vector<T>(const Vector<T>&);			// copy constructor
	~Vector<T>() { delete[] vData; }		// destructor

	template<unsigned dim>					// construct from array
	Vector<T>(const T(&vec)[dim]) : Vector<T>(dim)	// v inherits from Vector
	{
		if (vData != vec)	// Vector<double> vec({ -1, 7, 1});
			std::memcpy(vData, vec, dim * sizeof(T));
	}

	Vector<T>(const Matrix&) {};				// construct from Matrix

	Vector<T>& operator= (const Vector<T>&);// assignment operator=
	T& operator[] (unsigned) const;			// element access operator
	T& operator() (unsigned) const;			// the indexes are one-based, not zero based.

	unsigned Length() const { return vLen; }	// count elements
	Vector<T> Zeros();						// fill Vector with 0
	Vector<T> operator+ (T);				// add a scalar
	Vector<T> operator+ (const Vector<T>&);	// Vector addition
	Vector<T> operator- (const Vector<T>&);	// Vector subtraction
	Vector<T> operator* (T);				// multiplication by scalar
	unsigned maxPosition();					// maximum [] position

private:
	T* vData; unsigned vLen;
};


class Matrix
{
public:

	Matrix() : mRows(0), mCols(0), mData(nullptr) {};	// default constructor
	Matrix(const unsigned);					// constructor of a square matrix
	Matrix(const unsigned, const unsigned);	// constructor of a rectangular matrix
	Matrix Let(std::initializer_list<double>);
	Matrix(const Matrix&);					// copy constructor
	~Matrix() { delete[] mData; }			// destructor
	Matrix(const Vector<double>&);			// construct from Vector


	template<unsigned mRows, unsigned mCols>	// construct from array
	Matrix(const double(&arr)[mRows][mCols]) : Matrix(mRows, mCols)
	{
		std::memcpy(mData, arr, mRows * mCols * sizeof(double));
	}
	//	Matrix A((double[][3]){{ -1.6, 7, 13}, {5, 9.3, 2}}); enable array initialization

	Matrix& operator= (const Matrix&);		// assignment operator=

	double* operator[] (unsigned) const;	// element access operator
	double& operator()(unsigned, unsigned) const;//	the indexes are ONE-based, NOT zero based.
	unsigned rows() const { return mRows; }	// the number of rows
	unsigned cols() const { return mCols; }	// the number of columns

	Matrix operator+ (const Matrix&);		// matrix addition
	Matrix operator+ (double);				// add a multiple of the unit matrix
	Matrix operator- (const Matrix&);		// matrix subtraction
	Matrix operator* (double);				// multiplication by scalar
	Matrix operator* (const Matrix&);		// matrix multiplication
	Matrix operator% (const Matrix&);		// DREW elementwise mulitplication
	Matrix addScalar(double);				// DREW elementwise addition with scalar
	Matrix operator/ (Matrix);				// Solve A/b
	Vector<double> operator/ (Vector<double>);	// Solve A/b

	Matrix Fill(const double);

private:
	double* mData; unsigned mRows, mCols;
};


static Matrix Eye(const unsigned);			// identity matrix
Matrix Transpose(Matrix);					// transpose matrix
Matrix Inverse(Matrix, double);				// pseudoinverse
std::ostream& operator << (std::ostream&, const Matrix);


struct QResult	// container to hold result of QR decomposition
{
	Matrix Mat;		// the matrix Mat holds Q and R
	//unsigned Rank;	// the rank of the matrix //TODO changed to prevent warning
	unsigned Rank=0;	// the rank of the matrix
	Vector<double>    V1, V2;	// auxiliary vectors for the decomposition
	Vector<unsigned>  Permutation;		// permutation to pivot rows
};


template<typename T>	//	C++ lacks Square
const T Square(const T& a)
{
	return a * a;
}



template<typename T>		// copy constructor
inline Vector<T>::Vector(const Vector<T>& src)
{
	vLen = src.vLen;
	vData = new T[vLen];	// deep copy
	std::memcpy(vData, src.vData, vLen * sizeof(T));
}


template<>				// construct a Vector from Matrix
//inline Vector<double>::Vector<double>(const Matrix& src) //TODO changed to prevent warning
inline Vector<double>::Vector(const Matrix& src)
{
	vLen = src.rows() * src.cols();
	vData = new double[vLen];	// deep copy from Vector source
	std::memcpy(vData, src[0], vLen * sizeof(double));
}


template<typename T>		// assignment operator=
inline Vector<T>& Vector<T>::operator=(const Vector<T>& src)
{
	if (this != &src)		// protect against invalid self-assignment
	{
		if (vLen != src.vLen)
		{
			delete[] vData;				// destroy old stack?
			vLen = src.vLen;
			vData = new T[src.vLen];
		}	// deep copy
		std::memcpy(vData, src.vData, vLen * sizeof(T));
	}
	return *this;
}


template<typename T> 		// element access operator
inline T& Vector<T>::operator[](const unsigned _index) const
{
	return vData[_index];
}


template<typename T>		// element access operator
inline T& Vector<T>::operator()(const unsigned index) const
{
	return vData[index - 1];
}// the indexes are one-based, not zero based.


template<typename T>
inline Vector<T> Vector<T>::Zeros()//	set all entries to 0
{
	for (unsigned i = 0; i < vLen; i++)
		(*this)[i] = (T)0;
	return *this;
}


template<typename T>		// find the max [] position in the vector
inline unsigned Vector<T>::maxPosition()
{
	unsigned j = 0;
	for (unsigned i = 1; i < vLen; i++)
		if (vData[i] > vData[j]) j = i;
	return j;
}


template<typename T>
inline Vector<T> Vector<T>::operator+ (T lambda)	// add a scalar to this vector
{
	Vector<T> c(vLen);
	for (unsigned i = 0; i < vLen; i++)
		c.vData[i] = this->vData[i] + lambda;
	return c;
}


template<typename T>
inline Vector<T> Vector<T>::operator+ (const Vector<T>& other)	// add Vectors
{
	assert((vLen == other.Length()) && "Vector length mismatch. ");
	Vector c(vLen);
	for (unsigned i = 0; i < vLen; i++)
		c.vData[i] = this->vData[i] + other.vData[i];
	return c;
}


template<typename T>
inline Vector<T> Vector<T>::operator- (const Vector<T>& other)	// subtract Vectors
{
	assert((vLen == other.Length()) && "Vector length mismatch. ");
	Vector c(vLen);
	for (unsigned i = 0; i < vLen; i++)
		c.vData[i] = (*this).vData[i] - other.vData[i];
	return c;
}


template<typename T>
inline Vector<T> Vector<T>::operator* (T lambda)	// multiply this vector by scalar
{
	Vector c(vLen);
	for (unsigned i = 0; i < vLen; i++)
		c.vData[i] = lambda * this->vData[i];
	return c;
}


inline double Norm(const Vector<double> vec, const double p = 2)	// p-norm of a Vector
{
	if (p < 1) std::cout << "Norm: 1 > p= " << p << std::endl;
	double tmpNorm = 0;
	if (p == 1)			//	1-norm, weighted
	{
		for (unsigned i = 0; i < vec.Length(); i++)
			tmpNorm += fabs(vec[i]);
		return tmpNorm / vec.Length();
	}
	else if (p == 2)	//	2-norm, weighted
	{
		for (unsigned i = 0; i < vec.Length(); i++)
			tmpNorm += Square(vec[i]);
		return sqrt(tmpNorm / vec.Length());
	}
	else if (isinf(p))	//	inf-norm
	{
		for (unsigned i = 0; i < vec.Length(); i++)
			tmpNorm = std::max(tmpNorm, fabs(vec[i]));
		return tmpNorm;
	}
	else				//	all other weighted p-norms
	{
		for (unsigned i = 0; i < vec.Length(); i++)
			if (vec[i] != 0.) tmpNorm += pow(fabs(vec[i]), p);
		return pow(tmpNorm / vec.Length(), 1 / p);
	}
}


template<typename T>		// print
inline std::ostream& operator<< (std::ostream& os, const Vector<T>& src)
{
	os << "(" << src.Length() << "-vector):";
	for (unsigned i = 0; i < src.Length(); i++)
		os << " " << src[i];
	return os;
}

inline QResult decomposeQRQ(Matrix& A, double aTol = 1e-6)	//	QRQ 
{
	QResult QRQ; QRQ.Mat = A;			// in-place decomposition. check: Is this A?
	unsigned maxRank = std::min(QRQ.Mat.cols(), QRQ.Mat.rows());	// maximal rank;
	QRQ.Rank = 0;
	QRQ.Permutation = Vector<unsigned>(maxRank);
	for (unsigned k = 0; k < QRQ.Permutation.Length(); k++)
		QRQ.Permutation[k] = k;				// initialize P
	QRQ.V1 = Vector<double>(QRQ.Mat.rows());	// auxiliary vector

	Vector<double> normRow2(QRQ.Mat.rows());	// squares of the row norms
	for (unsigned i = 0; i < QRQ.Mat.rows(); i++)	// compute the row norms
	{
		normRow2[i] = 0.0;
		for (unsigned j = 0; j < QRQ.Mat.cols(); j++)// compute the norms
			normRow2[i] += Square(QRQ.Mat[i][j]);
	}

	for (unsigned k = 0; k < maxRank; k++)	// Householder iteration step
	{
		unsigned rMax = k;					// Pivot: find largest row
		for (unsigned i = k + 1; i < QRQ.Mat.rows(); i++)
			if (normRow2[i] > normRow2[rMax]) rMax = i;

		if (k != rMax)		//	row pivoting: swap rows
		{
			QRQ.Permutation[k] = rMax;
			std::swap(normRow2[k], normRow2[rMax]);
			for (unsigned j = 0; j < QRQ.Mat.cols(); j++)
				std::swap(QRQ.Mat[k][j], QRQ.Mat[rMax][j]);
		}
		if (normRow2[k] < aTol * aTol)	// numerical check
		{
			//std::cout << "QR: Rank deficient " << A.rows() << "x" << QRQ.Mat.cols() << " matrix. Rank: " << QRQ.Rank << std::endl;
			break;
		}
		++QRQ.Rank;	// increase the matrix rank

		double normTmp = sqrt(normRow2[k]);		// start Householder Q1
		if (QRQ.Mat[k][k] > 0)
		{
			QRQ.V1[k] = QRQ.Mat[k][k] + normTmp; QRQ.Mat[k][k] = -normTmp;
			normTmp = sqrt(2 * normTmp * QRQ.V1[k]);
		}
		else	// norm of Householder vector
		{
			QRQ.V1[k] = QRQ.Mat[k][k] - normTmp; QRQ.Mat[k][k] = normTmp;
			normTmp = -sqrt(-2 * normTmp * QRQ.V1[k]);
		}
		QRQ.V1[k] /= normTmp;		// normalize with V1 > 0
		for (unsigned j = k + 1; j < QRQ.Mat.cols(); j++)
			QRQ.Mat[k][j] /= normTmp;// normalize the reflector

		for (unsigned i = k + 1; i < QRQ.Mat.rows(); i++)	// Householder reflection-----
		{
			double sum = QRQ.V1[k] * QRQ.Mat[i][k];		// compute the inner product
			for (unsigned j = k + 1; j < QRQ.Mat.cols(); j++)
				sum += QRQ.Mat[k][j] * QRQ.Mat[i][j];

			sum *= 2; QRQ.Mat[i][k] -= sum * QRQ.V1[k];	// apply reflector
			for (unsigned j = k + 1; j < QRQ.Mat.cols(); j++)
				QRQ.Mat[i][j] -= sum * QRQ.Mat[k][j];
			normRow2[i] -= Square(QRQ.Mat[i][k]);
		}	// update all norms
	}

	unsigned Delta = QRQ.Mat.rows() - QRQ.Rank;	// rank deficiency
	if (Delta)		// apply smart Householder
	{
		QRQ.V2 = Vector<double>(QRQ.Rank);			// auxiliary vector
		for (unsigned k = QRQ.Rank - 1; k < -1; --k)	// start Householder Q2
		{
			double normTmp = 0;
			for (unsigned i = 0; i <= Delta; ++i)
				normTmp += Square(QRQ.Mat[k + i][k]);
			normTmp = sqrt(normTmp); double* s = &QRQ.Mat[k + Delta][k];
			if (*s > 0)			// sign of Householder transformation
			{
				QRQ.V2[k] = *s + normTmp; *s = -normTmp;
				normTmp = sqrt(2 * normTmp * QRQ.V2[k]);
			}
			else
			{
				QRQ.V2[k] = *s - normTmp; *s = normTmp;
				normTmp = -sqrt(-2 * normTmp * QRQ.V2[k]);
			}
			QRQ.V2[k] /= normTmp;				// now, V2 > 0
			for (unsigned i = 0; i < Delta; ++i)
				QRQ.Mat[k + i][k] /= normTmp;		// normalize the reflector


			for (unsigned j = k - 1; j < -1; --j)	// all Householder reflections-------
			{
				double sum = QRQ.V2[k] * QRQ.Mat[k + Delta][j];	// compute the inner product
				for (unsigned i = 0; i < Delta; ++i)
					sum += QRQ.Mat[k + i][k] * QRQ.Mat[k + i][j];

				sum *= 2; QRQ.Mat[k + Delta][j] -= QRQ.V2[k] * sum;	// apply reflector
				for (unsigned i = 0; i < Delta; ++i)
					QRQ.Mat[k + i][j] -= QRQ.Mat[k + i][k] * sum;
			}
		}
	}
	return QRQ;
}		/* return the decomposition	*/





inline Matrix SolveQRQ(QResult QRQ, Matrix b)
{
	assert((QRQ.Mat.rows() == b.rows()) && "SolveQRQ. Dimensions do not match.");
	Matrix xRes = Matrix(QRQ.Mat.cols(), b.cols());// the sizes of x and b differ!
	for (unsigned c = 0; c < b.cols(); c++)		// cover all columns of b
	{
		for (unsigned i = 0; i < QRQ.Permutation.Length(); i++)
			if (i != QRQ.Permutation[i])			// Step 1: Permutation P
				std::swap(b[i][c], b[QRQ.Permutation[i]][c]);

		unsigned Delta = QRQ.Mat.rows() - QRQ.Rank;
		if (Delta)	// Householder Q2, if applicable
			for (unsigned k = QRQ.Rank - 1; k < -1; --k)	// apply all Householders
			{
				double sum = QRQ.V2[k] * b[k + Delta][c];	// inner product Q2
				for (unsigned i = 0; i < Delta; ++i)
					sum += QRQ.Mat[k + i][k] * b[k + i][c];
				sum *= 2;
				b[k + Delta][c] -= QRQ.V2[k] * sum;	// apply Q2
				for (unsigned i = 0; i < Delta; ++i)
					b[k + i][c] -= QRQ.Mat[k + i][k] * sum;
			}

		for (unsigned i = 0; i < xRes.rows(); i++)	// Step 4: L^{-1}
			if (i < QRQ.Rank)				// compute xRes[i]
			{
				double* s = &xRes[i][c]; *s = b[i + Delta][c];
				for (unsigned j = 0; j < i; j++)	// substitute backwards
					*s -= QRQ.Mat[i + Delta][j] * xRes[j][c];
				*s /= QRQ.Mat[i + Delta][i];
			}
			else
				xRes[i][c] = 0;

		for (unsigned k = QRQ.Rank - 1; k < -1; k--)	// Step 5: all Householders, Q1
		{
			double sum = QRQ.V1[k] * xRes[k][c];	// the inner product
			for (unsigned j = k + 1; j < QRQ.Mat.cols(); j++)
				sum += QRQ.Mat[k][j] * xRes[j][c];

			sum *= 2; xRes[k][c] -= sum * QRQ.V1[k];	// apply Q1
			for (unsigned j = k + 1; j < QRQ.Mat.cols(); j++)
				xRes[j][c] -= sum * QRQ.Mat[k][j];
		}
	}
	return xRes;
}


inline Matrix Inverse(Matrix A, double aTol = 1e-6)	//	pseudoinverse
{
	return SolveQRQ(decomposeQRQ(A, aTol), Eye(A.rows()));
}
// class MatrixException
// {	public:
// 		std::string msg;
// 		MatrixException(std::string arg) : msg(arg) {}
// };


inline Matrix::Matrix(unsigned dim)					//	Construct a square matrix
{
	mRows = dim; mCols = dim; mData = new double[dim * dim];
}


inline Matrix::Matrix(unsigned dimR, unsigned dimC)	// construct a rectangular matrix
{
	mRows = dimR; mCols = dimC; mData = new double[dimR * dimC];
}


inline Matrix::Matrix(const Matrix& src)	// copy constructor
{
	mRows = src.mRows; mCols = src.mCols;
	mData = new double[mRows * mCols];		// deep copy
	std::memcpy(mData, src.mData, mRows * mCols * sizeof(double));
}


inline Matrix::Matrix(const Vector<double>& src)	// construct from Vector
{
	mRows = src.Length(); mCols = 1;
	mData = new double[mRows];	// deep copy from Vector source
	std::memcpy(mData, &src[0], mRows * sizeof(double));
}


inline Matrix& Matrix::operator=(const Matrix& src)	// assignment operator=
{
	if (this != &src)		// protect against invalid self-assignment
	{
		if (mRows * mCols != src.mRows * src.mCols)
		{
			delete[] mData;			// destroy old stack
			mData = new double[src.mRows * src.mCols];
		}
		mRows = src.mRows; mCols = src.mCols;	// deep copy
		std::memcpy(mData, src.mData, mRows * mCols * sizeof(double));
	}
	return *this;
}


inline double* Matrix::operator[](const unsigned _row) const 	// row access operator
{
	assert((_row < mRows) && "Matrix: Row out of range.");
	return &mData[mCols * _row];
}		//	row-major


inline double& Matrix::operator()(const unsigned row, const unsigned col) const 	//	element access operator
{
	assert((row > 0 && row <= mRows) && "Matrix: row out of range.");
	assert((col > 0 && col <= mCols) && "Matrix: col out of range.");
	return mData[(row - 1) * mCols + col - 1];
}	//	the indexes are one-based, not zero based.


inline Matrix Matrix::Fill(const double val = 0)	 //	fill this matrix with entries val
{
	for (unsigned i = 0; i < mRows; i++)
		for (unsigned j = 0; j < mCols; j++)
			(*this)[i][j] = val;
	return *this;
}


inline Matrix Eye(const unsigned mRows)	// provides the identity matrix
{
	Matrix c(mRows, mRows);
	for (unsigned i = 0; i < mRows; i++)
	{
		for (unsigned j = 0; j < mRows; j++) c[i][j] = 0;
		c[i][i] = 1;
	}
	return c;
}


inline Matrix Matrix::operator+ (double lambda)	// add lambda* identity
{
	Matrix c(mRows, mCols);
	unsigned k = 0;
	for (unsigned i = 0; i < mRows; i++)
		for (unsigned j = 0; j < mCols; j++)
		{
			c.mData[k] = this->mData[k]; if (i == j) c.mData[k] += lambda;
			k++;
		}
	return c;
}


inline Matrix Matrix::operator+ (const Matrix& other)	// add matrices
{
	assert((other.rows() == mRows && other.cols() == mCols) && "Matrices can't be added.");
	Matrix c(mRows, mCols);
	for (unsigned i = 0; i < mRows * mCols; i++)
		c.mData[i] = (*this).mData[i] + other.mData[i];
	return c;
}

inline Matrix Matrix::addScalar (double scalar)	// DREW elementwise multiply matrices
{
	Matrix c(mRows, mCols);
	for (unsigned i = 0; i < mRows * mCols; i++)
		c.mData[i] = (*this).mData[i] + scalar;
	return c;
}

inline Matrix Matrix::operator% (const Matrix& other)	// DREW elementwise multiply matrices
{
	assert((other.rows() == mRows && other.cols() == mCols) && "Matrices can't be elementwise multiplied.");
	Matrix c(mRows, mCols);
	for (unsigned i = 0; i < mRows * mCols; i++)
		c.mData[i] = (*this).mData[i] * other.mData[i];
	return c;
}


inline Matrix Matrix::operator- (const Matrix& other)	// subtract matrices
{
	assert((other.rows() == mRows && other.cols() == mCols) && "Matrices cannot be substracted.");
	Matrix diff(mRows, mCols);
	for (unsigned i = 0; i < mRows * mCols; i++)
		diff.mData[i] = (*this).mData[i] - other.mData[i];
	return diff;
}


inline Matrix Matrix::operator* (double lambda)	// multiplication by scalar
{
	Matrix c(mRows, mCols);
	for (unsigned i = 0; i < mRows; i++)
		for (unsigned j = 0; j < mCols; j++)
			c.mData[i * mCols + j] = lambda * this->mData[i * mCols + j];
	return c;
}


inline Matrix Transpose(Matrix A)					// transposition
{
	Matrix ATranspose(A.cols(), A.rows());
	for (unsigned i = 0; i < ATranspose.rows(); i++)
		for (unsigned j = 0; j < ATranspose.cols(); j++)
			ATranspose[i][j] = A[j][i];
	return ATranspose;
}


inline Matrix Matrix::operator* (const Matrix& other)	// multiply this matrix and another
{
	assert((mCols == other.rows()) && "Matrix multiplication: dimension mismatch.");
	Matrix c(mRows, other.cols());
	for (unsigned i = 0; i < mRows; i++)
		for (unsigned j = 0; j < other.cols(); j++)
		{
			double* s = &c[i][j]; *s = 0;
			for (unsigned k = 0; k < mCols; k++)
				*s += (*this).mData[i * mCols + k] * other.mData[k * other.mCols + j];
		}
	return c;
}


inline Matrix Matrix::operator/ (Matrix b)	// Solve A* x= b
{
	return SolveQRQ(decomposeQRQ(*this), b);
}


inline Vector<double> Matrix::operator/ (Vector<double> b)	// find A/ b
{
	return Vector<double>(SolveQRQ(decomposeQRQ(*this), Matrix(b)));
}


inline double Frobenius(const Matrix& mat)	// weighted Frobenius norm
{
	double* tmp = mat[0] + mat.rows() * mat.cols();
	double tmpNorm = 0;				// Hilbert-Schmidt
	for (double* s = mat[0]; s < tmp; s++)
		tmpNorm += Square(*s);
	return pow(tmpNorm / mat.rows() / mat.cols(), 1 / 2);
}


inline std::ostream& operator << (std::ostream& os, const Matrix m)	// output
{
	os << "(" << m.rows() << "x" << m.cols() << " matrix):";
	for (unsigned i = 0; i < m.rows(); i++)
	{
		os << std::endl;		// new row
		for (unsigned j = 0; j < m.cols(); j++)
			os << std::setw(11) << std::right << m[i][j] << " ";
	}
	return os;
}


// void Matrix::input()			//	manual input of matrix via console
// {	std::cout << "Enter the elements of the " << mRows << " x " << mCols << " Matrix: " << std::endl;
// 	for(unsigned i=0; i< mRows ; i++)
// 	{	std::cout << "row " << i+1 << ": ";
// 		for(unsigned j= 0; j< mCols; j++)
// 			std::cin >> mData[i* mCols+ j];}}


// Matrix input()					//	manual input of matrix via console
// {	unsigned	mRows, mCols;
// 	std::cout << std::endl << "Matrix input. Enter the dimension of Matrix: "; std::cin >> mRows >> mCols;
// 	Matrix c(mRows, mCols); c.input(); return c;}



#pragma GCC diagnostic pop
#endif		// matrix_H
