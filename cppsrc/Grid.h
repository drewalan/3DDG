/*
Grid.h
generated grid
Drew Murray
11/29/16
*/

//circular build protection
#ifndef GRID_H
#define GRID_H

//libraries
#include <vector>
using std::vector;

//class includes
#include "Block.h"
#include "Matrix.h"


#include <fstream>
using std::ofstream;

//forward declarations (only pointers used)
class Basis;
struct BCParam;
struct TestCase;

class Grid
{
	public:
		//constructors
		Grid(string filename,BCParam* BCP, unsigned short int order, unsigned int* numBlocks, unsigned int* numElem, float* elemWidth, float* startingCorner, Basis* basisPtr, float Rho0, float C0, float BoverA, float Mu_visc,bool evalAdvect,bool evalDiff,bool evalSponge);
		Grid(const Grid& g);

		//methods
		vector<vector<vector<Block>>> getBlocksCopy();
		void forwardEuler(float dt, Grid* oldGrid, ofstream& logFile);
		void RK33(float dt, Grid* oldGrid, Grid* k1, Grid* k2, ofstream& logFile);
		void communicate();
		string getLayout();
		void applyBoundaryConditions(float dt);//implemented in BCParam.cpp
		int* getNumBlocks(){return mNumBlocks;}
		int* getNumElements(){return mNumElements;}
		int getOrder(){return mOrder;}
		void periodicSelf();
		Matrix* getMass() { return &mMass; }
		Matrix* getStiff() { return &mStiff; }
		void readMatrices(string filename,int targetOrder);
		void setZero();
		string allDRhoToString(bool twoDOnly);



	//member variables
		vector<vector<vector<Block>>> mBlocks;
		Block mBlock;//for paralell, each "grid" has one block
	private:
		float mTime;
		int mRank;
		int mOrder;
		int mNumElements[3];
		int mNumBlocks[3];
		BCParam* BCParamPtr;
		Matrix mMass;
		Matrix mStiff;
		bool mEvalAdvect;
		bool mEvalDiff;
		bool mEvalSponge;

	//disable default constructors
	private:
		Grid();
		Grid operator=(const Grid& g);

};
#endif
