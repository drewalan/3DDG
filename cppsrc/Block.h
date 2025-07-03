/*
Block.h
generated block of elements
Drew Murray
1/17/17
*/

//circular build protection
#ifndef BLOCK_H
#define BLOCK_H

#define MAXBLOCKMSGSIZE 28672000 // must be at least NUMELEM^3 * (Order+1)^3 * 7 : 28672000 is for 4th order, 32x32x32 elements

//libraries
#include <vector>
using std::vector;
#include <sstream>
#include <iostream>
#include <mpi.h>


//Some compilers want the top version of SSTR, some want the bottom (also must be swapped in Element.h)
//#define SSTR( x ) ( std::ostringstream() << std::dec << x ).str()
#define SSTR( x ) static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

//class includes
#include "Element.h"

//forward declarations (only pointers used)
class Basis;
struct BCParam;
struct TestCase;

class Block
{
	public:
		//constructors
		Block(BCParam* BCPtr, unsigned short int order, int* globalPosition, unsigned int* numBlocks, int globalAddress, unsigned int* numElem, float* elemWidth, float* startingCorner, Basis* basisPtr, float Rho0, float C0, float BoverA, float diffusionConstant, Matrix* massPtr, Matrix* stiffPtr);
		Block(int dummyVar){int x= dummyVar; int y=x; x=y;}//dummy constructor to bypass initialization list requirements and create empty instance
		Block(const Block& g);
		Block operator=(const Block& g);
		~Block(){if (msgBuffer) delete[] msgBuffer; msgBuffer = nullptr;} //destructor needed  to free heap memory



		//methods
		string allDRhoToString(bool twoDOnly);
		vector<vector<vector<Element>>> getElementsCopy();
		void OperateUFromV(unsigned short int axis, unsigned short int mode, Block* oldBlock);
		void OperateSponge(unsigned short int axis, unsigned short int mode, Block* oldBlock);
		void OperateQFromDRho(unsigned short int axis, unsigned short int mode, Block* oldBlock);
		void OperateDRhoFromQ(unsigned short int axis, unsigned short int mode, Block* oldBlock);

		void selfMult(float scalar);
		void selfWeightedAdd(Block* other,float selfScalar,float otherScalar);

		void resetHaloBool(){for(int i=0;i<6;i++) haloAlreadyTraded[i]=false;}

		void initializeNeighbor(int i, int r){neighborExists[i]=true;neighborRank[i]=r;}

		void applyBoundaryConditions(BCParam* parameters, float t, float dt);//implemented in BCParam.cpp
		void sendFaceInfo(unsigned int dir, vector<MPI_Request> requests);//uni-directional, checks for neighbor existance
		void receiveFaceInfo(unsigned int dir, vector<MPI_Request> requests);//uni-directional, checks for neighbor existance
		
		
		string getLayout();
		unsigned int* getPosIdx(){return mGlobalPosition;}
		string getName(){return "B"+SSTR(mGlobalAddress)+":"+SSTR(mGlobalPosition[0])+\
							","+SSTR(mGlobalPosition[1])+","+SSTR(mGlobalPosition[2]);}


		int getRank(){if (maxRank>1) return mGlobalAddress; else return 0;}
		int getMaxRank(){return maxRank;}
		int getNeighborRank(int dir) const {return neighborRank[dir];}
		bool getNeighborExists(int i) const {return neighborExists[i];}
		void setNeighborExists(bool val, int i){neighborExists[i]=val;}
		void setZero();


	//member variables
	private:

	public:
		vector<vector<vector<Element>>> mElements;// TODO make private
	private:

		bool neighborExists[6] = {false, false, false, false, false, false};
		int neighborRank[6] = {-1, -1, -1, -1, -1, -1};
		int maxRank = -1;
		int mOrder = -1;
		unsigned int mGlobalPosition[3] = {9999, 9999, 9999 };
		int mGlobalAddress = -1;
		unsigned int mNumElements[3] = { 0, 0, 0 };
		Matrix* massPtr = NULL;
		Matrix* stiffPtr = NULL;

		//order i-1,i+1,j-1,j+1,k-1,k+1
		bool haloAlreadyTraded[6] = { false, false, false, false, false, false };
		float* msgBuffer = new float[MAXBLOCKMSGSIZE];
		

	//disable default constructors
	private:
		Block();

};
#endif
