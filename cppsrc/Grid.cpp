/*
Grid.cpp
generated grid
Drew Murray
11/29/16
*/


//libraries
#include <vector>
using std::vector;
#include <sstream>
using std::stringstream;
using std::istringstream;
#include <iostream>
using std::cout;
using std::cerr;
#include <string>
using std::string;
using std::to_string;
using std::stoi;
#include <fstream>
using std::ifstream;
#include <bitset>
using std::bitset;
#include <limits>
using std::numeric_limits;
#include <cmath>
using std::pow;
#include <cstring>
#include <mpi.h>

#define SSTR( x ) ( std::ostringstream() << std::dec << x ).str()

//class includes
#include "Grid.h"
#include "Node.h"
#include "Block.h"
#include "BCParam.h"


#include <fstream>
using std::ofstream;

//forward declarations (only pointers used)
class Basis;

//taken from http://stackoverflow.com/questions/8616573/convert-64-bit-binary-string-representation-of-a-double-number-back-to-double-nu
double bitstring_to_double(const char* p)
{
	unsigned long long x = 0;
	for (; *p; ++p)
	{
		x = (x << 1) + (*p - '0');
	}
	double d;
	memcpy(&d, &x, 8);
	return d;
#include <mpi.h>
}

//shockingly, an existing implementation doesn't seem to exist in C++
double hex_to_double(const string& hexstr)
{
	//convert to binary
	stringstream ss;
	ss << std::hex << hexstr;
	unsigned long long n;
	ss >> n;
	bitset<64> b(n);
	bitset<1> s(b.to_string().substr(0, 1));
	bitset<11> e(b.to_string().substr(1, 11));
	bitset<52> m(b.to_string().substr(12, 52));
	return bitstring_to_double(b.to_string().c_str());
}

//constructors BE SURE TO UPDATE ALL OF THEM IF MEMBER DATA IS CHANGED
Grid::Grid(const Grid& g) : mBlock(0)//dummy block constructor
{
	int size = g.mNumBlocks[0] * g.mNumBlocks[1] * g.mNumBlocks[2];
	if (size==1)mBlocks=g.mBlocks;
	else mBlock=g.mBlock;
	mStiff = g.mStiff;
	mMass = g.mMass;
	mTime=g.mTime;
	mOrder=g.mOrder;
	mRank=g.mRank;
	mEvalAdvect=g.mEvalAdvect;
	mEvalDiff=g.mEvalDiff;
	mEvalSponge=g.mEvalSponge;
	std::copy(g.mNumElements,g.mNumElements+3,mNumElements);
	std::copy(g.mNumBlocks,g.mNumBlocks+3,mNumBlocks);
	BCParamPtr=g.BCParamPtr;
	MPI_Comm_rank (MPI_COMM_WORLD, &mRank); // get current MPI-process ID. O, 1, ...
	// initialize neighbors in blocks

	if (size==1)
	{
		for(unsigned int i=0;i<mBlocks.size();i++)
		{
			for(unsigned int j=0;j<mBlocks[0].size();j++)
			{
				for(unsigned int k=0;k<mBlocks[0][0].size();k++)
				{
					if(i>0)
					{
						mBlocks[i][j][k].initializeNeighbor(4,mBlocks[i-1][j][k].getRank());
					}
					if(i<mBlocks.size()-1)
					{
						mBlocks[i][j][k].initializeNeighbor(5, mBlocks[i+1][j][k].getRank());
					}
					if(j>0)
					{
						mBlocks[i][j][k].initializeNeighbor(2,mBlocks[i][j-1][k].getRank());
					}
					if(j<mBlocks[0].size()-1)
					{
						mBlocks[i][j][k].initializeNeighbor(3,mBlocks[i][j+1][k].getRank());
					}
					if(k>0)
					{
						mBlocks[i][j][k].initializeNeighbor(0,mBlocks[i][j][k-1].getRank());
					}
					if(k<mBlocks[0][0].size()-1)
					{
						mBlocks[i][j][k].initializeNeighbor(1,mBlocks[i][j][k+1].getRank());
					}
				}
			}
		}
		if (BCParamPtr->IS_PERIODIC)
		{
			for(unsigned int i=0;i<mBlocks.size();i++)
			{
				for(unsigned int j=0;j<mBlocks[0].size();j++)
				{
					for(unsigned int k=0;k<mBlocks[0][0].size();k++)
					{
						if(i==0)
						{
							mBlocks[i][j][k].initializeNeighbor(4,mBlocks[mBlocks.size()-1][j][k].getRank());
						}
						if(i==mBlocks.size()-1)
						{
							mBlocks[i][j][k].initializeNeighbor(5, mBlocks[0][j][k].getRank());
						}
						if(j==0)
						{
							mBlocks[i][j][k].initializeNeighbor(2,mBlocks[i][mBlocks[0].size()-1][k].getRank());
						}
						if(j==mBlocks[0].size()-1)
						{
							mBlocks[i][j][k].initializeNeighbor(3,mBlocks[i][0][k].getRank());
						}
						if(k==0)
						{
							mBlocks[i][j][k].initializeNeighbor(0,mBlocks[i][j][mBlocks[0][0].size()-1].getRank());
						}
						if(k==mBlocks[0][0].size()-1)
						{
							mBlocks[i][j][k].initializeNeighbor(1,mBlocks[i][j][0].getRank());
						}
					}
				}
			}
		}
	}
	else
	{
		for (int i=0; i<6; i++) if(g.mBlock.getNeighborExists(i)) mBlock.initializeNeighbor(i,g.mBlock.getNeighborRank(i));
	}

}

Grid::Grid(string matrixFilename,BCParam* BCPtr, unsigned short int order, unsigned int* numBlocks, unsigned int* numElem, float* elemWidth, float* startingCorner, Basis* basisPtr, float Rho0, float C0, float BoverA, float diffusionConstant,bool evalAdvect,bool evalDiff,bool evalSponge) : mBlock(0)//dummy block constructor
{
	int size = numBlocks[0] * numBlocks[1] * numBlocks[2];
	MPI_Comm_rank (MPI_COMM_WORLD, &mRank); // get current MPI-process ID. O, 1, ...
	mTime=0;
	BCParamPtr=BCPtr;
	std::copy(numElem,numElem+3,mNumElements);
	std::copy(numBlocks,numBlocks+3,mNumBlocks);
	mOrder=order;
	mEvalAdvect=evalAdvect;
	mEvalDiff=evalDiff;
	mEvalSponge=evalSponge;
	readMatrices(matrixFilename, mOrder);
	int globalPosition[3]={0,0,0};
	float blockCorner[3]={0,0,0};
	//float coord[3]={0.,0.,0.};
	int num=0;
	if (size==1)
	{
		for(unsigned int i=0;i<numBlocks[0];i++)
		{
			vector<vector<Block>> plane;
			//coord[0]=i*elemWidth[0]+startingCorner[0];
			globalPosition[0] = i;
			blockCorner[0]=startingCorner[0]+elemWidth[0]*numElem[0]*i;
			for(unsigned int j=0;j<numBlocks[1];j++)
			{
				vector<Block> line;
				//coord[1]=j*elemWidth[1]+startingCorner[1];
				globalPosition[1] = j;
				blockCorner[1]=startingCorner[1]+elemWidth[1]*numElem[1]*j;
				for(unsigned int k=0;k<numBlocks[2];k++)
				{
					//coord[2]=k*elemWidth[2]+startingCorner[2];
					globalPosition[2] = k;
					blockCorner[2]=startingCorner[2]+elemWidth[2]*numElem[2]*k;
					line.push_back(Block(BCPtr, order, globalPosition, numBlocks, num, numElem, elemWidth, blockCorner, basisPtr, Rho0, C0, BoverA, diffusionConstant,getMass(), getStiff()));
					num++;
				}
				plane.push_back(line);
			}
			mBlocks.push_back(plane);
		}
	
		// initialize neighbors in blocks
		for(unsigned int i=0;i<mBlocks.size();i++)
		{
			for(unsigned int j=0;j<mBlocks[0].size();j++)
			{
				for(unsigned int k=0;k<mBlocks[0][0].size();k++)
				{
					if(i>0)
					{
						mBlocks[i][j][k].initializeNeighbor(4,mBlocks[i-1][j][k].getRank());
					}
					if(i<mBlocks.size()-1)
					{
						mBlocks[i][j][k].initializeNeighbor(5, mBlocks[i+1][j][k].getRank());
					}
					if(j>0)
					{
						mBlocks[i][j][k].initializeNeighbor(2,mBlocks[i][j-1][k].getRank());
					}
					if(j<mBlocks[0].size()-1)
					{
						mBlocks[i][j][k].initializeNeighbor(3,mBlocks[i][j+1][k].getRank());
					}
					if(k>0)
					{
						mBlocks[i][j][k].initializeNeighbor(0,mBlocks[i][j][k-1].getRank());
					}
					if(k<mBlocks[0][0].size()-1)
					{
						mBlocks[i][j][k].initializeNeighbor(1,mBlocks[i][j][k+1].getRank());
					}
				}
			}
		}
		for(unsigned int i=0;i<mBlocks.size();i++)
		{
			for(unsigned int j=0;j<mBlocks[0].size();j++)
			{
				for(unsigned int k=0;k<mBlocks[0][0].size();k++)
				{
					if(i==0 && BCPtr->IS_PERIODIC[0])
					{
						mBlocks[i][j][k].initializeNeighbor(4,mBlocks[mBlocks.size()-1][j][k].getRank());
					}
					if(i==mBlocks.size()-1 && BCPtr->IS_PERIODIC[0])
					{
						mBlocks[i][j][k].initializeNeighbor(5, mBlocks[0][j][k].getRank());
					}
					if(j==0 && BCPtr->IS_PERIODIC[1])
					{
						mBlocks[i][j][k].initializeNeighbor(2,mBlocks[i][mBlocks[0].size()-1][k].getRank());
					}
					if(j==mBlocks[0].size()-1 && BCPtr->IS_PERIODIC[1])
					{
						mBlocks[i][j][k].initializeNeighbor(3,mBlocks[i][0][k].getRank());
					}
					if(k==0 && BCPtr->IS_PERIODIC[2])
					{
						mBlocks[i][j][k].initializeNeighbor(0,mBlocks[i][j][mBlocks[0][0].size()-1].getRank());
					}
					if(k==mBlocks[0][0].size()-1 && BCPtr->IS_PERIODIC[2])
					{
						mBlocks[i][j][k].initializeNeighbor(1,mBlocks[i][j][0].getRank());
					}
				}
			}
		}
	}
	else//size>1
	{
		vector<vector<vector<int>>> dummyBlockNums;
		unsigned int i2=0, j2=0, k2=0;
		for(unsigned int i=0;i<numBlocks[0];i++)
		{
			vector<vector<int>> plane;
			//coord[0]=i*elemWidth[0]+startingCorner[0];
			globalPosition[0] = i;
			blockCorner[0]=startingCorner[0]+elemWidth[0]*numElem[0]*i;
			for(unsigned int j=0;j<numBlocks[1];j++)
			{
				vector<int> line;
				//coord[1]=j*elemWidth[1]+startingCorner[1];
				globalPosition[1] = j;
				blockCorner[1]=startingCorner[1]+elemWidth[1]*numElem[1]*j;
				for(unsigned int k=0;k<numBlocks[2];k++)
				{
					//coord[2]=k*elemWidth[2]+startingCorner[2];
					globalPosition[2] = k;
					blockCorner[2]=startingCorner[2]+elemWidth[2]*numElem[2]*k;
					if (num==mRank)
					{
						mBlock=Block(BCPtr, order, globalPosition, numBlocks, num, numElem, elemWidth, blockCorner, basisPtr, Rho0, C0, BoverA, diffusionConstant, getMass(),getStiff());
						i2=i;j2=j;k2=k;
					}
					line.push_back(num);
					num++;
				}
				plane.push_back(line);
			}
			dummyBlockNums.push_back(plane);
		}
		// initialize neighbors in blocks
		if(i2>0)
		{
			mBlock.initializeNeighbor(4,dummyBlockNums[i2-1][j2][k2]);
		}
		if(i2<dummyBlockNums.size()-1)
		{
			mBlock.initializeNeighbor(5, dummyBlockNums[i2+1][j2][k2]);
		}
		if(j2>0)
		{
			mBlock.initializeNeighbor(2,dummyBlockNums[i2][j2-1][k2]);
		}
		if(j2<dummyBlockNums[0].size()-1)
		{
			mBlock.initializeNeighbor(3,dummyBlockNums[i2][j2+1][k2]);
		}
		if(k2>0)
		{
			mBlock.initializeNeighbor(0,dummyBlockNums[i2][j2][k2-1]);
		}
		if(k2<dummyBlockNums[0][0].size()-1)
		{
			mBlock.initializeNeighbor(1,dummyBlockNums[i2][j2][k2+1]);
		}
		if(i2==0 && BCPtr->IS_PERIODIC[0])
		{
			mBlock.initializeNeighbor(4,dummyBlockNums[dummyBlockNums.size()-1][j2][k2]);
		}
		if(i2==dummyBlockNums.size()-1 && BCPtr->IS_PERIODIC[0])
		{
			mBlock.initializeNeighbor(5, dummyBlockNums[0][j2][k2]);
		}
		if(j2==0 && BCPtr->IS_PERIODIC[1])
		{
			mBlock.initializeNeighbor(2,dummyBlockNums[i2][dummyBlockNums[0].size()-1][k2]);
		}
		if(j2==dummyBlockNums[0].size()-1 && BCPtr->IS_PERIODIC[1])
		{
			mBlock.initializeNeighbor(3,dummyBlockNums[i2][0][k2]);
		}
		if(k2==0 && BCPtr->IS_PERIODIC[2])
		{
			mBlock.initializeNeighbor(0,dummyBlockNums[i2][j2][dummyBlockNums[0][0].size()-1]);
		}
		if(k2==dummyBlockNums[0][0].size()-1 && BCPtr->IS_PERIODIC[2])
		{
			mBlock.initializeNeighbor(1,dummyBlockNums[i2][j2][0]);
		}
	}

		
		

}

//methods
vector<vector<vector<Block>>> Grid::getBlocksCopy()
{
	vector<vector<vector<Block>>> ret = mBlocks;
	return ret;
}

//For debugging purposes, implements periodic BC without using MPI
void Grid::periodicSelf()
{
	int fromi = -1;
	int fromj = -1;
	int fromk = -1;
	bool send = false;
	for (int i = 0;i < mNumElements[0] + 2;i++)
	{
		for (int j = 0;j < mNumElements[1] + 2;j++)
		{
			for (int k = 0;k < mNumElements[2] + 2;k++)
			{
				fromi = i;
				fromj = j;
				fromk = k;
				send = false;
				if (i == 0)
				{
					fromi = mNumElements[0];
					if (BCParamPtr->IS_PERIODIC[0]) send = true;
				}
				else if (i == mNumElements[0] + 1)
				{
					fromi = 1;
					if (BCParamPtr->IS_PERIODIC[0]) send = true;
				}
				else if (j == 0)
				{
					fromj = mNumElements[1];
					if (BCParamPtr->IS_PERIODIC[1]) send = true;
				}
				else if (j == mNumElements[1] + 1)
				{
					fromj = 1;
					if (BCParamPtr->IS_PERIODIC[1]) send = true;
				}
				else if (k == 0)
				{
					fromk = mNumElements[2];
					if (BCParamPtr->IS_PERIODIC[2]) send = true;
				}
				else if (k == mNumElements[2] + 1)
				{
					fromk = 1;
					if (BCParamPtr->IS_PERIODIC[2]) send = true;
				}
				if (send)
				{
					mBlocks[0][0][0].mElements[i][j][k].setNodeData(&mBlocks[0][0][0].mElements[fromi][fromj][fromk]);
					mBlocks[0][0][0].mElements[i][j][k].setPhysicalConstants(&mBlocks[0][0][0].mElements[fromi][fromj][fromk]);
				}

			}
		}
	}
}

//For debugging purposes: get all DRho values in the entire grid into a string.  twoDOnly will make it ignore all layers in Z other than the first one
string Grid::allDRhoToString(bool twoDOnly)
{
	string ret;
	for (unsigned int i = 0;i < mBlocks.size();i++)
	{
		for (unsigned int j = 0;j < mBlocks[i].size();j++)
		{
			for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
			{
				if (!twoDOnly || k == 0)
				{
					ret = ret + "Block " + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "\n" + mBlocks[i][j][k].allDRhoToString(twoDOnly);
				}
			}
		}
	}
	return ret;
}

//Set state variables to zero in entire grid
void Grid::setZero()
{

	int size = mNumBlocks[0] * mNumBlocks[1] * mNumBlocks[2];
	if (size == 1)
	{
		for (unsigned int i = 0;i < mBlocks.size();i++)
		{
			for (unsigned int j = 0;j < mBlocks[i].size();j++)
			{
				for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
				{
					mBlocks[i][j][k].setZero();
				}
			}
		}
	}
	else//size>1
	{
		mBlock.setZero();
	}
}

//Calls all communication to sync state variables among processes
void Grid::communicate()
{
	int size = mNumBlocks[0] * mNumBlocks[1] * mNumBlocks[2];
	if (size == 1)
	{
		periodicSelf();


	}
	else//size>1
	{
		

		MPI_Barrier(MPI_COMM_WORLD);
		vector<MPI_Request> mpi_requests;
		mBlock.resetHaloBool();
		unsigned int myIdx=0;
		unsigned int numBlocks=0;
		bool skipPhase2=false;
		for (unsigned int dir = 0;dir< 6;dir++)//direction
		{
			skipPhase2=false;
			if (dir==0 || dir==1)
			{
				myIdx=mBlock.getPosIdx()[2];
				numBlocks=mNumBlocks[2];
			}
			else if (dir==2 || dir==3) 
			{
				myIdx=mBlock.getPosIdx()[1];
				numBlocks=mNumBlocks[1];
			}
			else //dir is 4 or 5
			{
				myIdx=mBlock.getPosIdx()[0];
				numBlocks=mNumBlocks[0];
			}
			if (myIdx%2==0)
			{
				mBlock.sendFaceInfo(dir,mpi_requests);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (myIdx%2==1)
			{
				mBlock.receiveFaceInfo(dir,mpi_requests);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			//if odd number of blocks AND PERIODIC, then the zeroth block will be 'even' and receive from an 'even' instead of an 'odd'
			if (numBlocks%2==1)
			{
				if(myIdx==numBlocks-1 && mBlock.getNeighborExists(dir) && dir%2==0)//dir being even means we are sending in the negative direction
				{
					mBlock.receiveFaceInfo(dir,mpi_requests);
					skipPhase2=true;
				}
				//the "else" is important for the numBlocks=1 case, in which case we only want to recv once (doesn't matter from which side)
				else if(myIdx==0 && mBlock.getNeighborExists(dir) && dir%2==1)//dir being odd means we are sending in the positive direction
				{
					mBlock.receiveFaceInfo(dir,mpi_requests);
					skipPhase2=true;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);

			//then the other half of processes send/receive
			if (myIdx%2==1)
			{
				mBlock.sendFaceInfo(dir,mpi_requests);
			}
			MPI_Barrier(MPI_COMM_WORLD);
			if (myIdx%2==0 and !skipPhase2)
			{
				mBlock.receiveFaceInfo(dir,mpi_requests);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}

void Grid::forwardEuler(float dt, Grid* oldGrid, ofstream& logFile)
{
	int size = mNumBlocks[0] * mNumBlocks[1] * mNumBlocks[2];
	mTime += dt;
	oldGrid->mTime += dt;

	setZero();
	oldGrid->communicate();


	if (size == 1)
	{
		for (unsigned int i = 0;i < mBlocks.size();i++)
		{
			for (unsigned int j = 0;j < mBlocks[i].size();j++)
			{
				for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
				{
					if (mBlocks[i][j][k].getRank() == mRank)
					{
						oldGrid->mBlocks[i][j][k].OperateQFromDRho(0, 0, &(oldGrid->mBlocks[i][j][k]));
					}
				}
			}
		}
		oldGrid->communicate();
		for (unsigned int i = 0;i < mBlocks.size();i++)
		{
			for (unsigned int j = 0;j < mBlocks[i].size();j++)
			{
				for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
				{
					if (mBlocks[i][j][k].getRank() == mRank)
					{
						mBlocks[i][j][k].OperateDRhoFromQ(0, 0, &(oldGrid->mBlocks[i][j][k]));
					}
				}
			}
		}
		for (unsigned int i = 0;i < mBlocks.size();i++)
		{
			for (unsigned int j = 0;j < mBlocks[i].size();j++)
			{
				for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
				{
					if (mBlocks[i][j][k].getRank() == mRank)
					{
						for (unsigned int axis = 0; axis < 3; axis++)
						{
							mBlocks[i][j][k].OperateUFromV(axis, 0, &(oldGrid->mBlocks[i][j][k]));
							mBlocks[i][j][k].OperateUFromV(axis, 1, &(oldGrid->mBlocks[i][j][k]));
						}
						mBlocks[i][j][k].OperateSponge(0, 0, &(oldGrid->mBlocks[i][j][k]));
						mBlocks[i][j][k].selfWeightedAdd(&(oldGrid->mBlocks[i][j][k]),dt, 1.0);

					}
				}
			}
		}
	}
	else //size >1
	{
		oldGrid->mBlock.OperateQFromDRho(0, 0, &(oldGrid->mBlock));
		oldGrid->communicate();
		mBlock.OperateDRhoFromQ(0, 0, &(oldGrid->mBlock));
		for (unsigned int axis = 0; axis < 3; axis++)
		{
			mBlock.OperateUFromV(axis, 0, &(oldGrid->mBlock));
			mBlock.OperateUFromV(axis, 1, &(oldGrid->mBlock));
		}
		mBlock.OperateSponge(0, 0, &(oldGrid->mBlock));
		mBlock.selfWeightedAdd(&(oldGrid->mBlock), dt, 1.0);
	}
}

void Grid::RK33(float dt, Grid* oldGrid, Grid* k1, Grid* k2, ofstream& logFile)
{

	int size = mNumBlocks[0]* mNumBlocks[1]* mNumBlocks[2];
	mTime += dt;
	oldGrid->mTime += dt;
	k1->setZero();
	k2->setZero();
	setZero();
	if (size == 1)
	{

		oldGrid->communicate();
		if (mEvalDiff)
		{
			for (unsigned int i = 0;i < mBlocks.size();i++)
			{
				for (unsigned int j = 0;j < mBlocks[i].size();j++)
				{
					for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
					{
						if (mBlocks[i][j][k].getRank() == mRank)
						{
							for (unsigned int axis = 0; axis < 3; axis++)
							{
								oldGrid->mBlocks[i][j][k].OperateQFromDRho(axis, 0, &(oldGrid->mBlocks[i][j][k]));
							}
						}
					}
				}
			}
			oldGrid->communicate();
			for (unsigned int i = 0;i < mBlocks.size();i++)
			{
				for (unsigned int j = 0;j < mBlocks[i].size();j++)
				{
					for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
					{
						if (mBlocks[i][j][k].getRank() == mRank)
						{
							for (unsigned int axis = 0; axis < 3; axis++)
							{
								k1->mBlocks[i][j][k].OperateDRhoFromQ(axis, 0, &(oldGrid->mBlocks[i][j][k]));
							}
						}
					}
				}
			}
		}
		for (unsigned int i = 0;i < mBlocks.size();i++)
		{
			for (unsigned int j = 0;j < mBlocks[i].size();j++)
			{
				for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
				{
					if (mBlocks[i][j][k].getRank() == mRank)
					{
						if(mEvalAdvect) for (unsigned int axis = 0; axis < 3; axis++)
						{
							k1->mBlocks[i][j][k].OperateUFromV(axis, 0, &(oldGrid->mBlocks[i][j][k]));//k1 = k1 + oper
							k1->mBlocks[i][j][k].OperateUFromV(axis, 1, &(oldGrid->mBlocks[i][j][k]));
						}
						if(mEvalSponge)k1->mBlocks[i][j][k].OperateSponge(0, 0, &(oldGrid->mBlocks[i][j][k]));
						k1->mBlocks[i][j][k].selfWeightedAdd(&(oldGrid->mBlocks[i][j][k]), dt, 1.0);//k1 = k1*dt + old*1


					}
				}
			}
		}

		k1->communicate();
		if (mEvalDiff)
		{
			for (unsigned int i = 0;i < mBlocks.size();i++)
			{
				for (unsigned int j = 0;j < mBlocks[i].size();j++)
				{
					for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
					{
						if (mBlocks[i][j][k].getRank() == mRank)
						{
							for (unsigned int axis = 0; axis < 3; axis++)
							{
								k1->mBlocks[i][j][k].OperateQFromDRho(axis, 0, &(k1->mBlocks[i][j][k]));
							}
						}
					}
				}
			}
			k1->communicate();
			for (unsigned int i = 0;i < mBlocks.size();i++)
			{
				for (unsigned int j = 0;j < mBlocks[i].size();j++)
				{
					for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
					{
						if (mBlocks[i][j][k].getRank() == mRank)
						{
							for (unsigned int axis = 0; axis < 3; axis++)
							{
								k2->mBlocks[i][j][k].OperateDRhoFromQ(axis, 0, &(k1->mBlocks[i][j][k]));
							}
						}
					}
				}
			}
		}

		for (unsigned int i = 0;i < mBlocks.size();i++)
		{
			for (unsigned int j = 0;j < mBlocks[i].size();j++)
			{
				for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
				{
					if (mBlocks[i][j][k].getRank() == mRank)
					{
						if (mEvalAdvect) for (unsigned int axis = 0; axis < 3; axis++)
						{
							k2->mBlocks[i][j][k].OperateUFromV(axis, 0, &(k1->mBlocks[i][j][k]));//k2 = k2 + oper
							k2->mBlocks[i][j][k].OperateUFromV(axis, 1, &(k1->mBlocks[i][j][k]));
						}
						if (mEvalSponge)k2->mBlocks[i][j][k].OperateSponge(0, 0, &(k1->mBlocks[i][j][k]));
						k2->mBlocks[i][j][k].selfWeightedAdd(&(k1->mBlocks[i][j][k]),dt*.25F,.25F);//k2 = k2*dt/4 + k1/4
						k2->mBlocks[i][j][k].selfWeightedAdd(&(oldGrid->mBlocks[i][j][k]),1.0F,.75F);//k2 = k2*1 + old*3/4

					}
				}
			}
		}
		
		k2->communicate();

		if (mEvalDiff)
		{
			for (unsigned int i = 0;i < mBlocks.size();i++)
			{
				for (unsigned int j = 0;j < mBlocks[i].size();j++)
				{
					for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
					{
						if (mBlocks[i][j][k].getRank() == mRank)
						{
							for (unsigned int axis = 0; axis < 3; axis++)
							{
								k2->mBlocks[i][j][k].OperateQFromDRho(axis, 0, &(k2->mBlocks[i][j][k]));
							}
						}
					}
				}
			}
			k2->communicate();
			for (unsigned int i = 0;i < mBlocks.size();i++)
			{
				for (unsigned int j = 0;j < mBlocks[i].size();j++)
				{
					for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
					{
						if (mBlocks[i][j][k].getRank() == mRank)
						{
							for (unsigned int axis = 0; axis < 3; axis++)
							{
								mBlocks[i][j][k].OperateDRhoFromQ(axis, 0, &(k2->mBlocks[i][j][k]));
							}
						}
					}
				}
			}
		}

		for (unsigned int i = 0;i < mBlocks.size();i++)
		{
			for (unsigned int j = 0;j < mBlocks[i].size();j++)
			{
				for (unsigned int k = 0;k < mBlocks[i][j].size();k++)
				{
					if (mBlocks[i][j][k].getRank() == mRank)
					{
						if (mEvalAdvect) for (unsigned int axis = 0; axis < 3; axis++)
						{
							mBlocks[i][j][k].OperateUFromV(axis, 0, &(k2->mBlocks[i][j][k]));//self = self + oper
							mBlocks[i][j][k].OperateUFromV(axis, 1, &(k2->mBlocks[i][j][k]));
						}
						if (mEvalSponge)mBlocks[i][j][k].OperateSponge(0, 0, &(k2->mBlocks[i][j][k]));
						mBlocks[i][j][k].selfWeightedAdd(&(k2->mBlocks[i][j][k]), dt * 2.0F/3.0F, 2.0F/3.0F);//self = self*dt*2/3 + k2*2/3
						mBlocks[i][j][k].selfWeightedAdd(&(oldGrid->mBlocks[i][j][k]), 1.0F, 1.0F/3.0F);//self = self*1 + old*1/3

					}
				}
			}
		}
		

	}
	else // size>1
	{
		oldGrid->communicate();

		if (mEvalDiff)
		{
			for (unsigned int axis = 0; axis < 3; axis++)
			{
				oldGrid->mBlock.OperateQFromDRho(axis, 0, &(oldGrid->mBlock));
			}
			oldGrid->communicate();
			for (unsigned int axis = 0; axis < 3; axis++)
			{
				k1->mBlock.OperateDRhoFromQ(axis, 0, &(oldGrid->mBlock));
			}
		}
		if (mEvalAdvect) for (unsigned int axis = 0; axis < 3; axis++)
		{
			k1->mBlock.OperateUFromV(axis, 0, &(oldGrid->mBlock));//k1 = k1 + oper
			k1->mBlock.OperateUFromV(axis, 1, &(oldGrid->mBlock));
		}
		if (mEvalSponge)k1->mBlock.OperateSponge(0, 0, &(oldGrid->mBlock));
		k1->mBlock.selfWeightedAdd(&(oldGrid->mBlock), dt, 1.0);//k1 = k1*dt + old*1
		k1->communicate();
		if (mEvalDiff)
		{
			for (unsigned int axis = 0; axis < 3; axis++)
			{
				k1->mBlock.OperateQFromDRho(axis, 0, &(k1->mBlock));
			}
			k1->communicate();
			for (unsigned int axis = 0; axis < 3; axis++)
			{
				k2->mBlock.OperateDRhoFromQ(axis, 0, &(k1->mBlock));
			}
		}
		if (mEvalAdvect) for (unsigned int axis = 0; axis < 3; axis++)
		{
			k2->mBlock.OperateUFromV(axis, 0, &(k1->mBlock));//k2 = k2 + oper
			k2->mBlock.OperateUFromV(axis, 1, &(k1->mBlock));
		}
		if (mEvalSponge)k2->mBlock.OperateSponge(0, 0, &(k1->mBlock));
		k2->mBlock.selfWeightedAdd(&(k1->mBlock), dt * .25F, .25F);//k2 = k2*dt/4 + k1/4
		k2->mBlock.selfWeightedAdd(&(oldGrid->mBlock), 1.0F, .75F);//k2 = k2*1 + old*3/4
		k2->communicate();
		if (mEvalDiff)
		{
			for (unsigned int axis = 0; axis < 3; axis++)
			{
				k2->mBlock.OperateQFromDRho(axis, 0, &(k2->mBlock));
			}
			k2->communicate();
			for (unsigned int axis = 0; axis < 3; axis++)
			{
				mBlock.OperateDRhoFromQ(axis, 0, &(k2->mBlock));
			}
		}
		if (mEvalAdvect) for (unsigned int axis = 0; axis < 3; axis++)
		{
			mBlock.OperateUFromV(axis, 0, &(k2->mBlock));//self = self + oper
			mBlock.OperateUFromV(axis, 1, &(k2->mBlock));
		}
		if (mEvalSponge)mBlock.OperateSponge(0, 0, &(k2->mBlock));
		mBlock.selfWeightedAdd(&(k2->mBlock), dt * 2.0F / 3.0F, 2.0F / 3.0F);//self = self*dt*2/3 + k2*2/3
		mBlock.selfWeightedAdd(&(oldGrid->mBlock), 1.0F, 1.0F / 3.0F);//self = self*1 + old*1/3
	}
}

string Grid::getLayout()
{
	string ret="";
	for(unsigned int i=0;i<mBlocks.size();i++)
	{
		for(unsigned int j=0;j<mBlocks[i].size();j++)
		{
			for(unsigned int k=0;k<mBlocks[i][j].size();k++)
			{
				ret+=mBlocks[i][j][k].getLayout()+"\n";
			}
		}
	}
	return ret;
}

void Grid::readMatrices(string filename, int targetOrder)
{
	if(mRank==0)cout << " - Creating Mass/Stiffness from matrices file\n";

	//Lines in the read-in file have a patterned order
	/*
		order
		M, M, ...
		M, M, ...
		S, S, ...
		S, S, ...
		order
		M, M, M, ...
		M, M, M, ...
		M, M, M, ...
		S, S, S, ...
		S, S, S, ...
		S, S, S, ...
		...
	*/
	ifstream file;
	string line;
	string tempString;
	unsigned char itr = 0;//which zero within an order
	unsigned char order = 0;
	int colNum = 0;

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
	mMass = Matrix(targetOrder + 1);
	mStiff = Matrix(targetOrder + 1);
	file.open(filename);
	if (file.is_open())
	{
		while (getline(file, line))
		{
			if (line.substr(0, 6) == "Order=")//next order
			{
				order = stoi(line.substr(6, line.length()-1));
			}
			else if (order == targetOrder)
			{
				istringstream iss(line);
				for (colNum = 0;colNum < order+1;colNum++)
				{
					iss >> tempString;
					if (itr < order + 1)mMass[itr][colNum] = hex_to_double(tempString);
					else mStiff[itr - (order + 1)][colNum] = hex_to_double(tempString);
				}
				itr += 1;
			}
			if (order > targetOrder) break;
		}
		file.close();
	}
	else cerr << "Unable to open file " << filename << "\n";
}

