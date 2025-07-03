/*
Block.cpp
generated block of elements
Drew Murray
1/17/17
*/



//libraries
#include <vector>
using std::vector;
#include <iostream>
using std::cout;
#include <algorithm>
using std::min;


#define loopinternalnodes(ngl) for (int __i = 0; __i < ngl; __i++) for (int __j = 0; __j < ngl; __j++) for (int __k = 0; __k < ngl; __k++)

#define MAXBLOCKMSGSIZE 28672000 // must be at least NUMELEM^3 * (Order+1)^3 * 7 : 28672000 is for 4th order, 32x32x32 elements

//class includes
#include "Block.h"
#include "Element.h"
#include "Node.h"
#include "BCParam.h"

//forward declarations (only pointers used)
class Basis;
struct BCParam;

//constructors BE SURE TO UPDATE ALL OF THEM IF MEMBER DATA IS CHANGED
	Block::Block(const Block& g)
	{
		mElements=g.mElements;
		maxRank=g.maxRank;
		mOrder = g.mOrder;
		massPtr = g.massPtr;
		stiffPtr = g.stiffPtr;
		std::copy(g.neighborExists,g.neighborExists+6,neighborExists);
		std::copy(g.mGlobalPosition,g.mGlobalPosition+3,mGlobalPosition);
		std::copy(g.mNumElements,g.mNumElements+3,mNumElements);
		std::copy(g.haloAlreadyTraded,g.haloAlreadyTraded+6,haloAlreadyTraded);

		mGlobalAddress=g.mGlobalAddress;
		for(unsigned int i=0;i<mElements.size();i++)
		{
			for(unsigned int j=0;j<mElements[0].size();j++)
			{
				for(unsigned int k=0;k<mElements[0][0].size();k++)
				{
					mElements[i][j][k].setBlockPtr(this);
				}
			}
		}
		// initialize neighbors in elements
		for(unsigned int i=1;i<mElements.size()-1;i++)
		{
			for(unsigned int j=1;j<mElements[0].size()-1;j++)
			{
				for(unsigned int k=1;k<mElements[0][0].size()-1;k++)
				{
					if(i>0)
					{
						mElements[i][j][k].initializeNeighbor(&mElements[i-1][j][k],4);
					}
					if(i<g.mElements.size()-1)
					{
						mElements[i][j][k].initializeNeighbor(&mElements[i+1][j][k],5);
					}
					if(j>0)
					{
						mElements[i][j][k].initializeNeighbor(&mElements[i][j-1][k],2);
					}
					if(j<g.mElements[0].size()-1)
					{
						mElements[i][j][k].initializeNeighbor(&mElements[i][j+1][k],3);
					}
					if(k>0)
					{
						mElements[i][j][k].initializeNeighbor(&mElements[i][j][k-1],0);
					}
					if(k<g.mElements[0][0].size()-1)
					{
						mElements[i][j][k].initializeNeighbor(&mElements[i][j][k+1],1);
					}
				}
			}
		}
	}

	Block Block::operator=(const Block& g)
	{
		mElements=g.mElements;
		maxRank=g.maxRank;
		mOrder = g.mOrder;
		massPtr = g.massPtr;
		stiffPtr = g.stiffPtr;
		std::copy(g.neighborExists,g.neighborExists+6,neighborExists);
		std::copy(g.mGlobalPosition,g.mGlobalPosition+3,mGlobalPosition);
		std::copy(g.mNumElements,g.mNumElements+3,mNumElements);
		std::copy(g.haloAlreadyTraded,g.haloAlreadyTraded+6,haloAlreadyTraded);


		mGlobalAddress=g.mGlobalAddress;
		//initialize upward pointer (needs different loop bounds than initializing neighbors below)
		for(unsigned int i=0;i<mElements.size();i++)
		{
			for(unsigned int j=0;j<mElements[0].size();j++)
			{
				for(unsigned int k=0;k<mElements[0][0].size();k++)
				{
					mElements[i][j][k].setBlockPtr(this);
				}
			}
		}
		// initialize neighbors in elements
		for(unsigned int i=1;i<mElements.size()-1;i++)
		{
			for(unsigned int j=1;j<mElements[0].size()-1;j++)
			{
				for(unsigned int k=1;k<mElements[0][0].size()-1;k++)
				{
					// Follow JFK's order
					if(i>0) mElements[i][j][k].initializeNeighbor(&mElements[i-1][j][k],4);
					if(i<g.mElements.size()-1) mElements[i][j][k].initializeNeighbor(&mElements[i+1][j][k],5);
					if(j>0) mElements[i][j][k].initializeNeighbor(&mElements[i][j-1][k],2);
					if(j<g.mElements[0].size()-1) mElements[i][j][k].initializeNeighbor(&mElements[i][j+1][k],3);
					if(k>0) mElements[i][j][k].initializeNeighbor(&mElements[i][j][k-1],0);
					if(k<g.mElements[0][0].size()-1) mElements[i][j][k].initializeNeighbor(&mElements[i][j][k+1],1);
				}
			}
		}
		return *this;
	}


	Block::Block(BCParam* BCPtr, unsigned short int order, int* globalPosition, unsigned int* numBlocks, int globalAddress, unsigned int* numElem, float* elemWidth, float* startingCorner, Basis* basisPtr, float Rho0, float C0, float BoverA, float diffusionConstant, Matrix* massPointer, Matrix* stiffPointer)
	{
		std::copy(globalPosition,globalPosition+3,mGlobalPosition);
		std::copy(numElem,numElem+3,mNumElements);

		MPI_Comm_size (MPI_COMM_WORLD, &maxRank);


		mOrder = order;
		mGlobalAddress=globalAddress;
		massPtr = massPointer;
		stiffPtr = stiffPointer;
		for(int i =0; i<6;i++)
		{
			neighborExists[i]=false;
			neighborRank[i]=-1;
			haloAlreadyTraded[i]=false;
		}
		int elemPosition[3]={0,0,0};
		int elemGlobalPosition[3]={0,0,0};
		float coord[3]={0.,0.,0.};	
		float sigmaU = 0.F;
		float sigmaDRho = 0.F;
		bool isGridEdge[6]= {false,false,false,false,false,false};
		float distanceToEdge[6] = {1000000000,1000000000,1000000000,1000000000,1000000000,1000000000 };//default values will never be the closest edge
		int edgeAxis = -1;
		float shortestDistanceToEdge = -1.0F;
		unsigned int num=(numElem[0]+2)*(numElem[1]+2)*(numElem[2]+2)*globalAddress;
		for(int i=0;i<int(numElem[0]+2);i++)
		{
			isGridEdge[4]=false;
			isGridEdge[5]=false;
			vector<vector<Element>> plane;
			coord[0]=(i-1)*elemWidth[0]+startingCorner[0];
			elemPosition[0] = i-1;
			elemGlobalPosition[0] = mGlobalPosition[0]*numElem[0]+i-1;//mGlobal is for the BLOCK
			if (BCPtr->spongeWidth[0] > 0)distanceToEdge[0] = abs(elemGlobalPosition[0])*elemWidth[0];
			if (BCPtr->spongeWidth[1] > 0)distanceToEdge[1] = abs(int(numElem[0]*numBlocks[0]) - abs(elemGlobalPosition[0]) - 1)*elemWidth[0];
			if (i==1&&globalPosition[0]==0)
			{
				isGridEdge[4]=true;
			}
			if (i==int(numElem[0])&&globalPosition[0]+1==signed(numBlocks[0]))
			{
				isGridEdge[5]=true;
			}
			for(int j=0;j<int(numElem[1]+2);j++)
			{
				isGridEdge[2]=false;
				isGridEdge[3]=false;
				vector<Element> line;
				coord[1]=(j-1)*elemWidth[1]+startingCorner[1];
				elemPosition[1] = j-1;
				elemGlobalPosition[1] = mGlobalPosition[1]*numElem[1]+j-1;//mGlobal is for the BLOCK
				if (BCPtr->spongeWidth[2] > 0)distanceToEdge[2] = abs(elemGlobalPosition[1])*elemWidth[1];
				if (BCPtr->spongeWidth[3] > 0)distanceToEdge[3] = abs(int(numElem[1]*numBlocks[1]) - abs(elemGlobalPosition[1]) - 1)*elemWidth[1];

				if (j==1&&globalPosition[1]==0)
				{
					isGridEdge[2]=true;
				}
				if (j==int(numElem[1])&&globalPosition[1]+1==signed(numBlocks[1]))
				{
					isGridEdge[3]=true;
				}
				for(int k=0;k<int(numElem[2]+2);k++)
				{
					isGridEdge[0]=false;
					isGridEdge[1]=false;
					coord[2]=(k-1)*elemWidth[2]+startingCorner[2];
					elemPosition[2] = k-1;
					elemGlobalPosition[2] = mGlobalPosition[2]*numElem[2]+k-1;//mGlobal is for the BLOCK
					if(BCPtr->spongeWidth[4]>0)distanceToEdge[4] = abs(elemGlobalPosition[2])*elemWidth[2];
					if(BCPtr->spongeWidth[5]>0)distanceToEdge[5] = abs(int(numElem[2]*numBlocks[2]) - abs(elemGlobalPosition[2]) - 1)*elemWidth[2];

					//find which direction is the closest edge to measure depth inside the sponge layer
					edgeAxis=0;
					shortestDistanceToEdge=distanceToEdge[0];
					for(int testAxis=1;testAxis<6;testAxis++)
					{
						if (distanceToEdge[testAxis] < shortestDistanceToEdge)
						{
							edgeAxis = testAxis;
							shortestDistanceToEdge=distanceToEdge[testAxis];
						}
					}
					//normalize shortestDistance to [0,1] as a fraction of the sponge layer's width
					shortestDistanceToEdge /= float(BCPtr->spongeWidth[edgeAxis])*elemWidth[int(edgeAxis*.5)];//edgeAxis must be converted from 'base 6' to 'base 3'

					//if shortest distance is to a 0-width sponge layer, then all layers are 0-width, since >0-width directions get set to shorter values than the default
					if (BCPtr->spongeWidth[edgeAxis] == 0)sigmaU = 0.F;
					else {
						sigmaU = fmax(0.0F,1.0F - shortestDistanceToEdge) * BCPtr->sigmaMax;
					}
					sigmaDRho = sigmaU;
					if (k==1&&globalPosition[2]==0)
					{
						isGridEdge[0]=true;
					}
					if (k==int(numElem[2])&&globalPosition[2]+1==signed(numBlocks[2]))
					{
						isGridEdge[1]=true;
					}
					line.push_back(Element(BCPtr, numElem, num, elemPosition, elemGlobalPosition, order, coord, elemWidth, basisPtr, Rho0, C0, BoverA, diffusionConstant, sigmaU, sigmaDRho, getRank(), this, isGridEdge, massPtr, stiffPtr));
					num++;
				}
				plane.push_back(line);
			}
			mElements.push_back(plane);
		}
		// initialize neighbors in elements
		for(unsigned int i=0;i<mElements.size();i++)
		{
			for(unsigned int j=0;j<mElements[0].size();j++)
			{
				for(unsigned int k=0;k<mElements[0][0].size();k++)
				{
					if (i>0) mElements[i][j][k].initializeNeighbor(&mElements[i - 1][j][k], 4);
					if (i<mElements.size() - 1) mElements[i][j][k].initializeNeighbor(&mElements[i + 1][j][k], 5);
					if (j>0) mElements[i][j][k].initializeNeighbor(&mElements[i][j - 1][k], 2);
					if (j<mElements[0].size() - 1) mElements[i][j][k].initializeNeighbor(&mElements[i][j + 1][k], 3);
					if (k>0) mElements[i][j][k].initializeNeighbor(&mElements[i][j][k - 1], 0);
					if (k<mElements[0][0].size() - 1) mElements[i][j][k].initializeNeighbor(&mElements[i][j][k + 1], 1);

				}
			}
		}

	}

//methods
	vector<vector<vector<Element>>> Block::getElementsCopy()
	{
		vector<vector<vector<Element>>> ret = mElements;
		return ret;
	}

	//get ALL Drho values in this block to a VERY LONG flattened string for debugging
	string Block::allDRhoToString(bool twoDOnly)
	{
		string ret;
		for (unsigned int i = 0;i < mElements.size();i++)
		{
			for (unsigned int j = 0;j < mElements[i].size();j++)
			{
				for (unsigned int k = 0;k < mElements[i][j].size();k++)
				{
					if (!twoDOnly || k==1)
					{
						ret = ret + "   Elem " + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + "\n" + mElements[i][j][k].allDRhoToString(twoDOnly);
					}
				}
			}
		}
		return ret;
	}

	//set all element state variables to zero
	void Block::setZero()
	{
		for (unsigned int i = 0;i < mElements.size();i++)
		{
			for (unsigned int j = 0;j < mElements[i].size();j++)
			{
				for (unsigned int k = 0;k < mElements[i][j].size();k++)
				{
					mElements[i][j][k].selfMult(0.0F);
				}
			}
		}
	}

	//multiply all element state variables by a scalar
	void Block::selfMult(float scalar)
	{
		for (unsigned int i = 1;i < mElements.size() - 1;i++)
		{
			for (unsigned int j = 1;j < mElements[i].size() - 1;j++)
			{
				for (unsigned int k = 1;k < mElements[i][j].size() - 1;k++)
				{
					mElements[i][j][k].selfMult(scalar);

				}
			}
		}
	}
	//applies formula: self = self * selfscalar + other * otherscalar where self
	void Block::selfWeightedAdd(Block* other, float selfScalar, float otherScalar)
	{
		for (unsigned int i = 1;i < mElements.size() - 1;i++)
		{
			for (unsigned int j = 1;j < mElements[i].size() - 1;j++)
			{
				for (unsigned int k = 1;k < mElements[i][j].size() - 1;k++)
				{
					mElements[i][j][k].selfWeightedAdd(&(other->mElements[i][j][k]),selfScalar,otherScalar);

				}
			}
		}
	}

	//All operators WRITE to themselves and READ from the oldBlock pointer.
	void Block::OperateUFromV(unsigned short int axis, unsigned short int mode, Block* oldBlock)
	{
		
		for (unsigned int i = 1;i < mNumElements[0] + 1;i++)
		{
			for (unsigned int j = 1;j < mNumElements[1] + 1;j++)
			{
				for (unsigned int k = 1;k < mNumElements[2] + 1;k++)
				{
					mElements[i][j][k].OperateUFromV(axis, mode, &(oldBlock->mElements[i][j][k]),false,false,4);
				}
			}
		}
	}
	//All operators WRITE to themselves and READ from the oldBlock pointer. EXCEPT QFromDRho which writes to Q as an intermediate variable (in oldBlock)
	void Block::OperateQFromDRho(unsigned short int axis, unsigned short int mode, Block* oldBlock)
	{

		for (unsigned int i = 1;i < mNumElements[0] + 1;i++)
		{
			for (unsigned int j = 1;j < mNumElements[1] + 1;j++)
			{
				for (unsigned int k = 1;k < mNumElements[2] + 1;k++)
				{
					mElements[i][j][k].OperateQFromDRho(axis, mode, &(oldBlock->mElements[i][j][k]));
				}
			}
		}
	}
	//All operators WRITE to themselves and READ from the oldBlock pointer.
	void Block::OperateDRhoFromQ(unsigned short int axis, unsigned short int mode, Block* oldBlock)
	{

		for (unsigned int i = 1;i < mNumElements[0] + 1;i++)
		{
			for (unsigned int j = 1;j < mNumElements[1] + 1;j++)
			{
				for (unsigned int k = 1;k < mNumElements[2] + 1;k++)
				{
					mElements[i][j][k].OperateDRhoFromQ(axis, mode, &(oldBlock->mElements[i][j][k]));
				}
			}
		}
	}

	//All operators WRITE to themselves and READ from the oldBlock pointer.
	void Block::OperateSponge(unsigned short int axis, unsigned short int mode, Block* oldBlock)
	{

		for (unsigned int i = 1;i < mNumElements[0] + 1;i++)
		{
			for (unsigned int j = 1;j < mNumElements[1] + 1;j++)
			{
				for (unsigned int k = 1;k < mNumElements[2] + 1;k++)
				{
					mElements[i][j][k].OperateSponge(axis, mode, &(oldBlock->mElements[i][j][k]));
				}
			}
		}
	}

	//return a long string showing the indices of the consituent elements
	string Block::getLayout()
	{

		string ret = getName()+"\n";

		for(signed int k=0;k<signed(mElements[0][0].size());k++)
		{
			ret+="   k="+SSTR(k-1)+"\n";
			for(unsigned int i=0;i<mElements.size();i++)
			{
				ret+="      ";
				for(unsigned int j=0;j<mElements[i].size();j++)
				{
					ret+=mElements[i][j][k].getName()+"   ";
				}
				ret+="\n";
			}
		}
		return ret;
	}

	void Block::sendFaceInfo(unsigned int dir, vector<MPI_Request> requests)//uni-directional, checks for neighbor existence first
	{
		//if(!neighborExists[dir])cout<<"ERROR!"<<getName()<<", dir="<<dir<<"\n";
		if(neighborExists[dir] && !haloAlreadyTraded[dir])
		{
			haloAlreadyTraded[dir]=true;

			int toRank = getNeighborRank(dir);
			//CALCULATE WHICH ARE ON FACE: there's some arithmatic magic to avoid ~2000 lines of if-statements and repeated code
			//sending so halos are not included
			int minI=1,minJ=1,minK=1;//minimum values on the face
			int maxI=mNumElements[0],maxJ=mNumElements[1],maxK=mNumElements[2];//max values on the face
			int useI=1,useJ=1,useK=1;//whether to use the corresponding values for the destination address
			int replaceI=0,replaceJ=0,replaceK=0;//what to use instead if the "use" is 0
			//sending so side furthest in dir (using JFK's order)
			if (dir==4){replaceI=maxI+1;maxI=minI;useI=0;}//-i			
			if (dir==5){minI=maxI;useI=0;}//+i
			if (dir==2){replaceJ=maxJ+1;maxJ=minJ;useJ=0;}//-j			
			if (dir==3){minJ=maxJ;useJ=0;}//+j
			if (dir==0){replaceK=maxK+1;maxK=minK;useK=0;}//-k			
			if (dir==1){minK=maxK;useK=0;}//+k
			//cout<<getRank()<<dir<<" "<<minI<<maxI<<minJ<<maxJ<<minK<<maxK<<"\n";



			//ensure no overlap in each batch of Comm.
			int messageTag=-1;

			//b/c the tag is based on the recipient, we need to find it's destination
			//the use vars are set to 0 to override that dimension with the face's value

		
			int itr = 0;
			int msgSize=(mOrder+1)*(mOrder+1)*(mOrder+1)*(maxI-minI+1)*(maxJ-minJ+1)*(maxK-minK+1)*7;
			for(int i=minI;i<=maxI;i++)
			{
				for(int j=minJ;j<=maxJ;j++)
				{
					for(int k=minK;k<=maxK;k++)
					{
						//if useX: then use the current X value, else use the replaceX value
						//the 10^x are to seperate the indices into chunks of 3 digits to pack into an int
						messageTag=i*useI*10000000+j*useJ*10000+k*useK*10+dir;
						//(1+useX)%2 switches 0 and 1
						messageTag+=replaceI*((1+useI)%2)*10000000+replaceJ*((1+useJ)%2)*10000+replaceK*((1+useK)%2)*10;
						itr=mElements[i][j][k].writeToBlockBuffer(msgBuffer,itr);
					}
				}
			}
			MPI_Request request;
			requests.push_back(request);
			MPI_Isend(msgBuffer, msgSize, MPI_FLOAT, toRank, dir, MPI_COMM_WORLD, &request);
		}
	}
	void Block::receiveFaceInfo(unsigned int dir, vector<MPI_Request> requests)//uni-directional, doesn't check for neighbor existence
	{
		//must inverse dir (using JFK's order)
		int fromDir = 2*(dir/2)+(dir+1)%2;
		if(neighborExists[fromDir])
		{

			int fromRank = getNeighborRank(fromDir);
			//CALCULATE WHICH ARE ON FACE
			//receiving so halos are included
			int minI=1,minJ=1,minK=1;
			int maxI=mNumElements[0],maxJ=mNumElements[1],maxK=mNumElements[2];
			//receiving so side furthest in fromDir
			if (fromDir==4)minI=maxI=0;//-i			
			if (fromDir==5)minI=maxI=maxI+1;//+i
			if (fromDir==2)minJ=maxJ=0;//-j			
			if (fromDir==3)minJ=maxJ=maxJ+1;//+j
			if (fromDir==0)minK=maxK=0;//-k			
			if (fromDir==1)minK=maxK=maxK+1;//+k
			int messageTag=-1;

			int itr = 0;
			int msgSize=(mOrder+1)*(mOrder+1)*(mOrder+1)*(maxI-minI+1)*(maxJ-minJ+1)*(maxK-minK+1)*7;
			MPI_Request request;
			requests.push_back(request);
			MPI_Status status;
			MPI_Recv(msgBuffer, msgSize, MPI_FLOAT, fromRank, dir, MPI_COMM_WORLD, &status);
			for(int i=minI;i<=maxI;i++)
			{
				for(int j=minJ;j<=maxJ;j++)
				{
					for(int k=minK;k<=maxK;k++)
					{	
						
						messageTag=i*10000000+j*10000+k*10+dir;
						itr=mElements[i][j][k].readFromBlockBuffer(msgBuffer,itr);
					}
				}
			}
		}
		
	}
