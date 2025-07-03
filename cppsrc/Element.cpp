/*
Element.cpp
element in DG
Drew Murray
6/30/16
*/

#define loopinternalnodes(ngl) for (int __i = 0; __i < ngl; __i++) for (int __j = 0; __j < ngl; __j++) for (int __k = 0; __k < ngl; __k++)


//class includes
#include "Element.h"
#include "Node.h"
#include "Basis.h"
#include "Block.h"
#include "Matrix.h"
#include "kernelFuncs.h"


//libraries
#include <iostream>
using std::cout;
#include <vector>
using std::vector;
#include <string>
using std::to_string;
#include <cmath>
#include <iomanip>

//constructors BE SURE TO UPDATE ALL OF THEM IF MEMBER DATA IS CHANGED
	Element::Element(const Element& e)
	{
		for (unsigned int iface = 0; iface<6; iface++) neighbors[iface] = NULL;


		mGlobalAddress=e.mGlobalAddress;
		mNodes=e.mNodes;
		std::copy(e.mLocalPosition,e.mLocalPosition+3,mLocalPosition);
		std::copy(e.mGlobalPosition,e.mGlobalPosition+3,mGlobalPosition);
		std::copy(e.neighbors,e.neighbors +6, neighbors);
		mOrder=e.mOrder;
		std::copy(e.mCorner,e.mCorner+3,mCorner);
		std::copy(e.mWidth,e.mWidth+3,mWidth);
		mRho0 = e.mRho0;
		mC0 = e.mC0;
		mBoverA = e.mBoverA;
		mSigmaU = e.mSigmaU;
		mSigmaDRho = e.mSigmaDRho;
		mDiffusionConstant = e.mDiffusionConstant;
		basisPtr = e.basisPtr;
		BCParamPtr=e.BCParamPtr;
		massPtr=e.massPtr;
		stiffPtr=e.stiffPtr;
		mRank = e.mRank;
		blockPtr = NULL;//force to use setBlockPtr()
		std::copy(e.mNumElem, e.mNumElem + 3, mNumElem);
		std::copy(e.mIsGridEdge, e.mIsGridEdge + 6, mIsGridEdge);
		std::copy(e.msgBuffer, e.msgBuffer + 6, msgBuffer);




		for(unsigned int i=0;i<mNodes.size();i++)
		{
			for(unsigned int j=0;j<mNodes[0].size();j++)
			{
				for(unsigned int k=0;k<mNodes[0][0].size();k++)
				{
					mNodes[i][j][k].setElemPtr(this);
				}
			}
		}
		
	}
	Element& Element::operator=(const Element& e)
	{
		for (unsigned int iface = 0; iface<6; iface++) neighbors[iface] = NULL;
			
		


		mGlobalAddress=e.mGlobalAddress;
		mNodes=e.mNodes;
		std::copy(e.mLocalPosition,e.mLocalPosition+3,mLocalPosition);
		std::copy(e.mGlobalPosition,e.mGlobalPosition+3,mGlobalPosition);
		std::copy(e.neighbors, e.neighbors + 3, neighbors);
		mOrder=e.mOrder;
		std::copy(e.mCorner,e.mCorner+3,mCorner);
		std::copy(e.mWidth,e.mWidth+3,mWidth);
		mRho0 = e.mRho0;
		mC0 = e.mC0;
		mBoverA = e.mBoverA;
		mSigmaU = e.mSigmaU;
		mSigmaDRho = e.mSigmaDRho;
		mDiffusionConstant = e.mDiffusionConstant;
		std::copy(e.mNumElem, e.mNumElem + 3, mNumElem);
		basisPtr = e.basisPtr;
		BCParamPtr=e.BCParamPtr;
		massPtr = e.massPtr;
		stiffPtr = e.stiffPtr;
		mRank = e.mRank;
		std::copy(e.mIsGridEdge, e.mIsGridEdge + 6, mIsGridEdge);
		std::copy(e.msgBuffer, e.msgBuffer + 6, msgBuffer);
		blockPtr = NULL;//force to use setBlockPtr()

		for(unsigned int i=0;i<mNodes.size();i++)
		{
			for(unsigned int j=0;j<mNodes[0].size();j++)
			{
				for(unsigned int k=0;k<mNodes[0][0].size();k++)
				{
					mNodes[i][j][k].setElemPtr(this);
				}
			}
		}
		return *this;
	}

	Element::Element(BCParam* BCPtr, unsigned int* numelem, unsigned int globalAddress, int* localPosition, int* globalPosition, unsigned short int order, float* corner, const float* width, Basis* paramBasisPtr, float Rho0, float C0, float BoverA, float diffusionConstant,float sigmaU, float sigmaDRho, int rank, Block* paramBlock, bool* pIsGridEdge, Matrix* massPointer, Matrix* stiffPointer)
	{
		for (unsigned int iface = 0; iface<6; iface++) neighbors[iface] = NULL;

	//simply move from parameters to member variables
		mGlobalAddress=globalAddress;
		std::copy(localPosition,localPosition+3,mLocalPosition);
		std::copy(globalPosition,globalPosition+3,mGlobalPosition);
		mOrder=order;
		std::copy(corner,corner+3,mCorner);
		std::copy(width,width+3,mWidth);
		basisPtr = paramBasisPtr;
		mRho0 = Rho0;
		mC0 = C0;
		mBoverA = BoverA;
		mSigmaU = sigmaU;
		mSigmaDRho = sigmaDRho;
		mDiffusionConstant = diffusionConstant;
		mRank = rank;
		blockPtr = paramBlock;
		BCParamPtr=BCPtr;
		std::copy(pIsGridEdge, pIsGridEdge + 6, mIsGridEdge);
		std::copy(numelem, numelem + 3, mNumElem);
		massPtr = massPointer;
		stiffPtr = stiffPointer;


		for(int i =0; i<6;i++)
		{
			neighbors[i]=NULL;
		}


	//temporary variables to initialize nodes
		float coord[3]={0.,0.,0.};
		unsigned short int elemCoord[3] = { 0, 0, 0 };
		float velocity[3]={0.,0.,0.};
		float Rho = 0;
		
	
	//coefficients of polynomial basis
		vector<double> xCoeff=*(basisPtr->getZeros(order));
		vector<double> yCoeff=*(basisPtr->getZeros(order));
		vector<double> zCoeff=*(basisPtr->getZeros(order));

		
	//create nodes
		for(unsigned int i = 0;i<(unsigned int)(order+1);i++)//comparison must be unsigned to prevent warnings
		{
			vector<vector<Node>> plane;//Oject oriented 2d array
			for(unsigned int j = 0;j<(unsigned int)(order+1);j++)
			{
				vector<Node> line;//Object oriented 1d array
				for(unsigned int k = 0;k<(unsigned int)(order+1);k++)
				{
					//coefficients are between [-1,1], they must be renormalized to [0,1]
					coord[0]=width[0]*((float)xCoeff[i]/2.F+.5F)+mCorner[0];
					coord[1]=width[1]*((float)yCoeff[j]/2.F+.5F)+mCorner[1];
					coord[2]=width[2]*((float)zCoeff[k]/2.F+.5F)+mCorner[2];
					elemCoord[0]=i;
					elemCoord[1]=j;
					elemCoord[2]=k;
					

					Node newNode(BCPtr, this, coord, elemCoord, velocity, Rho);
					line.push_back(newNode);
				}
				plane.push_back(line);
			}
			mNodes.push_back(plane);
		}
		if (isHalo())
		{
			mIsGridEdge[0]=false;
			mIsGridEdge[1]=false;
			mIsGridEdge[2]=false;
			mIsGridEdge[3]=false;
			mIsGridEdge[4]=false;
			mIsGridEdge[5]=false;
		}
	}

	string Element::allDRhoToString(bool twoDOnly)
	{
		std::stringstream ss;//using stringstream to convert float to string with selectable precision
		string ret;
		for (unsigned int i = 0;i < mNodes.size();i++)
		{
			for (unsigned int j = 0;j < mNodes[i].size();j++)
			{
				for (unsigned int k = 0;k < mNodes[i][j].size();k++)
				{

					if (!twoDOnly || k == 0)
					{
						ss << std::fixed << std::setprecision(15) << mNodes[i][j][k].getRho();
						ret = ret + "      Node " + std::to_string(i) + "," + std::to_string(j) + "," + std::to_string(k) + ": " + ss.str() + "\n";
						ss.str(string());//clear the stream
					}
				}
			}
		}
		return ret;
	}

	bool Element::isHalo()
	{
		bool ret = false;
		for(int i=0;i<3;i++)
		{
			if (mLocalPosition[i]==signed(mNumElem[i]) or mLocalPosition[i]==-1)
			{
				ret = true;
			}
		}
		return ret;
	}
	
	//copy state variables from another element into this one
	void Element::setNodeData(Element* el)
	{
		for(int i=0;i<mOrder+1;i++)
		{
			for(int j=0;j<mOrder+1;j++)
			{
				for(int k=0;k<mOrder+1;k++)
				{
					mNodes[i][j][k].mRho=el->mNodes[i][j][k].mRho;
					mNodes[i][j][k].mVelocity[0]=el->mNodes[i][j][k].mVelocity[0];
					mNodes[i][j][k].mVelocity[1]=el->mNodes[i][j][k].mVelocity[1];
					mNodes[i][j][k].mVelocity[2]=el->mNodes[i][j][k].mVelocity[2];
					mNodes[i][j][k].mQ[0] = el->mNodes[i][j][k].mQ[0];
					mNodes[i][j][k].mQ[1] = el->mNodes[i][j][k].mQ[1];
					mNodes[i][j][k].mQ[2] = el->mNodes[i][j][k].mQ[2];
				}
			}
		}
	}
	//copy physical constants from another element into this one
	void Element::setPhysicalConstants(Element* el)
	{
		mRho0 = el->mRho0;
		mC0 = el->mC0;
		mBoverA = el->mBoverA;
		mDiffusionConstant = el->mDiffusionConstant;
	}

	//multiply state variables in place by a scalar
	void Element::selfMult(float scalar)
	{
		for (int i = 0;i < mOrder + 1;i++)
		{
			for (int j = 0;j < mOrder + 1;j++)
			{
				for (int k = 0;k < mOrder + 1;k++)
				{
					mNodes[i][j][k].mRho *= scalar;
					mNodes[i][j][k].mVelocity[0] *= scalar;
					mNodes[i][j][k].mVelocity[1] *= scalar;
					mNodes[i][j][k].mVelocity[2] *= scalar;
					mNodes[i][j][k].mQ[0] *= scalar;
					mNodes[i][j][k].mQ[1] *= scalar;
					mNodes[i][j][k].mQ[2] *= scalar;
				}
			}
		}

	}
	//applies formula: self = self * selfscalar + other * otherscalar where self
	void Element::selfWeightedAdd(Element* other, float selfScalar,float otherScalar)
	{
		for (int i = 0;i < mOrder + 1;i++)
		{
			for (int j = 0;j < mOrder + 1;j++)
			{
				for (int k = 0;k < mOrder + 1;k++)
				{
					mNodes[i][j][k].mRho = mNodes[i][j][k].mRho * selfScalar + other->mNodes[i][j][k].mRho * otherScalar;
					mNodes[i][j][k].mVelocity[0] = mNodes[i][j][k].mVelocity[0] * selfScalar + other->mNodes[i][j][k].mVelocity[0] * otherScalar;
					mNodes[i][j][k].mVelocity[1] = mNodes[i][j][k].mVelocity[1] * selfScalar + other->mNodes[i][j][k].mVelocity[1] * otherScalar;
					mNodes[i][j][k].mVelocity[2] = mNodes[i][j][k].mVelocity[2] * selfScalar + other->mNodes[i][j][k].mVelocity[2] * otherScalar;
					mNodes[i][j][k].mQ[0] = mNodes[i][j][k].mQ[0] * selfScalar + other->mNodes[i][j][k].mQ[0] * otherScalar;
					mNodes[i][j][k].mQ[1] = mNodes[i][j][k].mQ[1] * selfScalar + other->mNodes[i][j][k].mQ[1] * otherScalar;
					mNodes[i][j][k].mQ[2] = mNodes[i][j][k].mQ[2] * selfScalar + other->mNodes[i][j][k].mQ[2] * otherScalar;
				}
			}
		}

	}

	
	//write node data into the block's MPI buffer, then return the index where the next element will start
	int Element::writeToBlockBuffer(float* buffer, int startIdx)
	{
		int offset=startIdx;
		int maxI=mOrder+1;
		int maxJ=mOrder+1;
		int maxK=mOrder+1;		
		for(int i=0;i<maxI;i++)
		{
			for(int j=0;j<maxJ;j++)
			{
				for(int k=0;k<maxK;k++)
				{
					buffer[offset]=mNodes[i][j][k].mRho;
					buffer[offset+1]=mNodes[i][j][k].mVelocity[0];
					buffer[offset+2]=mNodes[i][j][k].mVelocity[1];
					buffer[offset+3]=mNodes[i][j][k].mVelocity[2];
					buffer[offset+4]=mNodes[i][j][k].mQ[0];
					buffer[offset+5]=mNodes[i][j][k].mQ[1];
					buffer[offset+6]=mNodes[i][j][k].mQ[2];
					offset+=7;
				}
			}
		}
		return offset;
	}
	//read data out of the block's MPI buffer and into node data, then return the index where the next element will start
	int Element::readFromBlockBuffer(float* buffer, int startIdx)
	{
		int offset=startIdx;
		int maxI=mOrder+1;
		int maxJ=mOrder+1;
		int maxK=mOrder+1;		
		for(int i=0;i<maxI;i++)
		{
			for(int j=0;j<maxJ;j++)
			{
				for(int k=0;k<maxK;k++)
				{
					mNodes[i][j][k].mRho=buffer[offset];
					mNodes[i][j][k].mVelocity[0]=buffer[offset+1];
					mNodes[i][j][k].mVelocity[1]=buffer[offset+2];
					mNodes[i][j][k].mVelocity[2]=buffer[offset+3];
					mNodes[i][j][k].mQ[0]=buffer[offset+4];
					mNodes[i][j][k].mQ[1]=buffer[offset+5];
					mNodes[i][j][k].mQ[2]=buffer[offset+6];
					offset+=7;
					
				}
			}
		}
		return offset;
	}

	void Element::debugPrintElemAtCoord(string step, int x, int y, int z)
	{
		if (!isHalo()&&mGlobalPosition[0]==x&&mGlobalPosition[1]==y&&mGlobalPosition[2]==z)
		{
			cout<<"@ step "<<step<<", Coordinate=("<<mNodes[0][0][0].getCoordString()<<"), Rho=="<<mNodes[0][0][0].getRho()<<"\n";
		}
	}	
	
	//set the neighbor pointer (cannot be in constructor because half the time the neighbor won't exist yet!  Must run after all construction is complete)
	void Element::initializeNeighbor(Element* el, int i)
	{
		neighbors[i] = el;
	}


	string Element::getNodeCoords()
	{
		string output="";

		for(unsigned int i=0; i<mNodes.size(); i++)
		{
			for(unsigned int j=0; j<mNodes[i].size(); j++)
			{
				for(unsigned int k=0; k<mNodes[i][j].size(); k++)
				{
					output=output+mNodes[i][j][k].getCoordString()+"\n";
				}
			}
		}
		return output;
	}

	//get Rho values along a line through the given point in the direction of the given axis (e.g. in the +/- X direction, through y=1, z=2)
	vector<double> Element::getOutputAxisRho(int alongAxis, int* idx)
	{
		vector<double> ret;
		int coord[3];
		int iMin = idx[0], jMin = idx[1], kMin = idx[2];
		int iMax = idx[0] + 1, jMax = idx[1] + 1, kMax = idx[2] + 1;
		if (alongAxis == 0) { iMin = 0; iMax = mOrder + 1; };
		if (alongAxis == 1) { jMin = 0; jMax = mOrder + 1; };
		if (alongAxis == 2) { kMin = 0; kMax = mOrder + 1; };
		for (int i = iMin; i < iMax;i++)
		{
			coord[0] = i;
			for (int j = jMin; j < jMax;j++)
			{
				coord[1] = j;
				for (int k = kMin; k < kMax;k++)
				{
					coord[2] = k;
					ret.push_back(getOutputRho(coord));
				}
			}
		}
		return ret;
	}

	//get physical coordinates along the direction of the given axis (e.g. the x values; no y,z coords must be specified as the grid is rectilinear)
	vector<double> Element::getOutputAxisCoord(int alongAxis)
	{
		int idx[3] = { 0,0,0 };//assume parallel grid lines: e.g. x coords are the same for same x-index, can measure at yidx=zidx=0
		vector<double> ret;
		int coord[3];
		int iMin = idx[0], jMin = idx[1], kMin = idx[2];
		int iMax = idx[0] + 1, jMax = idx[1] + 1, kMax = idx[2] + 1;
		if (alongAxis == 0) { iMin = 0; iMax = mOrder + 1; };
		if (alongAxis == 1) { jMin = 0; jMax = mOrder + 1; };
		if (alongAxis == 2) { kMin = 0; kMax = mOrder + 1; };
		for (int i = iMin; i < iMax;i++)
		{
			coord[0] = i;
			for (int j = jMin; j < jMax;j++)
			{
				coord[1] = j;
				for (int k = kMin; k < kMax;k++)
				{
					coord[2] = k;
					ret.push_back(getOutputCoords(coord)[alongAxis]);
				}
			}
		}
		return ret;
	}
	//get vel values (of a provided dimension) along a line through the given point in the direction of the given axis (e.g. Vz in the +/- X direction, through y=1, z=2)
	vector<double> Element::getOutputAxisVelocity(int alongAxis, int velDir, int idx[3])
	{
		vector<double> ret;

		int coord[3];
		int iMin = idx[0], jMin = idx[1], kMin = idx[2];
		int iMax = idx[0] + 1, jMax = idx[1] + 1, kMax = idx[2] + 1;
		if (alongAxis == 0) { iMin = 0; iMax = mOrder + 1; };
		if (alongAxis == 1) { jMin = 0; jMax = mOrder + 1; };
		if (alongAxis == 2) { kMin = 0; kMax = mOrder + 1; };
		for (int i = iMin; i < iMax;i++)
		{
			coord[0] = i;
			for (int j = jMin; j < jMax;j++)
			{
				coord[1] = j;
				for (int k = kMin; k < kMax;k++)
				{
					coord[2] = k;
					ret.push_back(getOutputVelocity(coord)[velDir]);
				}
			}
		}
		return ret;
	}

	void Element::OperateUFromV(unsigned short int axis, unsigned short int mode, Element* oldElem, bool leftEdgeOverride, bool rightEdgeOverride, int velocityOverrideAxis)
	{
		//left and right edge overrides are for special BC that must NOT evaluate the flux in the normal way, instead the element treats itself as its own neighbor and uses its own values.  This is done by setting the neighbor pointer.

		//mode == 0 is U from DRho; mode == 1 is DRho from U

		//OFF axis the outer loop iterates over the whole spans (both of them) and the inner loop iterates exactly once
		int i2max = mOrder + 1, j2max = mOrder + 1, k2max = mOrder + 1, i2 = 0, j2 = 0, k2 = 0, prevDir=0, nextDir=0;
		
		int nodeIdx;//points at the nodal index of the chosen dimention that is being spanned to populate the dRho and velocity arrays
		//nodeIdx counts up each cycle of the inner loop regardless of it i, j, or k is driving it

		if (axis == 0) {i2max = 1; prevDir = 4; nextDir = 5; } //ON axis, the outer loop iterates exactly once and the inner loop iterates over the whole span
		else if (axis == 1) {j2max = 1; prevDir = 2; nextDir = 3; }
		else if (axis == 2) {k2max = 1; prevDir = 0; nextDir = 1; }
		else cout << "ERROR: OPERATOR HAS INVALID AXIS PARAMETER!\n";

		int coord[3];
		int otherCoord[3];
		Matrix Flux = Matrix(mOrder + 1, 1);
		Matrix dRho = Matrix(mOrder + 1, 1);
		Matrix vel = Matrix(mOrder + 1, 1);
		Matrix result = Matrix(mOrder + 1, 1);
		Flux.Fill(0.0);
		Element* other = NULL;
		float coeff1 = -123, coeff2 = -123, coeff3 = -123, coeff4 = -123, z1 = -123, z2 = -123, r1 = -123, r2 = -123;
		float edgeVelocityFilter = 1.0;
		

		float q = oldElem->mRho0;
		if (mode == 1) q = 1.0F / oldElem->mRho0;
		Matrix(*kernelFunc)(Matrix, Matrix, float, float, float) = momentum;
		float(*kernelFuncFlux)(float, float, float, float, float) = momentum;
		float(*kernelFuncFlux2)(float, float, float, float, float) = mass;
		if (mode == 1)
		{
			kernelFunc = mass;
			kernelFuncFlux = mass;
			kernelFuncFlux2 = momentum;
		}


		for (i2 = 0;i2 < i2max;i2++)
		{
			coord[0] = i2;
			otherCoord[0] = i2;
			for (j2 = 0;j2 < j2max;j2++)
			{
				coord[1] = j2;
				otherCoord[1] = j2;
				for (k2 = 0;k2 < k2max;k2++)
				{
					coord[2] = k2;
					otherCoord[2] = k2;

					//get flux at boundary with previous element while coord corresponds to nodeIdx=0
					other = oldElem->neighbors[prevDir];
					edgeVelocityFilter = 1.0;
					//for use by BC, may need to set 'other' to 'self' on left and/or right side and use zero instead of velocity
					if (leftEdgeOverride)
					{
						other = oldElem;
						if(axis==velocityOverrideAxis) edgeVelocityFilter = 0.0;
					}

					coord[axis] = 0;
					otherCoord[axis] = mOrder;
					z1 = oldElem->mRho0 * oldElem->mC0; z2 = other->mRho0 * other->mC0; r1 = oldElem->mRho0; r2 = other->mRho0;
					if (mode == 0)
					{
						coeff1 = z1 * r2 / (z1 + z2);
						coeff2 = z2 * r1 / (z1 + z2);
						coeff3 = z1 * z2 / r2 / (z1 + z2);
						coeff4 = -1.0F * z1 * z2 / r1 / (z1 + z2);
					}
					if (mode == 1)
					{
						coeff1 = z2 / r2 / (z1 + z2);
						coeff2 = z1 / r1 / (z1 + z2);
						coeff3 = r2 / (z1 + z2);
						coeff4 = -1.0F * r1 / (z1 + z2);
					}
					Flux[0][0] = kernelFuncFlux(other->getOutputVelocity(otherCoord)[axis] * edgeVelocityFilter, other->getOutputRho(otherCoord), other->mC0, other->mRho0, other->mBoverA) * coeff1;
					Flux[0][0] += kernelFuncFlux(oldElem->getOutputVelocity(coord)[axis] * edgeVelocityFilter, oldElem->getOutputRho(coord), oldElem->mC0, oldElem->mRho0, oldElem->mBoverA) * coeff2;
					Flux[0][0] += kernelFuncFlux2(other->getOutputVelocity(otherCoord)[axis] * edgeVelocityFilter, other->getOutputRho(otherCoord), other->mC0, other->mRho0, other->mBoverA) * coeff3;
					Flux[0][0] += kernelFuncFlux2(oldElem->getOutputVelocity(coord)[axis] * edgeVelocityFilter, oldElem->getOutputRho(coord), oldElem->mC0, oldElem->mRho0, oldElem->mBoverA) * coeff4;
					Flux[0][0] *= -1.;
					//populate DRho and U along the axis
					for (nodeIdx = 0; nodeIdx < mOrder + 1;nodeIdx += 1)
					{
						coord[axis] = nodeIdx;
						dRho[nodeIdx][0] = oldElem->getOutputRho(coord);
						vel[nodeIdx][0] = oldElem->getOutputVelocity(coord)[axis] * edgeVelocityFilter;
					}


					//get flux at boundary with next element after coord corresponds to nodeIdx=max
					other = oldElem->neighbors[nextDir];
					edgeVelocityFilter = 1.0;
					//for use by BC, may need to set 'other' to 'self' on left and/or right side and use zero instead of velocity
					if (rightEdgeOverride)
					{
						other = oldElem;
						if (axis == velocityOverrideAxis) edgeVelocityFilter = 0.0;
					}
					otherCoord[axis] = 0;
					z2 = other->mRho0 * other->mC0; r2 = other->mRho0;
					if (mode == 0)
					{
						coeff1 = z2 * r1 / (z1 + z2);
						coeff2 = z1 * r2 / (z1 + z2);
						coeff3 = z1 * z2 / r1 / (z1 + z2);
						coeff4 = -1.0F * z1 * z2 / r2 / (z1 + z2);
					}
					if (mode == 1)
					{
						coeff1 = z1 / r1 / (z1 + z2);
						coeff2 = z2 / r2 / (z1 + z2);
						coeff3 = r1 / (z1 + z2);
						coeff4 = -1.0F * r2 / (z1 + z2);
					}

					Flux[mOrder][0] = kernelFuncFlux(oldElem->getOutputVelocity(coord)[axis] * edgeVelocityFilter, oldElem->getOutputRho(coord), oldElem->mC0, oldElem->mRho0, oldElem->mBoverA) * coeff1;
					Flux[mOrder][0] += kernelFuncFlux(other->getOutputVelocity(otherCoord)[axis] * edgeVelocityFilter, other->getOutputRho(otherCoord), other->mC0, other->mRho0, other->mBoverA) * coeff2;
					Flux[mOrder][0] += kernelFuncFlux2(oldElem->getOutputVelocity(coord)[axis] * edgeVelocityFilter, oldElem->getOutputRho(coord), oldElem->mC0, oldElem->mRho0, oldElem->mBoverA) * coeff3;
					Flux[mOrder][0] += kernelFuncFlux2(other->getOutputVelocity(otherCoord)[axis] * edgeVelocityFilter, other->getOutputRho(otherCoord), other->mC0, other->mRho0, other->mBoverA) * coeff4;
					
					result = ((*massPtr) * q * oldElem->mWidth[axis]) / ((*stiffPtr) * q * kernelFunc(vel, dRho, oldElem->mC0, oldElem->mRho0, oldElem->mBoverA) - Flux);
					
					//store DRho and U along the axis
					for (nodeIdx = 0; nodeIdx < mOrder + 1;nodeIdx += 1)
					{
						coord[axis] = nodeIdx;
						if (mode == 0)
						{
							addVelocityComponent(coord, axis, (float)result[0][nodeIdx]);
						}
						else if(mode == 1)
						{
							addRho(coord, (float)result[0][nodeIdx]);
						}
					}
				}
			}
		}
	}


	void Element::OperateQFromDRho(unsigned short int axis, unsigned short int mode, Element* oldElem)
	{
		//OFF axis the outer loop iterates over the whole spans (both of them) and the inner loop iterates exactly once
		int i2max = mOrder + 1, j2max = mOrder + 1, k2max = mOrder + 1, i2 = 0, j2 = 0, k2 = 0, prevDir = 0, nextDir = 0;

		int nodeIdx;//points at the nodal index of the chosen dimention that is being spanned to populate the dRho and velocity arrays
		//nodeIdx counts up each cycle of the inner loop regardless of it i, j, or k is driving it

		if (axis == 0) { i2max = 1; prevDir = 4; nextDir = 5; } //ON axis, the outer loop iterates exactly once and the inner loop iterates over the whole span
		else if (axis == 1) { j2max = 1; prevDir = 2; nextDir = 3; }
		else if (axis == 2) { k2max = 1; prevDir = 0; nextDir = 1; }
		else cout << "ERROR: OPERATOR HAS INVALID AXIS PARAMETER!\n";

		int coord[3];
		int otherCoord[3];
		Matrix Flux = Matrix(mOrder + 1, 1);
		Matrix dRho = Matrix(mOrder + 1, 1);
		Matrix result = Matrix(mOrder + 1, 1);
		Flux.Fill(0.0);
		Element* other = NULL;

		for (i2 = 0;i2 < i2max;i2++)
		{
			coord[0] = i2;
			otherCoord[0] = i2;
			for (j2 = 0;j2 < j2max;j2++)
			{
				coord[1] = j2;
				otherCoord[1] = j2;
				for (k2 = 0;k2 < k2max;k2++)
				{
					coord[2] = k2;
					otherCoord[2] = k2;

					//get flux at boundary with previous element while coord corresponds to nodeIdx=0
					other = oldElem->neighbors[prevDir];
					coord[axis] = 0;
					otherCoord[axis] = mOrder;
					Flux[0][0] = other->getOutputRho(otherCoord);
					Flux[0][0] += oldElem->getOutputRho(coord);
					Flux[0][0] *= 0.5F;
					//populate DRho along the axis
					for (nodeIdx = 0; nodeIdx < mOrder + 1;nodeIdx += 1)
					{
						coord[axis] = nodeIdx;
						dRho[nodeIdx][0] = oldElem->getOutputRho(coord);
					}

					//get flux at boundary with next element after coord corresponds to nodeIdx=max
					other = oldElem->neighbors[nextDir];
					//for use by BC, may need to set 'other' to 'self' on left and/or right side and use zero instead of velocity
					otherCoord[axis] = 0;

					Flux[mOrder][0] = oldElem->getOutputRho(coord);
					Flux[mOrder][0] += other->getOutputRho(otherCoord);
					Flux[mOrder][0] *= -0.5F;
					result = ((*massPtr) * oldElem->mWidth[axis]) / ((*stiffPtr) * dRho + Flux);

					//store Q along the axis
					for (nodeIdx = 0; nodeIdx < mOrder + 1;nodeIdx += 1)
					{
						coord[axis] = nodeIdx;
						setQComponent(coord, axis, -1.0F*(float)result[0][nodeIdx]);
					}
				}
			}
		}
	}
	void Element::OperateDRhoFromQ(unsigned short int axis, unsigned short int mode, Element* oldElem)
	{
		//OFF axis the outer loop iterates over the whole spans (both of them) and the inner loop iterates exactly once
		int i2max = mOrder + 1, j2max = mOrder + 1, k2max = mOrder + 1, i2 = 0, j2 = 0, k2 = 0, prevDir = 0, nextDir = 0;

		int nodeIdx;//points at the nodal index of the chosen dimention that is being spanned to populate the dRho and velocity arrays
		//nodeIdx counts up each cycle of the inner loop regardless of it i, j, or k is driving it

		if (axis == 0) { i2max = 1; prevDir = 4; nextDir = 5; } //ON axis, the outer loop iterates exactly once and the inner loop iterates over the whole span
		else if (axis == 1) { j2max = 1; prevDir = 2; nextDir = 3; }
		else if (axis == 2) { k2max = 1; prevDir = 0; nextDir = 1; }
		else cout << "ERROR: OPERATOR HAS INVALID AXIS PARAMETER!\n";

		int coord[3];
		int otherCoord[3];
		Matrix Flux = Matrix(mOrder + 1, 1);
		Matrix Q = Matrix(mOrder + 1, 1);
		Matrix result = Matrix(mOrder + 1, 1);
		Flux.Fill(0.0);
		Element* other = NULL;

		for (i2 = 0;i2 < i2max;i2++)
		{
			coord[0] = i2;
			otherCoord[0] = i2;
			for (j2 = 0;j2 < j2max;j2++)
			{
				coord[1] = j2;
				otherCoord[1] = j2;
				for (k2 = 0;k2 < k2max;k2++)
				{
					coord[2] = k2;
					otherCoord[2] = k2;

					//get flux at boundary with previous element while coord corresponds to nodeIdx=0
					other = oldElem->neighbors[prevDir];
					coord[axis] = 0;
					otherCoord[axis] = mOrder;
					Flux[0][0] = other->getOutputQ(otherCoord)[axis];
					Flux[0][0] += oldElem->getOutputQ(coord)[axis];
					Flux[0][0] *= 0.5F;
					//populate Q along the axis
					for (nodeIdx = 0; nodeIdx < mOrder + 1;nodeIdx += 1)
					{
						coord[axis] = nodeIdx;
						Q[nodeIdx][0] = oldElem->getOutputQ(coord)[axis];
					}

					//get flux at boundary with next element after coord corresponds to nodeIdx=max
					other = oldElem->neighbors[nextDir];
					//for use by BC, may need to set 'other' to 'self' on left and/or right side and use zero instead of velocity
					otherCoord[axis] = 0;

					Flux[mOrder][0] = oldElem->getOutputQ(coord)[axis];
					Flux[mOrder][0] += other->getOutputQ(otherCoord)[axis];
					Flux[mOrder][0] *= -0.5F;
					result = ((*massPtr) * oldElem->mWidth[axis]) / (((*stiffPtr) * Q + Flux)* -1.0F * oldElem->getDiffusionConstant());

					//store DRho along the axis
					for (nodeIdx = 0; nodeIdx < mOrder + 1;nodeIdx += 1)
					{
						coord[axis] = nodeIdx;
						addRho(coord, (float)result[0][nodeIdx]);
					}
				}
			}
		}
	}

	void Element::OperateSponge(unsigned short int axis, unsigned short int mode, Element* oldElem)
	{
		string elemtest = getName();
		int coord[3] = {-1, -1, -1};
		for (int i = 0;i < mOrder+1;i++)
		{
			coord[0] = i;
			for (int j = 0;j < mOrder + 1;j++)
			{
				coord[1] = j;
				for (int k = 0;k < mOrder + 1;k++)
				{
					coord[2] = k;
					for (int itraxis = 0; itraxis < 3; itraxis++)
					{
						addVelocityComponent(coord, itraxis,-1.0F*oldElem->getOutputVelocity(coord)[itraxis]*getSigmaU());
					}
					addRho(coord, -1.0F*oldElem->getOutputRho(coord)* getSigmaDRho());
				}
			}
		}
	}
