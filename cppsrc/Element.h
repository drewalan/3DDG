/*
Element.h
element in DG
Drew Murray
6/30/16
*/
#define N_CANONICAL_FACES 6
#define N_DIMS 3
#define MAXMSGSIZE 7000
#define PI 3.14159265358979323846
//circular build protection
#ifndef ELEMENT_H
#define ELEMENT_H

//libraries
#include <vector>
#include <string>
using std::vector;
#include <sstream>

//Some compilers want the top version of SSTR, some want the bottom (also must be swapped in Block.h)
//#define SSTR( x ) ( std::ostringstream() << std::dec << x ).str()
#define SSTR( x ) static_cast< std::ostringstream & >(( std::ostringstream() << std::dec << x ) ).str()

//class includes
#include "Node.h"
#include "Basis.h"
#include "Block.h"
#include "Matrix.h"

//forward declarations (only pointers used)
class Block;
struct BCParam;

class Element
{
	public:

	//constructors/destructors/etc.
		Element(BCParam* BCPtr, unsigned int* numelem, unsigned int globalAddress, int* localPosition, int* globalPosition, unsigned short int order, float* corner, const float* width, Basis* basisPtr, float Rho0, float C0, float BoverA, float diffusionConstant, float sigmaU, float sigmaDRho, int rank, Block* paramBlock, bool* pIsGridEdge, Matrix* massPtr, Matrix* stiffPtr);
		~Element(){};
		Element(const Element& e);//copy constructor
		Element& operator=(const Element& g);//assignment operator



		
	//member functions

		//CRITICAL FUNCTIONS called from Block, all others are helper functions
		void initializeNeighbor(Element* el, int i);
		void sendFaceInfo(unsigned int dir, int tag, vector<MPI_Request> requests);
		void receiveFaceInfo(unsigned int dir, int tag, vector<MPI_Request> requests);
		int writeToBlockBuffer(float* buffer, int startIdx);
		int readFromBlockBuffer(float* buffer, int startIdx);
		void setNodeData(Element* el);
		void setPhysicalConstants(Element* el);
		string allDRhoToString(bool twoDOnly);

		bool isHalo();

		void debugPrintElemAtCoord(string step, int x, int y, int z);

		//returns unique identifier as string
		string getName(){
			return "EL" + SSTR(mGlobalAddress) + ":" + SSTR(mGlobalPosition[0]) + \
				"," + SSTR(mGlobalPosition[1]) + "," + SSTR(mGlobalPosition[2]);
		}
		string getLocalName(){
			return "(" + SSTR(mLocalPosition[0]) + \
				"," + SSTR(mLocalPosition[1]) + "," + SSTR(mLocalPosition[2])+")";
		}
		string getNodeCoords();
		int getRank(){return mRank;};
		bool isGridEdge(int dir){return mIsGridEdge[dir];}
		void setBlockPtr(Block* ptr){blockPtr=ptr;}

		void selfMult(float scalar);
		void selfWeightedAdd(Element* other, float selfScalar, float otherScalar);
		void OperateUFromV(unsigned short int axis, unsigned short int mode, Element* oldElem, bool leftEdgeOverride, bool rightEdgeOverride, int velocityOverrideAxis);
		void OperateSponge(unsigned short int axis, unsigned short int mode, Element* oldElem);
		void OperateQFromDRho(unsigned short int axis, unsigned short int mode, Element* oldElem);
		void OperateDRhoFromQ(unsigned short int axis, unsigned short int mode, Element* oldElem);

		void setVelocity(int idx[3], float* vel) { mNodes[idx[0]][idx[1]][idx[2]].setVelocity(vel); };
		void setVelocityComponent(int idx[3], int axis, float vel) { mNodes[idx[0]][idx[1]][idx[2]].setVelocityComponent(axis, vel); };
		void setQComponent(int idx[3], int axis, float vel) { mNodes[idx[0]][idx[1]][idx[2]].setQComponent(axis, vel); };
		void addVelocityComponent(int idx[3], int axis, float vel) { mNodes[idx[0]][idx[1]][idx[2]].addVelocityComponent(axis, vel); };
		void multVelocityComponent(int idx[3], int axis, float vel) { mNodes[idx[0]][idx[1]][idx[2]].multVelocityComponent(axis, vel); };
		void setRho(int idx[3], float rho) { mNodes[idx[0]][idx[1]][idx[2]].setRho(rho); };
		void addRho(int idx[3], float rho) { mNodes[idx[0]][idx[1]][idx[2]].addRho(rho); };
		void multRho(int idx[3], float rho) { mNodes[idx[0]][idx[1]][idx[2]].multRho(rho); };
		float getOutputRho(int idx[3]){ return mNodes[idx[0]][idx[1]][idx[2]].getRho(); };
		float* getOutputVelocity(int idx[3]){ return mNodes[idx[0]][idx[1]][idx[2]].getVelocity(); };
		float* getOutputQ(int idx[3]){ return mNodes[idx[0]][idx[1]][idx[2]].getQ(); };
		float* getOutputCoords(int idx[3]){ return mNodes[idx[0]][idx[1]][idx[2]].getCoord(); };
		vector<double> getOutputAxisRho(int alongAxis, int idx[3]);
		vector<double> getOutputAxisCoord(int alongAxis);
		vector<double> getOutputAxisVelocity(int alongAxis, int velDir, int idx[3]);

		unsigned short int getOrder(){ return mOrder; };
		float getRho0(){ return mRho0; };
		float getC0(){ return mC0; };
		float getBoverA(){ return mBoverA; };
		float getDiffusionConstant(){ return mDiffusionConstant; };
		float getSigmaU(){ return mSigmaU; };
		float getSigmaDRho(){ return mSigmaDRho; };
		void setRho0(float j) { mRho0 = j; };
		void setC0(float j) { mC0 = j; };
		void setBoverA(float j) { mBoverA = j; };
		void setSigmaU(float j) { mSigmaU = j; };
		void setSigmaDRho(float j) { mSigmaDRho = j; };
		void setDiffusionConstant(float j) { mDiffusionConstant = j; };
		void intializePhysicalParams(BCParam* BCParamPtr);


		void applyBoundaryConditions(BCParam* parameters, float t, float dt);

	private:

		unsigned short numNodes()//number of nodes in element
			{return (mOrder+1)*(mOrder+1)*(mOrder+1);};	

	public:
		vector<vector<vector<Node>>> mNodes;
	private:

		unsigned int mGlobalAddress;
		float mCorner[3];//left bottom back corner coordinates
		int mLocalPosition[3];//which # element it is in EACH dimension within a block
		int mGlobalPosition[3];//which # element it is in EACH dimension
		unsigned short int mOrder;
		int mRank;
		float mWidth[3];
		Basis* basisPtr;
		unsigned int mNumElem[3];
		Block* blockPtr;
		BCParam* BCParamPtr;
		Matrix* massPtr;
		Matrix* stiffPtr;
		bool mIsGridEdge[6];

		Element* neighbors[6];//use "haloing" so that the edge elements don't get called to check neighbors (null pointer)
		// order of neighbors:	k-1	k+1 j-1	j+1	i-1	i+1

		float mC0;
		float mRho0;
		float mBoverA;
		float mDiffusionConstant;
		float mSigmaU;
		float mSigmaDRho;

		float msgBuffer[MAXMSGSIZE];


	private:
		Element();//disable default constructor (to prevent accidental creation of defunct objects)
};
#endif
