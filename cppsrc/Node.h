/*
Node.h
Node in DG (NOTE: at the faces, both elements will have a node corresponding to the same place.
			Each element will have its own member variable node.)
Drew Murray
12/1/16
*/

//circular build protection
#ifndef NODE_H
#define NODE_H

//libraries
#include <vector>
using std::vector;
#include <string>
using std::string;
using std::to_string;

#include <iostream>
using std::cout;

#include <fstream>
using std::ofstream;


//forward declarations (only pointers used)
class Element;
class Basis;
struct BCParam;

class Node
{
	friend class Element;//allow element easier access to its nodes

	public:
		//constructor
		Node(BCParam* BCPtr, Element* elemParam, float* coord, unsigned short int* elemCoord, float* velocity, float Rho);
		//copy & assignment operator
		Node(const Node& n);
		Node& operator=(const Node& d);
		void setElemPtr(Element* e){elem = e;}
	

		//access functions
		string getCoordString();
		string getElemCoordString();
		bool isHalo();
		float checksum();
		void debugPrintNodeAtCoord(string step, int x, int y, int z);

		void setRho(float j);
		void addRho(float j) { mRho += j; };
		void multRho(float j) { mRho *= j; };
		float getRho(){ return mRho; };

		void setVelocity(float* j);
		void setVelocityComponent(int dim, float j) { mVelocity[dim] = j; }
		void addVelocityComponent(int dim, float j) { mVelocity[dim] += j; }
		void multVelocityComponent(int dim, float j) { mVelocity[dim] *= j; }
		float* getVelocity(){ return mVelocity; };
		float* getQ(){ return mQ; };
		void setQComponent(int dim, float j) { mQ[dim] = j; }

		float* getCoord(){ return mCoord; };


		//Added by JK, may not be used?
		typedef float(*pointer_to_arrays)[3]; //typedefs can make things more readable with these datatypes

		
		void InitializeRho(BCParam* BCPtr);
		void InitializeVelocities(BCParam* BCPtr);

	private:
		//member variables
		float mVelocity[3];
		float mRho;
		float mCoord[3];
		float mQ[3];//intermediate variable for LDG for diffusion
		unsigned short int mElemCoord[3];
		Element* elem;


		
	private:
		Node();//disable default constructor



};
#endif
