/*
Node.cpp
Node in DG
Drew Murray
12/1/16
*/

//class includes
#include "Node.h"
#include "Element.h"
#include "Basis.h"
#include "BCParam.h"
#include <cmath>

//libraries
#include <string>
using std::string;
using std::to_string;

#include <iostream>
using std::cout;
using std::cerr;

#include <cmath>
using std::abs;

#include <fstream>
using std::ofstream;


//constructors BE SURE TO UPDATE ALL OF THEM IF MEMBER DATA IS CHANGED
	Node::Node(const Node& n)
	{
		std::copy(n.mCoord,n.mCoord+3,mCoord);
		std::copy(n.mVelocity,n.mVelocity+3,mVelocity);
		std::copy(n.mQ,n.mQ +3, mQ);
		std::copy(n.mElemCoord,n.mElemCoord+3,mElemCoord);

		//std::copy cannot copy array of un-destructables (like c-arrays)
		//std::copy(&n.mMetric[0][0],&n.mMetric[0][0]+9,&mMetric[0][0]);//each element of the first array (an array in its own right) counts as 1

		elem=NULL;//force to use setElemPtr()
		mRho = n.mRho;

	}

	Node& Node::operator=(const Node& n)
	{
		std::copy(n.mCoord,n.mCoord+3,mCoord);
		std::copy(n.mVelocity,n.mVelocity+3,mVelocity);
		std::copy(n.mQ, n.mQ + 3, mQ);
		std::copy(n.mElemCoord,n.mElemCoord+3,mElemCoord);

		elem=NULL;//force to use setElemPtr()
		mRho = n.mRho;

		return *this;

	}

	Node::Node(BCParam* BCPtr, Element* elemParam, float* coord, unsigned short int* elemCoord, float* velocity, float Rho)
	{
		std::copy(coord,coord+3,mCoord);
		std::copy(velocity,velocity+3,mVelocity);
		mQ[0]=0.;mQ[1] = 0.;mQ[2] = 0.;
		std::copy(elemCoord,elemCoord+3,mElemCoord);

		elem=elemParam;
		mRho = Rho;

		InitializeRho(BCPtr);
		InitializeVelocities(BCPtr);
	}

//methods
	void Node::InitializeRho(BCParam* BCPtr){
		if (BCPtr->IC_X_GAUSS || BCPtr->IC_Y_GAUSS || BCPtr->IC_Z_GAUSS)
		{
			float dist = 0;
			if (BCPtr->IC_X_GAUSS)dist += (mCoord[0] - BCPtr->IC_center_x) * (mCoord[0] - BCPtr->IC_center_x);
			if (BCPtr->IC_Y_GAUSS)dist += (mCoord[1] - BCPtr->IC_center_y) * (mCoord[1] - BCPtr->IC_center_y);
			if (BCPtr->IC_Z_GAUSS)dist += (mCoord[2] - BCPtr->IC_center_z) * (mCoord[2] - BCPtr->IC_center_z);
			dist = sqrt(dist);
			setRho(BCPtr->IC_drho_amplitude * exp(-(dist) * (dist) / (1.F / BCPtr->IC_sigma)));
		}
	}
	void Node::InitializeVelocities(BCParam* BCPtr){
		if (BCPtr->IC_X_GAUSS || BCPtr->IC_Y_GAUSS || BCPtr->IC_Z_GAUSS)
		{
			float InitialVelocity[3];
			float dist = 0;
			if (BCPtr->IC_X_GAUSS)dist += (mCoord[0] - BCPtr->IC_center_x) * (mCoord[0] - BCPtr->IC_center_x);
			if (BCPtr->IC_Y_GAUSS)dist += (mCoord[1] - BCPtr->IC_center_y) * (mCoord[1] - BCPtr->IC_center_y);
			if (BCPtr->IC_Z_GAUSS)dist += (mCoord[2] - BCPtr->IC_center_z) * (mCoord[2] - BCPtr->IC_center_z);
			dist = sqrt(dist);
			InitialVelocity[0] = (BCPtr->IC_vx_amplitude * exp(-(dist) * (dist) / (1.F / BCPtr->IC_sigma)));
			InitialVelocity[1] = (BCPtr->IC_vy_amplitude * exp(-(dist) * (dist) / (1.F / BCPtr->IC_sigma)));
			InitialVelocity[2] = (BCPtr->IC_vz_amplitude * exp(-(dist) * (dist) / (1.F / BCPtr->IC_sigma)));
			setVelocity(InitialVelocity);
		}
	}
	
	bool Node::isHalo()
	{
		return elem->isHalo();
	}
	float Node::checksum()
	{
		float ret=mRho;
		for(unsigned int i=0;i<3;i++)
		{
			ret+=mVelocity[i];
		}
		return ret;
	}
	void Node::debugPrintNodeAtCoord(string step, int x, int y, int z)//debugging addition
	{
		cout<<step<<x<<y<<z<<mElemCoord[0]<<mElemCoord[1]<<mElemCoord[2]<<mCoord[0]<<mCoord[1]<<mCoord[2]<<"\n";
	}
	
	
	string Node::getCoordString()
	{
		return to_string((long double)mCoord[0])+","+to_string((long double)mCoord[1])+","+to_string((long double)mCoord[2]);
	}
	string Node::getElemCoordString()
	{
		return to_string((long double)mElemCoord[0])+","+to_string((long double)mElemCoord[1])+","+to_string((long double)mElemCoord[2]);
	}

	void Node::setVelocity(float* j)
	{
		std::copy(j, j + 3, mVelocity);
	}
	void Node::setRho(float j)
	{
		mRho = j;
	}

	


