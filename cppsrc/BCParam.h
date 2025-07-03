/*
BCParam.h
struct stores parameters for BC
Drew Murray
4/4/16
*/

//circular build protection
#ifndef BCPARAM_H
#define BCPARAM_H

//BE SURE TO UPDATE DEFAULT IN BCPARAM.CPP IF MORE MEMBER VARIABLES ARE ADDED
struct BCParam{

	//types (all false = no flux condition)
		bool IS_PERIODIC[3];

	//sponge layer
		unsigned int spongeWidth[6];
		float sigmaMax;

		bool IC_X_GAUSS;
		bool IC_Y_GAUSS;
		bool IC_Z_GAUSS;
		float IC_drho_amplitude;
		float IC_vx_amplitude;
		float IC_vy_amplitude;
		float IC_vz_amplitude;
		float IC_sigma;
		float IC_center_x;
		float IC_center_y;
		float IC_center_z;


}
extern const BCParamDefault;



#endif
