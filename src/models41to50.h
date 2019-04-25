#ifndef MODELS41TO50_H
#define MODELS41TO50_H

#include "CModelBase.h"

//! Tube shell model based on CModel33 with extra parameters to set the Ne
class CModel41 : public CModelBase
{
	public:
  	float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Tube shell model based on CModel41 with different axis orientation
//!
//! Transition for thickness at the foot-circle junction
//!
//! !!!! Changed the axis orientation to facilitate the positioning
//!      - Oyz is now the plane of the loop
//!      - Oz is the principal axis of the loop
class CModel42 : public CModelBase
{
	public:
  	float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Tube with an helix
class CModel43 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Helicoidal ribbon
//!
//! !!!! Changed the axis orientation to facilitate the positioning
//!      - Oyz is now the plane of the loop
//!      - Oz is the principal axis of the loop
//!      
//! pparam :
//! - [0] : rb : 1.5 [Rsun] : starting point of the structure
//! - [1] : alpha : 30*DTOR [rad] : angle between axis and feet
//! - [2] : rf : 10 [Rsun] : dist junction line/circle
//! - [3] : 1e5 : Ne min
//! - [4] : omegah : 0.33 [rad/Rsun] : pulsation of the helix
//! - [5] : helphase : PI/2 [rad] : phase of the helix
//! - [6] : xhalfthick : 2.5 [Rsun] : half thickness of the ribbon on X axis
//! - [7] : yhalfthick : 0.3 [Rsun] : half thickness of the ribbon on Y axis
class CModel44 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Tube shell model based on CModel42 but the shell skin has a gaussian profile.
//!
//! The edge factor of CModel42 is replaced by the sigma parameter (well, 
//! not exactly, see the code) of the gaussian profile.
//!
//! Axis orientation:
//!      - Oyz is now the plane of the loop
//!      - Oz is the principal axis of the loop
class CModel45 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
	
	protected:
		float rb,alpha,rf,ratio,nemin,wrapcoeff,h0,rs,rc,skinsigmain,skinsigmafr;
};
//! Tube shell model based on CModel45, with decreasing of density depending on height.
class CModel47 : public CModelBase
{
	public:
	  	float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);

	protected:
		float rb,alpha,rf,rs,rc,halfdist,ratio,neinit,wrapcoeff,ratio1,ratio2,stiffness,skinsigmain,skinsigmafr,ldecdens;
};
//! Spherically symmetric Ne, based on CModel14
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
class CModel49 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Slab model with modulation of the slab FO. Based on CModel38.
//!
//! Spherically symmetric Ne
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
//! Ne parameters passing for pminimizerXX.pro (XX=  )
//! - pparam[0..4] : coeff of the Ne polynomial
//! - pparam[5] : slab half angle
//! - pparam[6] : subscript max of the modulation profile (size - 1)
//! - pparam[7] : first element of the modulation profile
class CModel50 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};

#endif

