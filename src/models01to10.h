#ifndef _MODELS01TO10_H_
#define _MODELS01TO10_H_

#include <vector>
#include "Cvec.h"
#include "CModelBase.h"

//!M.Guhathakurta model. 
//!
//! From: "The Large-scale density structure of the solar corona and the heliospheric surrent sheet", ApJ, 458:817-831, 1996 Feb 20

class CModel01 : public CModelBase
{
	public:
	  float Density(const Cvec &v,float &temperature);
		//int getNbparam() {return 0;};
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Model for testing purpose
class CModel02 : public CModelBase
{
	public:
	  	float Density(const Cvec &v,float &temperature); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Model for testing purpose
class CModel03 : public CModelBase
{
	public:
	  	float Density(const Cvec &v,float &temperature); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Density cube from Y.M. with trilinear interpolation
class CModel04 : public CModelBase
{
	public:
		float Density(const Cvec &v,float &temperature); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
	protected:
		float* pmparam;
};
//! Same as Model 4 but without trilinear interpolation
class CModel05 : public CModelBase
{
	public:
	  float Density(const Cvec &v,float &temperature); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
	protected:
		float* pmparam;
};
//! Very simple model of streamer belt
class CModel06 : public CModelBase
{
	public:
	  	float Density(const Cvec &v,float &temperature); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Slab model of streamer based on Guhathakurta model.
class CModel07 : public CModelBase
{
	public:
	  	float Density(const Cvec &v,float &temperature); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};

//! Sphericaly symmetric Ne
//!
//! Use the Hayes model:
//!
//! Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//!
class CModel08 : public CModelBase
{
	public:
	  	float Density(const Cvec &v,float &temperature); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Guhatahakurta radial and Vibert shape.
//!
//! - Ne=Nradial(r) x Nshape(d)
//!
//! - Nshape(d)=exp(-(d/d0)^4)
//!
class CModel09 : public CModelBase
{
	public:
	  	float Density(const Cvec &v,float &temperature); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Sphericaly symmetric Ne
//!
//! use the Hayes model
//!
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-abs(d/d0))
//!
class CModel10 : public CModelBase
{
	public:
	  	float Density(const Cvec &v,float &temperature); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};

#endif
