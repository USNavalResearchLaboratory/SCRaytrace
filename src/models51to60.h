/** \file models51to60.h
 * \brief Model 51 to 60
 */

#ifndef MODELS51TO60_H
#define MODELS51TO60_H

#include "CModelBase.h"

//! \brief cme_demi_sphere: for comparison between Marseille and Arnaud's renderers
//!
//!  Created on 14/09/2005 by M. Burtin and F. Saez
class CModel51 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};


//! \brief Tube shell model based on CModel47. Tear drop cross section. Inner density not null.
class CModel52 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);

	protected:
		float rb,alpha,rf,rs,rc,halfdist,ratio,neinit,wrapcoeff,ratio1,ratio2,stiffness,skinsigmain,skinsigmafr,ldecdens;
};


//! \brief Tube shell model based on CModel45, different cross section radius calculation
class CModel53 : public CModelBase
{
	public:
		float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
	
	protected:
		float rb,alpha,rf,ratio,nemin,wrapcoeff,h0,rs,rc,skinsigmain,skinsigmafr,skinthresup,skinthresdw;
		Cvec cvo,cva;
};


//! \brief GCS model, as published in ApJs Volume 194, Issue 2, article id. 33, 6 pp. (2011)
class CModel54 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
	
	protected:
   float rb,alpha,h,kappa,k2,mk2,tan_delta,tan_alpha,nemin,thick,b,rho,rc1,skinsigmain,skinsigmafr,skinsigmaincut,skinsigmafrcut,neaxis;
   
   Cvec O,Vleg;

};


//! \brief Simple prominence material simulation based on CModel54. The density is placed following the loop axis.
class CModel55 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
    float rb,alpha,rf,ratio,exponent,k2,mk2,nemin,thick,rs,rs2,rc,rc1,rc2;
};


//! \brief Test of a bended streamer slab
class CModel56 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
    float h,s,thetam,theta0,hthick,hang,nefactor;
};


//! \brief Constant and uniform density and temperature
class CModel57 : public CModelBase
{
  public:
    float Density(const Cvec &v);
		float Density(const Cvec &v,float &temperature);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
    float unifdens;
		float uniftemp;
};


//! \brief Blob
class CModel58 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
    float blobdens,hwidth,lon,lat,alt;
    Cvec c; // center of the blob in cartesian coordinates
};


//! \brief Blob and tail (not sure this works)
class CModel59 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    float ComputeTime(const float &lon);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
    
  protected:
    float blobdens,hwidth,lon,lat,alt,pc[4],speed,timenucpos;
    Cvec c; // center of the blob in cartesian coordinates
};


//! \brief Tube shell model based on CModel54, but the density does not go down to 0 within the shell.
class CModel60 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
    
  protected:
    float rb,alpha,rf,ratio,k2,mk2,nemin,thick,h0,rs,rs2,rc,rc1,rc2,skinsigmain,skinsigmafr,neinside;

};

#endif
