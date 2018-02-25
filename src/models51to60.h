/***************************************************************************
 *  $Id: models51to60.h,v 1.7 2011-01-21 20:45:10 thernis Exp $
 *
 ****************************************************************************/
#ifndef MODELS51TO60_H
#define MODELS51TO60_H

#include "CModelBase.h"

//! cme_demi_sphere: for comparison between Marseille and Arnaud's renderers
//!
//!  Created on 14/09/2005 by M. Burtin and F. Saez
class CModel51 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Tube shell model based on CModel47. Tear drop cross section. Inner density not null.

class CModel52 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);

	protected:
		float rb,alpha,rf,rs,rc,halfdist,ratio,neinit,wrapcoeff,ratio1,ratio2,stiffness,skinsigmain,skinsigmafr,ldecdens;
};
//! Tube shell model based on CModel45, different cross section radius calculation
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


//! GCS model, as published in ApJs Volume 194, Issue 2, article id. 33, 6 pp. (2011)
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


//! Simple prominence material simulation based on CModel54. The density is placed following the loop axis.
class CModel55 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
    float rb,alpha,rf,ratio,exponent,k2,mk2,nemin,thick,rs,rs2,rc,rc1,rc2;
};


//! Test of a bended streamer slab
class CModel56 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
    float h,s,thetam,theta0,hthick,hang,nefactor;
};


//! Constant and uniform density and temperature
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

//! Blob
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

//! Blob and tail (not sure this works)
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


//! Tube shell model based on CModel54, but the density does not go down to 0 within the shell.
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

/*
* $Log: models51to60.h,v $
* Revision 1.7  2011-01-21 20:45:10  thernis
* Do some more clean up of model 54
*
* Revision 1.6  2011-01-21 15:07:32  thernis
* Clean up Model 54 code
*
* Revision 1.5  2009/02/09 20:51:18  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.4  2007/07/19 19:22:42  thernis
* Add model 60
*
* Revision 1.3  2007/07/10 21:19:01  thernis
* Implement very simple blob comet models. One is not finished...
*
* Revision 1.2  2007/05/14 17:19:41  thernis
* Add CVS id and log in all files
*
*/
