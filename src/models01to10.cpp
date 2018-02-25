/***************************************************************************
 *  $Id: models01to10.cpp,v 1.3 2009/02/09 20:51:12 thernis Exp $
 *
 ****************************************************************************/

#include "models01to10.h"
#include <cmath>
#include "constant.h"
#include "rtmiscfunc.h"



// -- density 1
// M.Guhathakurta model
float CModel01::Density(const Cvec &v,float &temperature) 
{
	float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
	
	if (r <= 1.05) return 0.;

	float alpha=0.*DTOR;
	float phi0=90.*DTOR;
//	float nch=1e10;
	
	float phi=0,theta;

	float cp[3]={0.14e7,8.02e7,8.12e7};
	float dp[3]={2.8,8.45,16.87};
	float ccs[3]={1.41e7,16.42e7,61.90e7};
	float dcs[3]={3.39,5.14,12.67};
	float gam[3]={16.3,10.0,43.20};
	float delt[3]={0.5,7.31,7.52};
	
	float np=0,ncs=0,w=0;
	float thetamg,nel;
	
	theta=asin(z/r);
	float ratio=x/(r*cos(theta));

	if (fabs(ratio) < 1) {
	  phi=acos(ratio);
	  if (y < 0) phi=2.*PI-phi;
	}

	
	for(int i=0;i<=2;i++) {
		np+=cp[i]*pow(r,-dp[i]);
		ncs+=ccs[i]*pow(r,-dcs[i]);
		w+=gam[i]*pow(r,-delt[i]);
	}

	thetamg=asin(-cos(theta)*sin(alpha)*sin(phi-phi0)+sin(theta)*cos(alpha));

	float expo=pow(thetamg/(w*DTOR),2);
	if (expo > 1E1) nel=np; else nel=np+(ncs-np)*exp(-expo);

	return nel;
}
void CModel01::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","M.Guhathakurta model, frozen parameters.","",""));	
	return;
}




// -- density 2
// Model for testing purpose
float CModel02::Density(const Cvec &v,float &temperature) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.05) return 0.;

  float nel=0;

//  float rhosquare=x*x+y*y;

  if (z > 0.) {
    float powexp=(x*x+y*y)/(0.5*0.5);
    if (powexp < 1e1) nel=exp(-r/1.)*exp(-powexp);
  }
  return nel;
}
void CModel02::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","Model for testing purpose","",""));	
	return;
}



// -- density 3
// Model for testing purpose
float CModel03::Density(const Cvec &v,float &temperature) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.05) return 0.;
	
  float nel=0;

  //  float theta=0,phi;  
  //phi=acos(z/r);
  //float ratio=x/(r*sin(phi));
  //if (abs(ratio) < 1.) {
  //  theta=acos(ratio);
  //  if (y < 0.) theta=2.*PI-theta;
  // }

  float phi,theta;

  cvcoord(x,y,z,&r,&theta,&phi);

  //nel=(1./r)*exp(-pow(2.*sin(4*theta),2)/2.)*exp(-pow(phi-PI/2,2)/0.1);
  if (r > 3 && r < 6 && phi < 0.5 && phi > -0.5 && theta >= 0 && theta <= PI/2) nel=1e4;

  //nel=(1./r)*exp(-pow(phi-PI/2+0.05*cos(theta),2)/0.01);

  return nel;
}
void CModel03::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","Model for testing purpose.","",""));	
	return;
}

// -- density 4
// Density cube from Y.M. with trilinear interpolation
float CModel04::Density(const Cvec &v,float &temperature) 
{
return(densyming2(v.v[0],v.v[1],v.v[2],0.,
		pmparam,pmparam+3,pmparam+6,pmparam+9));
}
void CModel04::initParam(float* pparam)
{
this->pmparam=pparam;
}
void CModel04::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Density cube from Y-M.Wang, trilinear interpolation in the cube.","",""));	
	vp.push_back(moddefparam("mp","","",""));	
	vp.push_back(moddefparam("","load_dens,dens,rco,phico,thetaco","Load PFSS default density cube.",""));	
	vp.push_back(moddefparam("","mp=[rco,phico,thetaco,reform(dens,n_elements(dens))]","Format the array of parameters.",""));	
	vp.push_back(moddefparam("","dens=0b","Free memory.",""));	

	return;
}

// -- density 5
//! Same as Model 4 but without trilinear interpolation
//!
float CModel05::Density(const Cvec &v,float &temperature) {
   return(densymingnn(v.v[0],v.v[1],v.v[2],0.,
		      pmparam,pmparam+3,pmparam+6,pmparam+9));
 }
void CModel05::initParam(float* pparam)
{
this->pmparam=pparam;
}
void CModel05::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Density cube from Y-M.Wang, no interpolation in the cube.","",""));	
	vp.push_back(moddefparam("modparam","","",""));	
	vp.push_back(moddefparam("","load_dens,dens,rco,phico,thetaco","Load PFSS default density cube.",""));	
	vp.push_back(moddefparam("","modparam=[rco,phico,thetaco,reform(dens,n_elements(dens))]","Format the array of parameters.",""));	
	vp.push_back(moddefparam("","dens=0b","Free memory.",""));	

	return;
}


// -- density 6
// Very simple model of streamer belt
float CModel06::Density(const Cvec &v,float &temperature) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.05) return 0.;
	
  float d0=0.3; // half thickness of the streamer

  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&theta,&phi);

  float neradial=1/r;

  float phi0=0.5*cos(2*theta);
  float d=r*sin(phi-phi0);
  float neshape=exp(-pow(d/d0,4));

  float nel=neradial*neshape;

  return nel;
}
void CModel06::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","","",""));	
	return;
}

// -- density 7
// Slab model of streamer 
// based on Guhathakurta model.
//
float CModel07::Density(const Cvec &v,float &temperature) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.05) return 0.;

  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  if (phi < 0 or phi > (180.*DTOR)) return 0.;
  //  if (abs(theta) > (60./DTOR)) return 0.;


//  float cp[3]={0.14e7,8.02e7,8.12e7};
//  float dp[3]={2.8,8.45,16.87};
  float ccs[3]={1.41e7,16.42e7,61.90e7};
  float dcs[3]={3.39,5.14,12.67};
  float gam[3]={16.3*DTOR,10.0*DTOR,43.20*DTOR};
  float delt[3]={0.5,7.31,7.52};

  float ncs=0,w=0;//np=0

  for(int i=0;i<=2;i++) {
    //np+=cp[i]*pow(r,-dp[i]);
    ncs+=ccs[i]*pow(r,-dcs[i]);
    w+=gam[i]*pow(r,-delt[i]);
  }

  // use only the curent sheet model
  float ee=pow(theta/w,2);

  if (ee > 1E1) return 0.;

  float nel=ncs*exp(-ee);
  //float nel=np;

  return nel;
}
void CModel07::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","Slab model of streamer, based on Guhathakurta model.","",""));	
	return;
}




// -- density 8
// Sphericaly symmetric Ne
float CModel08::Density(const Cvec &v,float &temperature) {
//  float x=v.v[0],y=v.v[1],z=v.v[2];
  float r=v.norm();
  if (r <= 1.05) return 0.;

  
  float c[4]={-70544.822   ,    1206889.8  ,    -4264689.2    ,   9676682.1};
  float coef[4]={-1,-2,-3,-4};
  float nel=0;

  for(int i=0;i<=3;i++) {
    nel+=c[i]*pow(r,coef[i]);
 }

  return nel;
}
void CModel08::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Uniform corona, Hayes radial.","",""));	
	return;
}


// -- density 9
// Guhata. radial and Vibert shape
float CModel09::Density(const Cvec &v,float &temperature) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.05) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  if (phi < 0 or phi > (180.*DTOR)) return 0.;
  //  if (abs(theta) > (60./DTOR)) return 0.;


//  float cp[3]={0.14e7,8.02e7,8.12e7};
//  float dp[3]={2.8,8.45,16.87};
  float ccs[3]={1.41e7,16.42e7,61.90e7};
  float dcs[3]={3.39,5.14,12.67};
//  float gam[3]={16.3*DTOR,10.0*DTOR,43.20*DTOR};
//  float delt[3]={0.5,7.31,7.52};

  // -- streamer half thickness
  float d0=0.2; 

//  float np=0,ncs=0,w=0;
  float ncs=0;
  // -- radial term
  for(int i=0;i<=2;i++) {
    //np+=cp[i]*pow(r,-dp[i]);
    ncs+=ccs[i]*pow(r,-dcs[i]);
    //w+=gam[i]*pow(r,-delt[i]);
  }

  // -- shape term
  //float ee=pow((z/d0),4);
  float ee=fabs(z/d0);
  if (ee > 1E1) return 0.;

  float nel=ncs*exp(-ee);

  return nel;
  
}
void CModel09::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","Guhata. radial and Vibert shape","",""));	
	return;
}




// -- density 10
// Sphericaly symmetric Streamer
float CModel10::Density(const Cvec &v,float &temperature) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.5) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  if (phi < 1.0472 or phi > 2.0944) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //  float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
  // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun

  float coef[4]={-1,-2,-3,-4};
  //float gam[3]={16.3*DTOR,10.0*DTOR,43.20*DTOR};
  //float delt[3]={0.5,7.31,7.52};
  //w+=gam[i]*pow(r,-delt[i]);

  float nel=0;

  // -- streamer half thickness
  float d0=0.5; 

  for(int i=0;i<=3;i++) {
    nel+=c[i]*pow(r,coef[i]);
  }
  
  //nel=1./r;

  float ee=z/d0;
  ee*=ee;
  if (ee > 1E1) return 0.;

  nel*=exp(-ee);

  return nel;
}
void CModel10::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","Sphericaly symmetric streamer, Hayes radial, exp(-abs(d/d0)) orthoradial","",""));	
	return;
}



/*
* $Log: models01to10.cpp,v $
* Revision 1.3  2009/02/09 20:51:12  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.2  2007/05/14 17:19:41  thernis
* Add CVS id and log in all files
*
*/
