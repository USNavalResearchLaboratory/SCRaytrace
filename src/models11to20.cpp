
#include "models11to20.h"
#include <cmath>
#include "Cvec.h"
#include "constant.h"
#include "rtmiscfunc.h"

// -- density 11
// Streamer belt simulation with source surface field map.
float CModel11::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.55) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- param for the c. sheet map
  int sang=360,slat=181; // one point per degree

  // ---- get the distance from the point
  //      to the nearest neighbors
  // -- nearest neighbor seeking range
  int srang=15,srlat=15; // in degrees

  float thetannpos,phinnpos,dist,val;

  int nnok=wherenn(pnsheetmap,sang,slat,srang,srlat,(phi*RADEG),(theta*RADEG),&phinnpos,&thetannpos,&dist,&val);

  if (nnok == 0) return 0.;

  // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};
  float coef[4]={-1,-2,-3,-4};
  float nel=0;
  // -- streamer half thickness
  float w0=0,u0=0; 

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    nel+=c[i]*powr;
    w0+=cw[i]*powr;
    u0+=cu[i]*powr;
    //nel+=c[i]*pow(r,coef[i]);
 }

  float dcross=w0/(2*u0);

  //theta=r*sin(dist*DTOR); // ca c'est tres bourrin !
  theta=dist*DTOR; // ca c'est tres bourrin !

  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }

  //  nel*=val*exp(-ee);

  return nel;
}
void CModel11::initParam(float* pparam)
{
this->pnsheetmap=pparam;
}
void CModel11::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Streamer belt simulation with source surface field map, fixed resolution.","",""));	
	vp.push_back(moddefparam("crot","1912L","Carrington rotation number",""));	
	vp.push_back(moddefparam("modparam","","",""));	
	vp.push_back(moddefparam("","rdtxtmagmap,nsheetmap,crot=$crot","Load PFSS map from WSO web site.",""));	
	vp.push_back(moddefparam("","modparam=reform(nsheetmap,n_elements(nsheetmap))","Format the array of parameters.",""));	
	
	return;
}



// -- density 12
// Simple neutral sheet: sin(theta)
float CModel12::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.05) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- compute distance point-neutral sheet
  float phioffset=10*DTOR;
  float phinbstep=30;
  float phiscanstep=2*phioffset/phinbstep;
    float phiscan=phi-phioffset;
  float mindist2=1e20,dist2,tmp2;
  for (int i=0;i<phinbstep;i++,phiscan+=phiscanstep) {
    tmp2=theta-0.78*sin(phiscan);
    dist2=(phi-phiscan)*(phi-phiscan)+tmp2*tmp2;
    if (dist2 < mindist2) mindist2=dist2;
  }

  // -- streamer half thickness
  float d0=0.1; 
  float ee=fabs(r*sin(sqrt(mindist2))/d0);
  if (ee > 1E1) return 0.;

  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  float coef[4]={-1,-2,-3,-4};
  float nel=0;


  for(int i=0;i<=3;i++) {
    nel+=c[i]*pow(r,coef[i]);
 }


  nel*=exp(-ee);

  return nel;
}
void CModel12::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","Simple neutral sheet: sin(theta)","",""));	
	return;
}



// -- density 13
// Streamer slab model
float CModel13::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.05) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  if (phi < 1.0472 or phi > 2.0944) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
                                                // see 7 May 2004
  // *2 correction fudge factor (don't know why that works !)
  // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
 //  float cw[4]={17.79/4.71,-85.93/4.71,138.32/4.71,391.45/4.71}; // valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float coef[4]={-1,-2,-3,-4};
  float nel=0;

  // -- streamer half thickness
  float w0=0; 

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    nel+=c[i]*powr;
    w0+=cw[i]*powr;
 }

  float ee=(theta*theta)/w0;
  if (ee > 1E1) return 0.;

  nel*=exp(-ee);  

  return nel;
}
void CModel13::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x6;
	vp.push_back(moddefparam("","Sphericaly symmetric streamer, Hayes radial, exp(-d^2/w0(r)) orthoradial","",""));	
	return;
}


// -- density 14
// Streamer slab model
// see 11 May 2004
float CModel14::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  if (phi < 1.0472 or phi > 2.0944) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
                                                // see 7 May 2004
  // *2 correction fudge factor (don't know why that works !)
   // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
//  float cw[4]={17.79/4.71,-85.93/4.71,138.32/4.71,391.45/4.71}; // valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};

  float coef[4]={-1,-2,-3,-4};
  float nel=0;

  // -- streamer half thickness
  float w0=0,u0=0; 
  

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    nel+=c[i]*powr;
    w0+=cw[i]*powr;
    u0+=cu[i]*powr;
 }

  float dcross=w0/(2*u0);

  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }

 
  return nel;
}
void CModel14::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Streamer slab model, Hayes radial, exp(-d^2/w0(r)) and exp(-abs(d/u0(r))) orthoradial.","",""));	
	return;
}


// -- density 15
// Streamer slab model
// see 11 May 2004
float CModel15::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  if ((phi < 1.0472) || (phi > 2.0944)) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
                                                // see 7 May 2004
  // *2 correction fudge factor (don't know why that works !)
   // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
//  float cw[4]={17.79/4.71,-85.93/4.71,138.32/4.71,391.45/4.71}; // valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};

  float coef[4]={-1,-2,-3,-4};
  float nel=0;

  // -- streamer half thickness
  float w0=0,u0=0; 
  

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    nel+=c[i]*powr;
    w0+=cw[i]*powr;
    u0+=cu[i]*powr;
 }

  float dcross=w0/(2*u0);

  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }

  // ---- sinusoidal modulation function
  // the *pmodul should contain 120 sample 
  // representing the angles from 60 to 120 deg.,
  // with a range dynamic from 0 to 1
  // 

  int pos=(int) (-119.+113.637*phi);

  if (pos < 0) pos=0;
  if (pos > 119) pos=119;

  nel*=*(pmodul+pos);

  return nel;
}
void CModel15::initParam(float *pparam)
{
this->pmodul=pparam;
}
void CModel15::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Streamer slab model, Hayes radial, exp(-d^2/w0(r)) and exp(-abs(d/u0(r))) orthoradial, FO modulation.","",""));	
	vp.push_back(moddefparam("pm","","",""));	
	vp.push_back(moddefparam("","pm=(cos(findgen(120)/2)+1.1)/2.1","Example of sinusoidal modulation of the Ne.",""));	
	return;
}


// -- density 16
// Cylinder model
float CModel16::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.5 || z <= 1.1) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  float dist=sqrt(x*x+y*y);

  //if (phi < 1.0472 or phi > 2.0944) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
  // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
  float coef[4]={-1,-2,-3,-4};
  //float gam[3]={16.3*DTOR,10.0*DTOR,43.20*DTOR};
  //float delt[3]={0.5,7.31,7.52};
  //w+=gam[i]*pow(r,-delt[i]);

  float nel=0;

  // -- streamer half thickness
  float d0=0.2; 

  for(int i=0;i<=3;i++) {
    nel+=c[i]*pow(r,coef[i]);
 }

  //nel=1./r;

  float ee=dist/d0;
  ee*=ee;
  if (ee > 1E1) return 0.;

  nel*=exp(-ee);

  return nel;
}
void CModel16::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Cylinder model, Hayes radial, exp(-abs(d/d0)) orthoradial","",""));	
	return;
}


// -- density 17
// Cylinder model with parameters passing
float CModel17::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.5 || z <= 1.1) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  float dist=sqrt(x*x+y*y);

  //if (phi < 1.0472 or phi > 2.0944) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
  float coef[4]={-1,-2,-3,-4};
  //float gam[3]={16.3*DTOR,10.0*DTOR,43.20*DTOR};
  //float delt[3]={0.5,7.31,7.52};
  //w+=gam[i]*pow(r,-delt[i]);

  float nel=0;

  // -- streamer half thickness
  float d0=0.4; 
  nel+=(*pparam)*pow(r,coef[0]);
  nel+=(*(pparam+1))*pow(r,coef[1]);
  nel+=(*(pparam+2))*pow(r,coef[2]);
  nel+=(*(pparam+3))*pow(r,coef[3]);
  
  float ee=dist/d0;
  ee*=ee;
  if (ee > 1E1) return 0.;

  nel*=exp(-ee);

  return nel;
}
void CModel17::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Cylinder model with parameters passing.","",""));	
	vp.push_back(moddefparam("c","[1.34e4,1.15e6,-6.022e6,5.577e7]","Coefficients for the Hayes model.",""));	
	return;
}



// -- density 18
// Streamer slab model with parameter passing.
// see 11 May 2004
// see 14 June 2004
float CModel18::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  if (phi < 1.0472 or phi > 2.0944) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
//  float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
                                                // see 7 May 2004
  // *2 correction fudge factor (don't know why that works !)
 //  float cw[4]={17.79/4.71,-85.93/4.71,138.32/4.71,391.45/4.71}; // valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};

  float coef[4]={-1,-2,-3,-4};
  float nel=0;

  // -- streamer half thickness
  float w0=0,u0=0; 
  
  float *pp;
  pp=pparam;

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    //    nel+=c[i]*powr;
    nel+=(*(pp++))*powr;
    w0+=cw[i]*powr;
    u0+=cu[i]*powr;
 }

  float dcross=w0/(2*u0);

  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }

 
  return nel;
}
void CModel18::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Streamer slab model with parameter passing.","",""));	
	vp.push_back(moddefparam("c","[1.34e4,1.15e6,-6.022e6,5.577e7]","Coefficients for the Hayes model.",""));	
	return;
}



// -- density 19
// Sphericaly symmetric Ne
float CModel19::Density(const Cvec &v) {
  float r=v.norm();
  if (r <= 1.05) return 0.;

  // -- coeff for the equatorial hole density model
  float c1=5.27e6;
  float d1=3.30;
  float c2=3.54e6;
  float d2=5.80;

  float nel=c1*pow(r,-d1)+c2*pow(r,-d2);

  return nel;
}
void CModel19::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Saito equatorial hole density.","",""));
	return;
}


// -- density 20
// Sphericaly symmetric Ne
float CModel20::Density(const Cvec &v) {
  float r=v.norm();
  if (r <= 1.05) return 0.;

  // -- coeff for the polar regions (hole) density model
  float c1=3.15e6;
  float d1=4.71;
  float c2=1.60e6;
  float d2=3.01;

  float nel=c1*pow(r,-d1)+c2*pow(r,-d2);

  return nel;
}
void CModel20::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Saito polar regions (hole) density","",""));	
	return;
}

