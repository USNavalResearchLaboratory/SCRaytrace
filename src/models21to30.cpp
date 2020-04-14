/***************************************************************************
 *  $Id: models21to30.cpp,v 1.8 2009/03/17 14:44:43 thernis Exp $
 *
 ****************************************************************************/

#include "models21to30.h"
#include <cmath>
#include <vector>
#include "Cvec.h"
#include "constant.h"
#include "rtmiscfunc.h"

// -- density 21
// Sphericaly symmetric Ne
//
// Use the Saito background (equator) density model
float CModel21::Density(const Cvec &v) {
  float r=v.norm();
  if (r <= 1.05) return 0.;

  // -- coeff for the background (equator) density model
  float c1=1.36e6;
  float d1=2.14;
  float c2=1.68e8;
  float d2=6.13;

  float nel=c1*pow(r,-d1)+c2*pow(r,-d2);

  return nel;
}
void CModel21::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Saito background (equator) density.","",""));
	return;
}

// -- density 22
// Same as density 14 but with a modified value of the angular sector
// see 11 May 2004
float CModel22::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  //if (phi < 1.0472 or phi > 2.0944) return 0.;
  //  if (phi < (45*DTOR) or phi > (135*DTOR)) return 0.;
  if (phi < (45*DTOR) or phi > (135*DTOR)) return 0.;
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
void CModel22::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Streamer slab model, with a modified value of the angular sector compare to model 14.","",""));	
	return;
}



// -- density 23
// Same as density 18 but with a modified value of the angular sector
// see 11 May 2004
// see 14 June 2004
float CModel23::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  float sechang=*(pparam+4);

  if (phi < (PI/2-sechang) or phi > (PI/2+sechang)) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
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
void CModel23::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel23::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Streamer slab model, with a modified value of the angular sector compare to model 18.","",""));	
	vp.push_back(moddefparam("c","[1.34e4,1.15e6,-6.022e6,5.577e7,0.698]","Coefficients for the Hayes model.",""));	
	return;
}



// -- density 24
// Cylinder model, on axis Ox
float CModel24::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.5 || x <= 1.1) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  float dist=sqrt(z*z+y*y);

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
void CModel24::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0x4;
	vp.push_back(moddefparam("","Cylinder model, Hayes radial, exp(-abs(d/d0)) orthoradial, axis Ox","",""));	
	return;
}

// -- density 25
// Datacube passed in pparam pointer
float CModel25::Density(const Cvec &v) {
  // -- nothing inside the Sun !
  if (v.norm() <= 1.01) return 0.;

  // ---- compute position of the voxel in the data cube
  //      return 0 if out of the cube
	// -- X
	float xf=(pparam[3]+v.v[0]/voxsize);
  if (xf < 0 || xf > (pparam[0])) return 0;
  long xp=long(floor(double(xf)));

  // -- Y
  float yf=(pparam[4]+v.v[1]/voxsize);
  if (yf < 0 || yf > (pparam[1])) return 0;
  long yp=long(floor(double(yf)));

  // -- Z
  float zf=(pparam[5]+v.v[2]/voxsize);
  if (zf < 0 || zf > (pparam[2])) return 0;
  long zp=long(floor(double(zf)));

  // -- position of the voxel:
  long offset=7 + xp + yp*sx + zp*sx*sy;

  return pparam[offset];
}
void CModel25::initParam(float *pparam)
{

this->pparam=pparam;
// -- get the voxel side size in rsun
voxsize=pparam[6];
sx=long(pparam[0]);
sy=long(pparam[1]);

}
void CModel25::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Density cube.","",""));	
	vp.push_back(moddefparam("cubefile","getenv('RT_PATH')+get_delim()+'testcube.fts'","Filename of the density cube",""));	
	vp.push_back(moddefparam("mp","","",""));	
	vp.push_back(moddefparam("","mp=loaddenscube($cubefile)","Load the density cube.",""));	

	return;
}


// -- density 26
// Datacube passed in pparam pointer.
float CModel26::Density(const Cvec &vv) {

  // -- compute position of the voxel in the data cube
  //    return 0 if out of the cube
  // ---- X
  float xf=( (vv[0]/voxsize) + pparam[3]);
  if (xf < 0 || xf >=pparam[0]) return 0;
  float xs,t;
	long xp;
	if (xf <= 0.5) 
  {
    xs=0.;xp=0;t=0.;
  } 
  else if (xf >= (pparam[0]-0.5)) 
  {
    xs=pparam[0]-2;xp=long(xs);t=1.;
  } 
  else 
  {
    xs=xf-0.5;xp=long(xs);t=(xs-float(xp));
  }

  // ---- Y
  float yf=((vv[1]/voxsize) + pparam[4]);
  if (yf < 0 || yf >=pparam[1]) return 0;
  float ys,u;
	long yp;
	if (yf <= 0.5) 
  {
    ys=0.;yp=0;u=0.;
  } 
  else if (yf >= (pparam[1]-0.5)) 
  {
    ys=pparam[1]-2;yp=long(ys);u=1.;
  } 
  else 
  {
    ys=yf-0.5;yp=long(ys);u=(ys-float(yp));
  }

  // ---- Z
  float zf=((vv[2]/voxsize) + pparam[5]);
  if (zf < 0 || zf >=pparam[2]) return 0;
  float zs,v;
	long zp;
	if (zf <= 0.5) 
  {
    zs=0.;zp=0;v=0.;
  } 
  else if (zf >= (pparam[2]-0.5)) 
  {
    zs=pparam[2]-2;zp=long(zs);v=1.;
  } 
  else 
  {
    zs=zf-0.5;zp=long(zs);v=(zs-float(zp));
  }

	float *pcube;
	pcube=pparam+7;

  // -- trilinear interpolation
  float neinterp=trilininterp(t,u,v,xp,yp,zp,pcube,sx,sy);

  return neinterp;
}
void CModel26::initParam(float *pparam)
{
this->pparam=pparam;

voxsize=pparam[6];
sx=long(pparam[0]);
sy=long(pparam[1]);
sxsy=sx*sy;

}
void CModel26::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Density cube, trilinear interpolation for smoothing.","",""));	
	vp.push_back(moddefparam("cubefile","getenv('RT_PATH')+get_delim()+'testcube.fts'","Filename of the density cube",""));	
	vp.push_back(moddefparam("mp","","",""));	
	vp.push_back(moddefparam("","mp=loaddenscube($cubefile)","Load the density cube.",""));	
	
	return;
}

//------------------------------------------------
// -- density 27
// Spherical Shell Ne density
float CModel27::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2];//,r=v.norm();

  float r0 = pparam[0]; // 8; //height of center sphere
  float d = pparam[1];  // 2; //radius of sphere
  float d0 = pparam[2]; //.3; //1/2 thickness of shell
  float nemin = pparam[3]; // 1e5
  //float nemax = 1e6;

  float rr = sqrt(x*x+y*y+(z-r0)*(z-r0)); //-r0

  if (z < r0) return 0.;
    
  //float nel=0;
  if (fabs(rr-d) <= d0){
    return nemin;// + (nemax-nemin)*(z-r0-d-d0)/(2*d+2*d0);
  }

  return 0. ;
}
void CModel27::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel27::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Spherical Shell.","",""));	
	vp.push_back(moddefparam("r0","8.","height of center sphere","Rsun"));	
	vp.push_back(moddefparam("d","2.","radius of sphere","Rsun"));	
	vp.push_back(moddefparam("d0","0.3","1/2 thickness of shell","Rsun"));	
	vp.push_back(moddefparam("nemin","1e5","Ne","cm^-3"));	

	return;
}




//-----------------------------------------------
// -- density 28
// Cylindrical Shell Ne density
float CModel28::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2];//,r=v.norm();

  float r0  = pparam[0]; //8;   //height of center sphere
  float d   = pparam[1]; //2;   //radius of sphere
  float len = pparam[2]; //2;   //length of cylender in each direction
  float d0  = pparam[3]; //0.3; //1/2 thickness of shell

  float nemin = 1;
  float nemax = 1;

  float rr = sqrt(x*x+(z-r0)*(z-r0));

  float nel=0;

  if (fabs(y) <= len)
    if (fabs(rr-d) <= d0) {
      nel = nemin + (nemax-nemin) * (z-r0-d-d0)/(2*d+2*d0);
      //ofstream myFile("output22.txt",ios::app);
      //myFile<<x<<" "<<y<<" "<<z<<" "<<nel<<endl;  
      //myFile.close();
    }

  return nel;
}
void CModel28::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel28::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Cylindrical Shell.","",""));	
	vp.push_back(moddefparam("r0","8.","height of axis","Rsun"));	
	vp.push_back(moddefparam("d","2.","radius of shell","Rsun"));	
	vp.push_back(moddefparam("len","2.","Lenght of the cylinder","Rsun"));	
	vp.push_back(moddefparam("d0","0.3","Half thickness of the shell","Rsun"));	

	return;
}

//-----------------------------------------------
// -- density  29
// Bow Shock Model by Angelos Vourlidas
float CModel29::Density(const Cvec &v) {

  //pmparam=d,s,h,d0,nel0

  float x=v.v[0],y=v.v[1],z=v.v[2];//,r=v.norm();

  float d = pparam[0]; //1;
  float s = pparam[1]; //1.5;
  float h = pparam[2]; //1;?
  float d0 = pparam[3]; //.1
  float nel0 =pparam[4];//1e5
  //s = 1/h;
  //s/=2;
  //s+=2;
  //d+=h/5;

  float zz = h - (d/s)*pow((sqrt(x*x+y*y)/d),s);
  
 float nel = 0;

  
  if (fabs(zz-z) <= d0) {
    nel = nel0;
  }

  return nel;
}
void CModel29::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel29::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Bow Shock Model (A.Vourlidas)","",""));	
	vp.push_back(moddefparam("d","1.","","Rsun"));	
	vp.push_back(moddefparam("s","1.5","","Rsun"));	
	vp.push_back(moddefparam("h","1","Height of the shell","Rsun"));	
	vp.push_back(moddefparam("d0","0.1","Thickness on the shell","Rsun"));	
        vp.push_back(moddefparam("nel0","1e5","density","cm^-3"));
	return;
}


// -- density 30
// Streamer belt simulation with source surface field map:
float CModel30::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.55) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- param for the neutral sheet map
  int sang=(int) *pparam,slat=(int) *(pparam+1);

  // ---- get the distance from the point
  //      to the nearest neighbors
  // -- nearest neighbor seeking range
  int srang=15,srlat=15; // in pix

  // -- point to the neutral sheet map
  float *pnsheetmap;
  pnsheetmap=pparam+2;

  float thetannpos,phinnpos,dist,val;

  int nnok=wherennsmoothed(pnsheetmap,sang,slat,srang,srlat,(phi*RADEG),(theta*RADEG),&phinnpos,&thetannpos,&dist,&val);

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
  //float ee=(theta*theta)/0.005;
  //float ee=theta/0.01;
  //if (ee > 1E1) return 0;
  //nel*=exp(-ee);  

  
  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=val*exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=val*exp(-ee);

  }
  
  //nel*=val*exp(-ee);

  return nel;
}
void CModel30::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel30::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Streamer belt simulation with source surface field map, user resolution.","",""));	
	vp.push_back(moddefparam("crot","1912L","Carrington rotation number",""));	
	vp.push_back(moddefparam("modparam","","",""));	
	vp.push_back(moddefparam("","rdtxtmagmap,nsheetmap,crot=$crot","Load PFSS map from WSO web site.",""));	
	vp.push_back(moddefparam("","modparam=reform(360,181,nsheetmap,n_elements(nsheetmap))","Format the array of parameters.",""));	
	
	return;
}


/*
* $Log: models21to30.cpp,v $
* Revision 1.8  2009/03/17 14:44:43  thernis
* Fix bugs in model 26
*
* Revision 1.7  2009/02/09 20:51:14  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.6  2008/09/23 14:08:12  thernis
* - add new testing models, - implement integration in front of the instrument
*
* Revision 1.5  2008/08/18 18:16:15  thernis
* Change comments in model 25
*
* Revision 1.4  2008/08/18 17:22:11  ontivero
* update model 29
*
* Revision 1.3  2007/11/20 20:46:32  thernis
* Fix rounding problems in model 25
*
* Revision 1.2  2007/05/14 17:19:41  thernis
* Add CVS id and log in all files
*
*/
