//
// C++ Interface: rtmiscfunc
//
// Description: Miscellaneous functions useful for the raytracing and the models
//
// $Id: rtmiscfunc.h,v 1.3 2010-09-01 15:34:45 thernis Exp $
//

#ifndef RTMISCFUNC_H
#define RTMISCFUNC_H


#include <iostream>
#include "constant.h"
#include "Cvec.h"
#include "Cbasis.h"
#include "ModelPosition.h"


using namespace std;




//! Search index if a value in a 1D array:
//! The array must be sorted in a increasing order
//!< @param xin value to search in the array
//!< @param xvec the array to search into, sorted in increasing order
//!< @param nbsamp the size of the array
//!< @param flagout returns -1 if xin si too low, +1 if too high, 0 if within range
//!< @param idshift idshift between 0 and 1, gives the relative position of xin compare to the contiguous array elements where it's been found
//!< @return the position index of the array element just lt xin
inline unsigned int searchnearestindex02(const float &xin,const float *xvec,const unsigned int &nbsamp,short &flagout,float &idshift)
{
// ---- check if too low
if (xvec[0]>xin)
  {
  flagout=-1;
  return 0;
  }
// ---- check if too high
if (xvec[nbsamp-1]<xin)
  {
  flagout=1;
  return 0;
  }

flagout=0;

unsigned int i=0;
if (xvec[i]!=xin) {
	while ((xvec[i]<=xin) && (i < (nbsamp))) i++;
	i--;
}

idshift=0.;
idshift=(xin-xvec[i])/(xvec[i+1]-xvec[i]);

return i;
}





//! Search index if a value in a 1D array:
//! The array must be sorted in a increasing order
inline unsigned int searchnearestindex(const float &xin,float *xvec,const unsigned int &nbsamp)
{
unsigned int i=0;
while ((xvec[i]<xin) && (i < (nbsamp))) i++;
return i;
}


//! Nearest neighbor 1D interp
inline float nearestneighbor1dinterp(const float &xin,float *xvec,const unsigned int &nbsamp,float *yvec)
{
unsigned int i=searchnearestindex(xin,xvec,nbsamp);
float yout=0.;
if (i>0 && i<nbsamp) {
	yout=yvec[i-1]+(xin-xvec[i-1])*(yvec[i]-yvec[i-1])/(xvec[i]-xvec[i-1]);
  } else if (i==0) {
    yout=yvec[0];
  } else {
    yout=yvec[nbsamp-1];
  }
  return yout;
}







//! Change cartesian coordinates to polar coordinates
inline void cvcoord ( const float &x,const float &y,
                const float &z,
                float *pr,float *pphi,float *ptheta )
{
    *pr=sqrt ( x*x+y*y+z*z );
    *ptheta=asin ( z/ *pr );
    float ratio=fabs ( x ) / ( *pr * cos ( *ptheta ) );
    *pphi=0;

    if ( fabs ( ratio ) < 1 )
    {
        *pphi=acos ( ratio );
    }
    if ( x < 0 )
        *pphi=PI-*pphi;
    if ( y < 0 )
        *pphi=TWOPI-*pphi;

    return;
}

//! Change cartesian coordinates to polar coordinates
inline void cvcoord ( const Cvec &p,float *pr,float *pphi,float *ptheta )
{
    return cvcoord ( p.v[0],p.v[1],p.v[2],pr,pphi,ptheta );
}

//! Change cartesian coordinates to polar coordinates
inline void cvcoord ( const float &x,const float &y,
                float *pr,float *pphi )
{
    *pr=sqrt ( x*x+y*y );
    float ratio=fabs ( x ) / ( *pr );
    *pphi=0;

    if ( fabs ( ratio ) < 1 )
    {
        *pphi=acos ( ratio );
    }
    if ( x < 0 )
        *pphi=PI-*pphi;
    if ( y < 0 )
        *pphi=TWOPI-*pphi;

    return;
}




//! Convert from Carrington to Cartesian coordinates
inline Cvec carrington2cart(const float &lon,const float &lat,const float &height)
{
  Cvec c;
  
  c.v[0]=height*sin(lat);
  c.v[1]=height*sin(lon)*cos(lat);
  c.v[2]=height*cos(lon)*cos(lat);
  
  return c;
}

//! Convert from Cartesian to Carrington coordinates
inline void cart2carrington(const Cvec &c,float *lonlatcar)
{
  
  lonlatcar[2]=c.mag();
  lonlatcar[1]=asin(c.v[0]/lonlatcar[2]);
  lonlatcar[0]=atan2(c.v[1],c.v[2]);
  
}

//! Convert from Cartesian to r lon lat: Ox: rotation axis of the Sun, Oz: points to (lon,lat)=(0,0)
inline void cart2rlonlat(const Cvec &v,float &r,float &lon,float &lat)
{
r=v.mag();
lon=atan2(-v[1],v[2]);
while (lon < 0.) lon+=TWOPI;
lat=asin(v[0]/r);
}






//! Apply WCS projection to the rrr angle
template< typename T > 
inline T applyprojection(const T &rrr,const int &projtypecode,const T &pv2_1) {
  T rout=rrr;
  switch (projtypecode) {
    case 1 : break; // ARC 
    case 2 : {rout=tan(rrr);break;} // TAN
    case 3 : {rout=sin(rrr);break;} // SIN
    case 4 : rout=((pv2_1+1)*sin(rrr))/(pv2_1+cos(rrr)); // AZP
  }
  return rout;
}

//! Apply inverse WCS projection to the rrr angle
inline float applyinverseprojection(const float &rrr,const int &projtypecode,const float &pv2_1) {
  float rout=rrr;
  
  switch (projtypecode) {
    case 1 : break; // ARC 
    case 2 : {rout=atan(rrr);break;} // TAN
    case 3 : {rout=asin(rrr);break;} // SIN
    case 4 : {
      float r=rrr/(pv2_1+1);
      rout=(PI/2.)-(atan2(float(1.),float(r))-asin(r*pv2_1/sqrt(r*r+1.))); // AZP
    }
  }
  return rout;
}








//! Trilinear inperpolation in a density cube
float trilininterp(float x,float y,float z,
		   float xint[2],float yint[2],float zint[2],
		   float cube[2][2][2]);
float trilininterp(const float &t,const float &u,const float &v,
										const unsigned int &xi0,
										const unsigned int &yi0,
										const unsigned int &zi0,
										float *cube,
										const unsigned int &sx,const unsigned int &sy);


//! Seeks the nearest neighbor points in a neutral sheet map 
int wherenn(float* pnsheetmap,int slon,int slat,
	int srlon,int srlat,
	float lon,float lat,
	float* plon,float* plat,float* pdist,float* pval);


//! Seeks the nearest neighbor points in a neutral sheet map with a smoothing
int wherennsmoothed(float* pnsheetmap,int slon,int slat,
	int srlon,int srlat,
	float lon,float lat,
	float* plon,float* plat,float* pdist,float* pval);


//! Return the pixel value at the requested position on a PFSS map
int getPosOnSSMap(float *pnsheetmap,int sang,int slat,
	float phi,float theta,float& dist,float& val);

//! Get the density out of a cube in polar coord and with trilinear interp.
float densyming2(const float x,const float y,const float z,const float phi0,
		 const float *prco,const float *pphico,const float *pthetaco,
		 float *pdens);

//! Get the density out of a cube in polar coord, use nearest neighbor.
float densymingnn(const float x,const float y,const float z,const float phi0,
		 const float *prco,const float *pphico,const float *pthetaco,
		  float *pdens);



//! Change coordinate to density
inline Cvec ChangetoDensityCoord(const Cbasis &nps,const Cvec &vs)
{
return nps.u * (vs - nps.o);
}

//! Change coordinate to density
inline Cvec ChangetoDensityCoord(const Cbasis &nps,const Cvec &vs,const Cbasis &necoord2)
{
return (necoord2.u * (ChangetoDensityCoord(nps,vs) - necoord2.o) + necoord2.o);
}

//! Change coordinate to density
inline Cvec ChangetoDensityCoord(const ModelPosition &m,const Cvec &vs)
{
return (m.rotation.u * (ChangetoDensityCoord(m.modelbasis,vs) - m.rotation.o) + m.rotation.o - m.translation);
}


//! Change coordinate from density to abs
inline Cvec ChangeDensityCoordtoAbs(const ModelPosition &m,const Cvec &v)
{
return (m.modelbasis.ui*((m.rotation.ui*(v+m.translation-m.rotation.o))+m.rotation.o)+m.modelbasis.o);
}




#endif



// $Log: rtmiscfunc.h,v $
// Revision 1.3  2010-09-01 15:34:45  thernis
// Clean up
//
// Revision 1.2  2009/04/13 21:02:57  thernis
// - Implement ChangeCoordtoDensity function.
//
// Revision 1.1  2009/02/09 20:46:32  thernis
// - Clean up the code
// - Change CModel::Density prototype
// - Update documentation
// - Implement multi-threading using boost thread
// - Add new models
//
