// $Id: sun.h,v 1.2 2010-09-17 15:15:11 thernis Exp $

#ifndef SUN_H
#define SUN_H

#include "Cmat.h"
#include "Cvec.h"
#include "Cbasis.h"

class Scene;

//!Parameters of the Sun and Thomson scattering computation function
class Sun{
public:
    Sun();
    ~Sun();
private:
	void calcConstfactor() {constfactor=1.24878E-25/(1-u/3)*RSUN_CM;}
public:
	float getConstfactor() {return constfactor;}

//! Compute the Thomson scattering coefficients for a given position in space
//! \param r distance to the center of the Sun
//! \param rho impact parameter
//! \param btotcoeff returns the coefficient for the total brightness
//! \param bpolcoeff returns the coefficient for the polarized brightness
inline void getThomsonCoeff(const float &r,const float &rho,float &btotcoeff,float &bpolcoeff)
{
  double sinomega=1./double(r);
  double sinsquareomega=sinomega*sinomega;
  double cossquareomega=1-sinsquareomega;
  double cosomega=sqrt(cossquareomega);
    
  double logterm=log((1.+sinomega)/cosomega);
    
  double a=cosomega*sinsquareomega;
  double b=-1./8.*(1.-3.*sinsquareomega-cossquareomega*((1.+3.*sinsquareomega)/sinomega)*logterm);
    
  double c=(4./3.)-cosomega-(cosomega*cossquareomega)/3.;
  double d=(1./8.)*(5.+sinsquareomega-cossquareomega*((5.-sinsquareomega)/sinomega)*logterm);
    
    // ---- sum in the pixel
  double rhooverr=rho/r;
    // the polarized brightness
  bpolcoeff=float((a+u*(b-a))*rhooverr*rhooverr);
  btotcoeff=float((2*(c+u*(d-c))-bpolcoeff));

};
    float getRadius();
    float getLimbDarkening();
    void setLimbDarkening(float limbdarkeningcoeff=0.7);
    Cbasis* getPosition();
    
	friend class Scene;
		
protected:
  double u; //!< limb darkening
  float radius; //!< radius of the Sun
  float constfactor; //!< constant multiplicative factor of the Thomson scattering integral
  Cbasis *basis;
};



#endif
