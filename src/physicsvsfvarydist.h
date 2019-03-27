// $Id$

#ifndef PHYSICSVSFVARYDIST_H
#define PHYSICSVSFVARYDIST_H

#include <string>

#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"

/**
Volume Scattering Function, for F-corona and zodiacal light, varying with distance as described in Lamy & Perrin 1986.

Caveat: The VSF is valid down to 5 deg elongation from the Sun center.

*/
class PhysicsVSFVaryDist : public PhysicsBase
{
public:
    PhysicsVSFVaryDist();

    ~PhysicsVSFVaryDist();

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout);

    
    //! Returns the constant factor of the integration
    void getConstFactors(float &btf,float &bpf,float &nef, float rho);
    
    
    void initDensityModel(const Cvec &obsPos);

        
  private:
    static const std::string filename;    //!> filename of the VSF
    float *ang, *vsf;
    int nbsamp;
};





#endif
