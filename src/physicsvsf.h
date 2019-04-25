/** \file physicsvsf.h
 * \brief Implements the Volume Scattering Function physics, for the modeling of the F corona.
 */

#ifndef PHYSICSVSF_H
#define PHYSICSVSF_H

#include <string>

#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"

/**
Volume Scattering Function, for F-corona and zodiacal light
*/
class PhysicsVSF : public PhysicsBase
{
public:
    PhysicsVSF();

    ~PhysicsVSF();

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout);

    void getConstFactors(float &btf,float &bpf,float &nef, float rho);
    
    
  private:
    static const std::string filename;    //!> filename of the VSF
    float *ang,*vsf;
    int nbsamp;
};





#endif
