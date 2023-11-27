/** \file physicsvariablevsf.h
 * \brief Implements the a Volume Scattering Function that can be modified by changing the .dat file specified by filename.
 */

#ifndef PHYSICSVARIABLEVSF_H
#define PHYSICSVARIABLEVSF_H

#include <string>

#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"

/** \brief Volume Scattering Function, for F-corona and zodiacal light, that can be modified as previously mentioned.

Caveat: The VSF is valid down to 5 deg elongation from the Sun center.

*/
class PhysicsVariableVSF : public PhysicsBase
{
public:
    PhysicsVariableVSF();

    ~PhysicsVariableVSF();

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
