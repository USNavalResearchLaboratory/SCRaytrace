/** \file physicsisotropic.h
 * \brief Implement a simple isotropic physics. For testing purpose.
 */

#ifndef PHYSICSISOTROPIC_H
#define PHYSICSISOTROPIC_H

#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"

/** \brief Isotropic scattering, for testing purpose
*/
class PhysicsIsotropic : public PhysicsBase
{
public:
    PhysicsIsotropic() {physicsName="Isotropic Scattering";};

    ~PhysicsIsotropic() {};

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout);

    void getConstFactors(float &btf,float &bpf,float &nef, float rho);

    
    
  private:
    
};





#endif
