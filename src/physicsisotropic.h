// $Id$

#ifndef PHYSICSISOTROPIC_H
#define PHYSICSISOTROPIC_H

#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"

/**
Isotropic scattering for testing purpose
*/
class PhysicsIsotropic : public PhysicsBase
{
public:
    PhysicsIsotropic() {physicsName="Isotropic Scattering";};

    ~PhysicsIsotropic() {};

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout);

    void getConstFactors(float &btf,float &bpf,float &nef);

    
    
  private:
    
};





#endif
