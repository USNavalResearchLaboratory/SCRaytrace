/** \file physicsthomson.h
 * \brief Implement Thomson scattering physics.
 */

#ifndef PHYSICSTHOMSON_H
#define PHYSICSTHOMSON_H

#include "physicsbase.h"
#include "Cvec.h"

/**
Thomson scattering physics implementation
*/
class PhysicsThomson : public PhysicsBase
{
public:
    PhysicsThomson() {physicsName="Thomson Scattering";};

    ~PhysicsThomson() {};

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout);

    void getConstFactors(float &btf,float &bpf,float &nef, float rho);
    
};

#endif
