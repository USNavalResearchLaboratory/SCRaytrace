/** \file physicsfcor.h
 * \brief Implements the physics for the F corona
 */

#ifndef PHYSICSFCOR_H
#define PHYSICSFCOR_H

#include "physicsbase.h"
#include "Cvec.h"

/**
F corona physics implementation
*/
class PhysicsFCor : public PhysicsBase
{
public:
    PhysicsFCor() {physicsName="F Corona";};

    ~PhysicsFCor() {};

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout);

    void getConstFactors(float &btf,float &bpf,float &nef, float rho);

};

#endif
