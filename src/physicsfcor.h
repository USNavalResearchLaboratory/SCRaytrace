// $Id: physicsfcor.h,v 1.1 2010-09-01 19:54:30 thernis Exp $


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
