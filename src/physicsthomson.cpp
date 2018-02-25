// $Id: physicsthomson.cpp,v 1.2 2010-09-17 15:22:56 thernis Exp $

#include "physicsthomson.h"
#include "rtmiscfunc.h"
#include "Cvec.h"
#include "scene.h"

bool PhysicsThomson::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout)
{

    neout=pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition,vs));

    if (neout <= MIN_ELECTRON_DENSITY) return 1;

    if (pparentscene->neonly) return 0;

    float btotcoeff,bpolcoeff;
    pparentscene->csun.getThomsonCoeff(r,rho,btotcoeff,bpolcoeff);

    btout=neout*btotcoeff;
    bpout=neout*bpolcoeff;

    return 0;
}

void PhysicsThomson::getConstFactors(float &btf,float &bpf,float &nef)
{
    btf = pparentscene->csun.getConstfactor() * pparentscene->los.ds;
    bpf = btf;
    nef = RSUN_CM * pparentscene->los.ds;
}
