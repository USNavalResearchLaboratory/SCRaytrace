
#include "physicsthomson.h"
// #include "rtmiscfunc.h"
#include "Cvec.h"
#include "scene.h"

bool PhysicsThomson::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout)
{

    neout = pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition, vs));

    if (neout <= MIN_ELECTRON_DENSITY) return 1;

    if (pparentscene->neonly) return 0;

    float btotcoeff,bpolcoeff;
    this->getThomsonCoeff(r,rho,btotcoeff,bpolcoeff);

    btout = neout * btotcoeff;
    bpout = neout * bpolcoeff;

    return 0;
}

void PhysicsThomson::getConstFactors(float &btf,float &bpf,float &nef, float rho)
{
    

    
    btf = constfactor * pparentscene->los.ds;
    bpf = btf;
    nef = RSUN_CM * pparentscene->los.ds;
}




