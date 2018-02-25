// $Id$

#include "physicsisotropic.h"
#include "rtmiscfunc.h"
#include "Cvec.h"
#include "scene.h"

bool PhysicsIsotropic::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &density)
{
    // -- compute density at point vs
    density = pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition,vs));
    
    // -- compute distance observer - point of LOS
    float x = (vs - pparentscene->obs.o).mag();
    
    // -- Compute distance and angle of vs with respect to observer and Sun
//     float distObs2Sun = (pparentscene->obs.o).mag();
//     float sinTheta = rho / distObs2Sun;
//     float theta = asin(sinTheta);
//     float cosTheta = cos(theta);
//     float phi = atan2(cosTheta - x / distObs2Sun, sinTheta);
    
    // -- compute scattering: phase function is constant
    float integrand = density / (x * x + r * r);

    btout = integrand;
    bpout = integrand;

    return 0;
}

void PhysicsIsotropic::getConstFactors(float &btf,float &bpf,float &nef)
{
    float ds = pparentscene->los.ds / 4 / PI;
    btf=ds;
    bpf=ds;
    nef=ds;
}
