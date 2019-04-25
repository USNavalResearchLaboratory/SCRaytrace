
#include "physicsmie.h"
#include "rtmiscfunc.h"
#include "Cvec.h"
#include "scene.h"


bool PhysicsMie::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btot,float &bpol,float &density)
{
    // ---- compute density
    density = pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition,vs));

    // -- compute distance observer - point of LOS
    Cvec vx = vs - pparentscene->obs.o;
    float x = vx.mag();
    
    // ---- compute scattering angle
    float cosAlpha = pscal(vs,vx) / (r * x);
    if (cosAlpha > 1.) cosAlpha = 1.;
    if (cosAlpha < -1.) cosAlpha = -1.;
    float phi = acos(-cosAlpha);
    
    unsigned int idx;
    idx = (unsigned int) (float(SIZELOOKUPTABLE-1) * phi / PI);
    
    if (idx > (SIZELOOKUPTABLE-1)) idx = SIZELOOKUPTABLE-1;
    
    // -- assume Rsun = 1, and Lsun = 1
    float integrand = density / ( 1. + r * r );

    // ---- Formulas from Van de Hulst, see also Stokes parameters I Q U V
    btot = (i1table[idx] + i2table[idx]) * integrand;       // I
    bpol = (i2table[idx] - i1table[idx]) * integrand;       // Q

    
    return 0;
}



void PhysicsMie::getConstFactors(float &btf,float &bpf,float &nef, float rho)
{
    float ds = pparentscene->los.ds;
     btf = PI * 0.5 * ds / (k * k);
     bpf = btf;
     nef = RADEG;
}
