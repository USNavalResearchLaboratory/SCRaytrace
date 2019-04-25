
#include <fstream>

#include "config.h"
#include "physicsvsf.h"
#include "rtmiscfunc.h"
#include "Cvec.h"
#include "scene.h"

// const std::string PhysicsVSF::filename = SCRAYTRACE_DATA_DIR "/binsilicate_n100_s1_5_trimmed.dat";
const std::string PhysicsVSF::filename = SCRAYTRACE_DATA_DIR "/VSFLamyPerrinPlaneofSymInterp06.dat";


PhysicsVSF::PhysicsVSF()
{
  physicsName="VSF for zodiacal light simulation";  
  printvar(physicsName);
  
  // ---- open file containing the LamyPerrin VSF
  //      See Lamy and Perrin, A&A 163, pp269-286 (1986).
   ifstream file;
   file.open( filename.c_str(), ios_base::in | ios_base::binary);

//    file.open("/home/thernis/work/mie/VSFLamyPerrinPlaneofSymInterp06.dat",ios_base::in |ios_base::binary);
  
  // -- First record is the size of the data
  int *pnbsamp;
  pnbsamp = new int;
  file.read((char*)pnbsamp, sizeof(int));
  printvar(*pnbsamp);
  
  nbsamp = *pnbsamp;
  
  delete pnbsamp;
  
  // -- create lookup table
  ang = new float[nbsamp];
  vsf = new float[nbsamp];
  
  // -- now read the lookup table data
  file.read((char*)ang, sizeof(float) * nbsamp);
  file.read((char*)vsf, sizeof(float) * nbsamp);
  
  file.close();

  printvar(ang[0]);
  printvar(vsf[0]);
  
  
  
}



PhysicsVSF::~PhysicsVSF()
{
    delete[] ang;
    delete[] vsf;
}




bool PhysicsVSF::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &density)
{
    // -- compute density at point vs
    density = pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition,vs));
    
    // ---- compute elongation
    // -- compute distance observer - point of LOS
    Cvec vl = vs - pparentscene->obs.o;
    float l = vl.mag();

    float cosAlpha = pscal(vs,vl) / (r * l);
    if (cosAlpha > 1.) cosAlpha = 1.;
    if (cosAlpha < -1.) cosAlpha = -1.;
    float theta = acos(-cosAlpha);
    
    
    // -- get VSF in lookup table
    unsigned int idx;
    idx = (unsigned int) (float(nbsamp-1) * theta / PI);
    if (idx > (nbsamp-1)) idx = nbsamp-1;
    
    btout = density * vsf[idx] / ( r * r );
    bpout = vsf[idx];
     
    return 0;
}

/*!
 * The VSF used is in units of cm^-1.sr^-1. The RSUN_CM converts it into Rsun^-1.sr^-1 units.
 * 
 * The solid angle of the Sun viewed from the element of volume is ~PI Rsun^2 / r^2. 
 * Since the distance are computed in Rsun, the Rsun^2 is equal to 1.
 * 
 * The 1361 factor corrsponds to the solar flux, in W/m^2, received at 1 AU.
 * 
 * See also: "Light Scattering by Solar System Dust: Image Reconstruction of the Lunar Sunrise Sketches Drawn by the Apollo 17 Crew", Niklas Siipola, Master's thesis, University of Oulu, Spring 2017
 * 
 * 
 */
void PhysicsVSF::getConstFactors(float &btf,float &bpf,float &nef, float rho)
{
    float ds = pparentscene->los.ds;
    btf = ds * RSUN_CM * RSUN * RSUN * 1361;
    bpf = ds;
    nef = ds;
}


