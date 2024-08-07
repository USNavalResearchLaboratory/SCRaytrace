
#include <fstream>

#include "config.h"
#include "physicsvariablevsf.h"
#include "rtmiscfunc.h"
#include "Cvec.h"
#include "scene.h"


// const std::string PhysicsVariableVSF::filename = SCRAYTRACE_DATA_DIR "/VSFvariable.dat"; //NOTE: Make Sure to have this .dat file written before trying to create this object
const std::string PhysicsVariableVSF::filename = SCRAYTRACE_DATA_DIR "/VSFvariable_New_Lamy.dat"; //NOTE: Make Sure to have this .dat file written before trying to create this object


PhysicsVariableVSF::PhysicsVariableVSF()
{
  physicsName="VSF for Variable VSF Physics Model.";  
  printvar(physicsName);
  
  printvar(RSUN)
  
  // ---- open file containing the chosen VSF
   ifstream file;
   file.open( filename.c_str(), ios_base::in | ios_base::binary);

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

//   ofstream debug_log; //(DEBUG)
//   debug_log.open("/Users/smidt/Documents/PythonSol/Avery/debug_log.txt", ios_base::out); //(DEBUG)
//   streambuf* oldErrStrm = cerr.rdbuf(debug_log.rdbuf()); //(DEBUG)
//   std::cerr << nbsamp << std::endl; //(DEBUG)
//   std::cerr << "\n" << std::endl; //(DEBUG)
//   for(int i = 0; i < nbsamp; i++){ //(DEBUG)
//     std::cerr << "ang[" << i << "]: " << ang[i] << std::endl; //(DEBUG)
//   } //(DEBUG)
//   std::cerr << "\n" << std::endl; //(DEBUG)
//   for(int i = 0; i < nbsamp; i++){ //(DEBUG)
//     std::cerr << "vsf[" << i << "]: " << vsf[i] << std::endl; //(DEBUG)
//   } //(DEBUG)
//   cerr.rdbuf(oldErrStrm); //(DEBUG)
//   debug_log.close(); //(DEBUG)
  
}



PhysicsVariableVSF::~PhysicsVariableVSF()
{
    delete[] ang;
    delete[] vsf;
}


// void PhysicsVSFVaryDist::initDensityModel(const Cvec &vlosabs)
// {
//     // Transform obsPos from abs to density coordinate system
//     Cvec vlos_inDens;
//     vlos_inDens = ChangetoDensityCoord(pparentscene->modelposition, vlosabs);
//     pparentscene->pmod->initDensityConstFactors(vlos_inDens);
// }



bool PhysicsVariableVSF::computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &density)
{
    // -- compute density at point vs, transforming vs to a frame centered at the Sun; the vector from the Sun to the point
    density = pparentscene->pmod->Density(ChangetoDensityCoord(pparentscene->modelposition, vs));
    
    // ---- compute elongation
    // -- compute distance observer - point of LOS
    Cvec vl = vs - pparentscene->obs.o;
    float l = vl.mag();

    float cosAlpha = pscal(vs, vl) / (r * l);
    if (cosAlpha > 1.) cosAlpha = 1.;
    if (cosAlpha < -1.) cosAlpha = -1.;
    // float theta = acos(-cosAlpha);
    float theta = PI - acos(cosAlpha);

    // //For Debugging forward/backscattering START
    // // btout = density; bpout = 1;

    // // return 0;

    // if(theta < PI/2 && theta >= 0.0) { btout = density * 1; bpout = 1; }
    // else if(density != 0.0){ btout = density * 0.5; bpout = 0.5; }
    // else { btout = 0; bpout = 0;}

    // return 0;
    // //END DEBUG
    
    
    // -- get VSF in lookup table
    unsigned int idx;
    idx = (unsigned int) (float(nbsamp-1) * theta / PI);
    if (idx > (nbsamp-1)) idx = nbsamp-1;
    
    btout = density * pow(pparentscene->obs.o.mag() / ONEAU_RSUN, -0.3) * vsf[idx] / ( r * r );
    bpout = vsf[idx];
    
    // std::cerr << "Btout: " << btout << "; Density: " << density << std::endl; //(DEBUG)
//     cout << "idx : " << idx << ", theta [deg] : " << theta * 180 / PI << ", vsf : " << vsf[idx] << endl;
     
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
void PhysicsVariableVSF::getConstFactors(float &btf,float &bpf,float &nef, float rho)
{
   
    float ds = pparentscene->los.ds;
    btf = ds * RSUN_CM * RSUN * RSUN * 1361;
    bpf = ds;
    nef = ds;
}


