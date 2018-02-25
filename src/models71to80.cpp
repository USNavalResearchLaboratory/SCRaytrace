//
// $Id: models71to80.cpp,v 1.1 2010-09-01 19:49:05 thernis Exp $
//

#include "models71to80.h"
#include <cmath>
#include <string>
#include <fstream>
#include "config.h"

float CModel71::Density(const Cvec &v)
{
    float r=v.norm();
    if (r <= 1.05) return 0.;

    // -- see article SolPhy 183: 165-180, 1998
    /*  float a=3.3e5;
      float b=4.1e6;
      float c=8.0e7;*/

    float nel=3.3e5*pow(r,-2)+4.1e6*pow(r,-4)+8.0e7*pow(r,-6);

    return nel;
}
void CModel71::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase=0x4;
    vp.push_back(moddefparam("","Leblanc, Dulk, Bougeret model","",""));
    return;
}


// Density 72: Mie scattering testing. Density varying in 1 / r.
float CModel72::Density(const Cvec &v)
{
    if (v.mag() < dustfreelimit) return 0.;
    float dens;
    dens = C / v.mag();
    return dens;
}

// Inititialization of the parameters
void CModel72::initParam(float* pparam)
{
    C = pparam[0];               // -- constant factor
    dustfreelimit = pparam[1];   // -- dust free zone limit in Rsun
}

void CModel72::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Mie scattering testing", "", ""));
    vp.push_back(moddefparam("C", "1.", "Density", "particules/cm^3"));
    vp.push_back(moddefparam("dustfreelimit", "10.", "Dust free zone limit", "Rsun"));
    return;
}




// Density 73: Fan model
float CModel73::Density(const Cvec &v)
{
  // ---- compute BetaPi, the elevation above the ecliptic plane
  //      which is assumed to be the (O,y,z) plane.
  float betapi = atan2(v[0],sqrt(v[1] * v[1] + v[2] * v[2]));
  float absSinBetaPi = fabs(sin(betapi));
  // -- Lamy Perrin 1986
  //   float dens = C / v.mag() * exp(-3.5 * pow(absSinBetaPi, 0.775 + 0.624 * absSinBetaPi));
  
  // -- Leinart 1976
  float dens = C / v.mag() * exp(-2.6 * pow(absSinBetaPi,1.3));
  
  // -- Leinart 1978
  //   float dens = C * pow(v.mag(), -1.3) * exp(-2.1 * absSinBetaPi);

  if (v.mag() <= dustFreeLimit) dens *= decreaseFactor;
  
  return dens;
}

// Inititialization of the parameters
void CModel73::initParam(float* pparam)
{
    C = pparam[0];               // -- constant factor
    dustFreeLimit = pparam[1];   // -- dust free zone limit in Rsun
    decreaseFactor = pparam[2];  // -- factor decrease 
}

void CModel73::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Mie scattering testing", "", ""));
    vp.push_back(moddefparam("C", "1.", "Density", "particules/cm^3"));
    vp.push_back(moddefparam("dustFreeLimit", "5.", "Dust free zone limit", "Rsun"));
    vp.push_back(moddefparam("decreaseFactor", "0.9", "Density decrease in the dust free zone", ""));
    return;
}





// Density 74: Leinart fan model (1976), with onion decrease bellow 15 Rsun
float CModel74::Density(const Cvec &v)
{
  // ---- compute BetaPi, the elevation above the ecliptic plane
  //      which is assumed to be the (O,y,z) plane.
  float betapi = atan2(v[0],sqrt(v[1] * v[1] + v[2] * v[2]));
  float absSinBetaPi = fabs(sin(betapi));
  
  // -- Leinart 1976
  float r = v.mag();
  float dens = C / r * exp(-2.6 * pow(absSinBetaPi,1.3));
  
  // -- attenuate depending on the depth of the onion layer
  if (r > attenRsun[0]) return dens;
  unsigned int i=0;
  while (r <= attenRsun[i]) i++;
  return dens * attenPcnt[i-1];
}

// Inititialization of the parameters
void CModel74::initParam(float* pparam)
{
    C = pparam[0];               // -- constant factor
}


const float CModel74::attenRsun[] = {15., 13, 12, 11, 10,  9,  8,  7,  6, 5 , 0.};
const float CModel74::attenPcnt[] = {0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1, 0., 0.};


void CModel74::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Mie scattering with onion dust free zone", "", ""));
    vp.push_back(moddefparam("C", "1.", "Density", "particules/cm^3"));
    return;
}






// Density 75: Cylinder model
float CModel75::Density(const Cvec &v)
{ 
  // The cylinder axis is X
  
  // -- return 0 if outside the width
  if (fabs(v[0]) > SWidth) return 0.;
  
  // -- compute distance point to cylinder axis
  float r = (Cvec(0., v[1], v[2] - height)).mag();
//   float r = (Cvec(0., v[1], v[2])).mag();
  
  if (r > r_out || r < r_in) return 0.;
  
  return dens;
}

// Inititialization of the parameters
void CModel75::initParam(float* pparam)
{
    r_out  = pparam[0];     // Outer radius
    r_in   = pparam[1];     // inner radius
    SWidth = pparam[2];     // Semi-width
    height = pparam[3];     // height of the cylinder center
    dens   = pparam[4];     // density
}


void CModel75::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Cylinder", "", ""));
    vp.push_back(moddefparam("r_out", "1.", "Outer radius", "Rsun"));
    vp.push_back(moddefparam("r_in", "0.8", "Inner radius", "Rsun"));
    vp.push_back(moddefparam("SWidth", "0.5", "Semi-width", "Rsun"));
    vp.push_back(moddefparam("height", "5.", "Height of the center", "Rsun"));
    vp.push_back(moddefparam("dens", "1.", "Density", "electron/cm^3"));
    return;
}



// Density 76: Torus model
float CModel76::Density(const Cvec &v)
{ 
  // The torus is in the X, Z plane.
  // Height is in direction of Z axis
  
//   Cvec OC = Cvec(0, 0, h-R);
//   Cvec OH = Cvec(v[0], 0, v[2]);
  
//   if (v[0] < 0 ) return 100.;
  
  
//   Cvec CH = Cvec(v[0], 0, v[2] - height + R);
  
//   float gamma = acos( pscal(Cvec(0,0,1), CH) / CH.norm() );
//   float gamma = acos( CH[2] / CH.norm() );
  float gamma = atan2(v[0],v[2] - height + R);
//   return gamma;
  
  
  // -- return 0. if outside the angular range
//   if (fabs(gamma) > alphasur2) return 0.;
  if (fabs(gamma) > alpha) return 0.;

  Cvec MP = Cvec(   R * sin(gamma), 
                    0, 
                    height - R + R * cos(gamma)) - v;
  
  float normMP = MP.norm();
  
  if (normMP > r_out || normMP < r_in) return 0.; //else return dens;
  
  return dens;
}

// Inititialization of the parameters
void CModel76::initParam(float* pparam)
{
//     l      = pparam[0];     // axis lenght
    alpha  = pparam[0];     // semi-angular width
    R      = pparam[1];     // radius
    r_out  = pparam[2];     // minor outer radius
    r_in   = pparam[3];     // minor inner radius
    height = pparam[4];     // height of the axis center
    dens   = pparam[5];     // density
    
//     alphasur2 = l / R / 2.;
    
//     std::cout << "alphasur2 : " << alphasur2 << std::endl;
}


void CModel76::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Torus", "", ""));
//     vp.push_back(moddefparam("l", "0.5", "Axis lenght", "Rsun"));
    vp.push_back(moddefparam("alpha", "0.5", "Semi angular width", "rad"));
    vp.push_back(moddefparam("R", "2.", "Radius", "Rsun"));
    vp.push_back(moddefparam("r_out", "0.5", "Outer minor radius", "Rsun"));
    vp.push_back(moddefparam("r_in", "0.3", "Inner minor radius", "Rsun"));
    vp.push_back(moddefparam("height", "2.", "Semi-width", "Rsun"));
    vp.push_back(moddefparam("dens", "10000.", "Density", "electron/cm^3"));
    return;
}





// Density 77: F corona with Kobayashi dust enhancement
float CModel77::Density(const Cvec &v)
{
  // ---- compute BetaPi, the elevation above the ecliptic plane
  //      which is assumed to be the (O,y,z) plane.
  float betapi = atan2(v[0],sqrt(v[1] * v[1] + v[2] * v[2]));
  float absSinBetaPi = fabs(sin(betapi));
  
  // -- Leinart 1976
  float r = v.mag();
  float dens = C / r * exp(-2.6 * pow(absSinBetaPi,1.3));
  
  // -- multiply by enhancement factor 
  float ef;
  ef = this->getEnhanceFactor(r);
  
  return dens * ef;
}


// Inititialization of the parameters
void CModel77::initParam(float* pparam)
{
    C = pparam[0];               // -- constant factor
    
    // modelid:
    //  0 for silicate
    //  1 for carbon
    //  2 for Rowan Robinson & May
    //  3 for Nesvorny
    //  4 for Liou
    modelid = (unsigned int) pparam[1];
    
    if (modelid > 4) modelid = 4;
    
    std::cout << "init param of CModel77: enhancement model id " << modelid << std::endl;
      
    std::cout << "reading file : " << filename[modelid] << std::endl;

      // ---- open file the requested density enhancement
   std::ifstream file;
   file.open(filename[modelid].c_str(), std::ios_base::in | std::ios_base::binary);
  
  // -- First record is the size of the data
  int *pnbsamp;
  pnbsamp = new int;
  file.read((char*)pnbsamp, sizeof(int));
//   printvar(*pnbsamp);
  
  nbsamp = *pnbsamp;
  
  delete pnbsamp;
  
  // -- create lookup table
  rsun = new float[nbsamp];
  enhanceFactor = new float[nbsamp];
  
  // -- now read the lookup table data
  file.read((char*)rsun, sizeof(float) * nbsamp);
  file.read((char*)enhanceFactor, sizeof(float) * nbsamp);
  
  file.close();

//   printvar(ang[0]);
//   printvar(vsf[0]);

  slope = (rsun[nbsamp-1] - rsun[0]) / float(nbsamp);
  oord = rsun[0];
    
}



void CModel77::checkData()
{
  std::cout << "rsun[1] : " << rsun[1] << std::endl;
  std::cout << "enhanceFactor[1] : " << enhanceFactor[1] << std::endl;
}


float CModel77::getEnhanceFactor(float r)
{
  int i;
  i = int((r - oord) / slope);
  
  if (i < 0) return 0.; 
  if (i >= nbsamp) return 1.;
  
  return enhanceFactor[i];
}


// -- destructor
CModel77::~CModel77()
{
  delete[] rsun, enhanceFactor;
}


const std::string CModel77::filename[] = {
    SCRAYTRACE_DATA_DIR  "/binsilicate_n100_s1_5_trimmed.dat",
    SCRAYTRACE_DATA_DIR  "/bincarbon_n100_s1_5_trimmed.dat",
    SCRAYTRACE_DATA_DIR  "/bindustmix_Rowan-Robinson12_sil_62.336_carb_37.664.dat",
    SCRAYTRACE_DATA_DIR  "/bindustmix_Nesvorny10_sil_52.5_carb_47.5.dat",
    SCRAYTRACE_DATA_DIR  "/bindustmix_Liou95_sil_66.65_carb_33.35.dat"};


void CModel77::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "F corona with Kobayashi enhancement", "", ""));
    vp.push_back(moddefparam("C", "1.", "Density", "particules/cm^3"));
    vp.push_back(moddefparam("modelid", "0", "0: silicate, 1: carbon", ""));
    return;
}



// Density 78: football model
float CModel78::Density(const Cvec &v)
{ 

    float dens;
    dens = 0.;
  
    return dens;
}

// Inititialization of the parameters
void CModel78::initParam(float* pparam)
{
//    alpha  = pparam[0];     // 

}


void CModel78::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Football", "", ""));
    return;
}





