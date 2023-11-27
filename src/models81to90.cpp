
#include "models81to90.h"
#include "constant.h"
#include <cmath>
#include <string>
#include <fstream>
#include "config.h"
#include "scene.h"

// Density 81: Variable VSF Model for density function modeled as a superellipse
float CModel81::Density(const Cvec &v)
{    
    // -- Express obs.o in the density reference frame
    Cvec obsPos_inNe = ChangetoDensityCoord(pparentscene->modelposition, pparentscene->obs.o);
    
    // -- Compute vector OV
    Cvec OV = v - obsPos_inNe;

    float y = v[1];
    float yabs = fabs(y);
    float x0abs = fabs(v[0]);
    float alpha = asin(x0abs/v.norm());
    float r = (a*b)/(pow(pow(b, n)*pow(cos(alpha), n) + pow(a, n)*pow(sin(alpha), n), 1/n));
    
    float dens = (1/v.norm()) * r * C;

    
    if (v.mag() <= dustFreeStart) {
        dens *= (v.mag() - dustFreeEnd) / (dustFreeStart - dustFreeEnd);
        if (dens < 0.) dens = 0.;
    }
    return dens;
}


// Inititialization of the parameters
void CModel81::initParam(float* pparam)
{
    C = pparam[0];               // -- constant factor
    dustFreeStart = pparam[1];   // -- Start of the dust free zone, in Rsun
    dustFreeEnd = pparam[2];     // -- End of the dust free zone (inner radial dist, where N_dust=0), in Rsun
    n = pparam[3];               // -- rectangleness of the superellipse
    a = pparam[4];               // -- semi-major axis of the superellipse
    b = pparam[5];               // -- semi-minor axis of the superellipse
}

void CModel81::dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Fan model from Lamy Perrin 1986", "", ""));
    vp.push_back(moddefparam("C", "1.", "Density", "particules/cm^3"));
    vp.push_back(moddefparam("dustFreeLimit", "5.", "Dust free zone limit", "Rsun"));
    vp.push_back(moddefparam("decreaseFactor", "0.9", "Density decrease in the dust free zone", ""));
    return;
}


// Density 82: 3D Rectangular Prism Model For Debugging
float CModel82::Density(const Cvec &v)
{    
    // -- Express obs.o in the density reference frame
    Cvec obsPos_inNe = ChangetoDensityCoord(pparentscene->modelposition, pparentscene->obs.o);
    
    // -- Compute vector OV
    Cvec OV = v - obsPos_inNe;

    if(v[0] > a){
        return 0;
    }

    if(v[1] > b){
        return 0;
    }

    if(v[2] > c){
        return 0;
    }

    return 1;
}


// Inititialization of the parameters
void CModel82::initParam(float* pparam)
{
    a = pparam[0];
    b = pparam[1];
    c = pparam[2];
}

void CModel82::dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Fan model from Lamy Perrin 1986", "", ""));
    vp.push_back(moddefparam("C", "1.", "Density", "particules/cm^3"));
    vp.push_back(moddefparam("dustFreeLimit", "5.", "Dust free zone limit", "Rsun"));
    vp.push_back(moddefparam("decreaseFactor", "0.9", "Density decrease in the dust free zone", ""));
    return;
}


// Density 83: Variable VSF Model for density function modeled as a superellipse and the semi-major axis is logistic with respect to the distance from the Sun
float CModel83::Density(const Cvec &v)
{    
    // -- Express obs.o in the density reference frame
    Cvec obsPos_inNe = ChangetoDensityCoord(pparentscene->modelposition, pparentscene->obs.o);
    
    // -- Compute vector OV
    Cvec OV = v - obsPos_inNe;

    //The commented code below was used for testing various functions for logistic growth rates j and k to see which worked best for LASCO and SECCHI data.

    // float k = v.norm() / 215;
    // float k = obsPos_inNe.norm() / 215;
    // float k = 0.15;
    // float k = obsPos_inNe.norm();
    // float k = 2 / (1 + exp(-0.6*(obsPos_inNe.norm() - 210)));
    float k = 1 / (1 + exp(-0.07*(obsPos_inNe.norm() - 230)));
    float j = k / 2;
    // float k = 0.3;
    // float j = 0.1;
    // std::cout << "k: " << k << std::endl;
    // exit(0);
    // float a = a0 / (1 + exp(-k * (sqrt(v[0]*v[0] + v[1]*v[1]) - midpoint)));
    float a = a0 / (1 + exp(-k * (v.norm() - midpoint)));
    float b = b0 / (1 + exp(-j * (v.norm() - midpoint)));

    float y = v[1];
    float yabs = fabs(y);
    float x0abs = fabs(v[0]);
    float alpha = asin(x0abs/v.norm());
    float r = (a*b)/(pow(pow(b, n)*pow(cos(alpha), n) + pow(a, n)*pow(sin(alpha), n), 1/n));
    
    float dens = (1/v.norm()) * r * C;

    
    if (v.mag() <= dustFreeStart) {
        dens *= (v.mag() - dustFreeEnd) / (dustFreeStart - dustFreeEnd);
        if (dens < 0.) dens = 0.;
    }
    return dens;
}


// Inititialization of the parameters
void CModel83::initParam(float* pparam)
{
    C = pparam[0];               // -- constant factor
    dustFreeStart = pparam[1];   // -- Start of the dust free zone, in Rsun
    dustFreeEnd = pparam[2];     // -- End of the dust free zone (inner radial dist, where N_dust=0), in Rsun
    n = pparam[3];               // -- rectangleness of the superellipse
    a0 = pparam[4];              // -- semi-major axis of the superellipse
    b0 = pparam[5];              // -- semi-minor axis of the superellipse
    midpoint = pparam[6];        // -- midpoint constant of logistic function for the semi-major axis
}

void CModel83::dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase)
{
    flagcase = 0;
    vp.push_back(moddefparam("", "Fan model from Lamy Perrin 1986", "", ""));
    vp.push_back(moddefparam("C", "1.", "Density", "particules/cm^3"));
    vp.push_back(moddefparam("dustFreeLimit", "5.", "Dust free zone limit", "Rsun"));
    vp.push_back(moddefparam("decreaseFactor", "0.9", "Density decrease in the dust free zone", ""));
    return;
}
