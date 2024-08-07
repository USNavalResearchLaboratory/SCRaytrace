/** \file models81to90.h
 * \brief Model 71 to 80
 */

#ifndef MODELS81TO90_H
#define MODELS81TO90_H

#include "CModelBase.h"
#include <string>


//! \brief F-corona and zodiacal light dust density model with variable VSF for density function modeled as a super ellipse
class CModel81 : public CModelBase
{
public:
    float Density(const Cvec &v);
//     void initDensityConstFactors(const Cvec &vlos_inDens);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float C;                 //!> density constant factor
    float dustFreeStart;     //!> Start of the dust free zone, in Rsun
    float dustFreeEnd;       //!> End of the dust free zone (inner radial dist, where N_dust=0), in Rsun
    float n;                 //!> Shape rectangleness of the superellipse
    float a;                 //!> Semi-major axis of the superellipse
    float b;                 //!> Semi-minor axis of the superellipse
};


//! \brief 3D Rectanglular Prism Model for debugging
class CModel82 : public CModelBase
{
public:
    float Density(const Cvec &v);
//     void initDensityConstFactors(const Cvec &vlos_inDens);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float a;                 //!> X-extent
    float b;                 //!> Y-extent
    float c;                 //!> Z-extent
};


//! \brief F-corona and zodiacal light dust density model with variable VSF for density function modeled as a super ellipse with the semi-major axis being a logisitc function of dist. to . This model was fit on LASCO C3 MMB data.
class CModel83 : public CModelBase
{
public:
    float Density(const Cvec &v);
//     void initDensityConstFactors(const Cvec &vlos_inDens);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float C;                 //!> density constant factor
    float dustFreeStart;     //!> Start of the dust free zone, in Rsun
    float dustFreeEnd;       //!> End of the dust free zone (inner radial dist, where N_dust=0), in Rsun
    float midpoint;          //!> The midpoint constant for the logistic function of the semi-major axis
};


//! \brief F-corona and zodiacal light dust density model with variable VSF for density function modeled as a super ellipse with the semi-major axis being a logisitc function of dist. to Sun. This model was fit on the Lamy et al. 2022
// composite F-Corona map.
class CModel84 : public CModelBase
{
public:
    float Density(const Cvec &v);
//     void initDensityConstFactors(const Cvec &vlos_inDens);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float C;                 //!> density constant factor
    float dustFreeStart;     //!> Start of the dust free zone, in Rsun
    float dustFreeEnd;       //!> End of the dust free zone (inner radial dist, where N_dust=0), in Rsun
    float midpoint;          //!> The midpoint constant for the logistic function of the semi-major axis
};


#endif
