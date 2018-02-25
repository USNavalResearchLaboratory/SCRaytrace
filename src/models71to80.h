//
// $Id: models71to80.h,v 1.1 2010-09-01 19:49:05 thernis Exp $
//
#ifndef MODELS71TO80_H
#define MODELS71TO80_H

#include "CModelBase.h"
#include <string>


//! Leblanc, Dulk, Bougeret electron density model: SolPhy 183: 165-180, 1998
class CModel71 : public CModelBase
{
public:
    float Density(const Cvec &v);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};


//! Density Mie scattering testing: density varies in C / r
class CModel72 : public CModelBase
{
public:
    float Density(const Cvec &v);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float C;        //!> density constant factor
    float dustfreelimit;    //!> dust free zone limit in Rsun

};


//! F-corona and zodiacal light dust density model (Leinart 1976 fan model)
class CModel73 : public CModelBase
{
public:
    float Density(const Cvec &v);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float C;                 //!> density constant factor
    float dustFreeLimit;     //!> dust free zone limit in Rsun
    float decreaseFactor;    //!> density decrease factor in the dust free zone
};


//! F-corona and zodiacal light dust density model (Leinart 1976 fan model), with onion shape decrease bellow 15 Rsun
class CModel74 : public CModelBase
{
public:
    float Density(const Cvec &v);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float C;                 //!> density constant factor
    static const unsigned int SIZEATTENTABLE=11;      //!> onion layer table size
    static const float attenRsun[SIZEATTENTABLE];     //!> onion layer height
    static const float attenPcnt[SIZEATTENTABLE];     //!> attenuation 

};


//! Cylinder model
class CModel75 : public CModelBase
{
public:
    float Density(const Cvec &v);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
protected:
    float r_out;        //!> Outer radius
    float r_in;         //!> Inner radius
    float SWidth;       //!> Semi-width
    float height;       //!> Height of the cylinder center
    float dens;         //!> density
};

//! Torus model
class CModel76 : public CModelBase
{
public:
    float Density(const Cvec &v);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
protected:
//    float l;            //!> Torus axial lenght
    float alpha;        //!> Semi angular width
    float R;            //!> Major radius
    float r_out;        //!> Minor outer radius
    float r_in;         //!> Minor inner radius
    float height;       //!> Height of the axis center
    float dens;         //!> density

//     float alphasur2;
};


//! F-corona and zodiacal light dust density model (Leinart 1976 fan model), with dust enhancement from Kobayashi (Icarus 201 (2009) pp395-405), and others
class CModel77 : public CModelBase
{
public:
  
    ~CModel77();

    float Density(const Cvec &v);
    void initParam(float* pparam);
    void checkData();
    float getEnhanceFactor(float r);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp, int& flagcase);
protected:
    float C;                 //!> density constant factor
    unsigned int modelid;             //!> 0 for silicate, 1 for carbon, 2 Rowan Robinson & May, 3 for Nesvorny, 4 for Liou

    
    static const std::string filename[5];           //!> filenames for the model enhancement profiles
    
    float *rsun, *enhanceFactor;
    unsigned int nbsamp;
    float slope, oord;
};


//! Football model
class CModel78 : public CModelBase
{
public:
    float Density(const Cvec &v);
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);

protected:
//     float alpha;        //!> Semi angular width
//     float R;            //!> Major radius
//     float r_out;        //!> Minor outer radius
//     float r_in;         //!> Minor inner radius
//     float height;       //!> Height of the axis center
//     float dens;         //!> density


};




#endif
