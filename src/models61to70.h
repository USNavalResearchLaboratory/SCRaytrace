/** \file models61to70.h
 * \brief Model 61 to 70
 */

#ifndef MODELS61TO70_H
#define MODELS61TO70_H

#include "CModelBase.h"
#include "models11to20.h"
#include "models21to30.h"


//! \brief Full spherical shell, cartesian coordinates
class CModel61 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
    
  protected:
  float dens; // -- constant Ne for the shell
  float radius; // -- radius of the sphere
  float thickness; // -- Shell thickness
  Cvec cntr; // -- center
  float radminusthick;

};


//! \brief Spherical blob, spherical coordinates
class CModel62 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
    
  protected:
    float dens,radius,thickness,lon,lat,alt;
    Cvec c; // center of the blob in cartesian coordinates
    float radminusthick;
};


//! \brief Full spherical shell
class CModel63 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
  float dens; // -- constant Ne for the shell
  float radius; // -- radius of the sphere
  float thickness; // -- Shell thickness
  float radminusthick;
};


//! \brief Hollow cylinder
class CModel64 : public CModelBase
{
  public:
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
  float hlead,k,b,alpha,dens,thickness,mk,h0;
  int flag;
};


//! \brief Rod filled with density enhancements
//!
//! - The rod is along the Z axis. 
//! - Some cylindrical density enhancements are present within the rod.
//! - Use of Saito equatorial within the rod
//! - Use of Saito polar outside the rod
class CModel65 : public CModelBase
{
  public:
    ~CModel65();
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
 CModel20 *ModelSaitoPolar; //!< Background Ne
 CModel21 *ModelSaitoEquat; //!< Rod Ne
 float RodRadius,strandradius; 
 unsigned int nbstrands;
 float* pstrandspos;
};


//! \brief Wavy neutral sheet with archimedian spiral: Jokipii model, ApJ 1981.
class CModel66 : public CModelBase
{
  public:
    ~CModel66();
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
float alpha; //!< Wave amplitude
float sinalpha; //!< sin of the wave amplitude
float phi0; //!< phase shift
float Omegasun; //!< [rad.s^-1] angular rotation velocity of the Sun
float Vwind; //!< [m.s^-1] wind velocity
float OmegasunoverVwind; //!< Omegasun / Vwind
float SheetThickness; //!< Thickness of the simulated sheet [Rsun]
CModel21 *ModelSaitoEquat; //!< Radial Ne model

};


//! \brief CCMC irregular spherical density cube.
class CModel67 : public CModelBase
{
  public:
   	float Density(const Cvec &v,float &temperature); 
		float Density(const Cvec &v);
    void initParam(float* pparam);
    //void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);

  private:
	unsigned int sr;		//!< size of r
	unsigned int slon;	//!< size of lon
	unsigned int slat;	//!< size of lat
	float *r;						//!< points to first element of r
	float *lon;					//!< first elem of lon
	float *lat;					//!< first elem of lat
	float *nele;				//!< first elem of the electron density cube
	float *temp;				//!< first elem of the temperature cube

};


//! \brief Archimedean spiral: r = a + b * theta
class CModel68 : public CModelBase
{
  public:
    ~CModel68();
    float Density(const Cvec &v); 
    void initParam(float* pparam);
    void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
float a;
float Vwind;
float SheetThickness; //!< Thickness of the simulated sheet [Rsun]
CModel21 *ModelSaitoEquat; //!< Radial Ne model



//float alpha; //!< Wave amplitude
//float sinalpha; //!< sin of the wave amplitude
//float phi0; //!< phase shift
//float Omegasun; //!< [rad.s^-1] angular rotation velocity of the Sun
//float Vwind; //!< [m.s^-1] wind velocity
//float OmegasunoverVwind; //!< Omegasun / Vwind
//float SheetThickness; //!< Thickness of the simulated sheet [Rsun]
//CModel21 *ModelSaitoEquat; //!< Radial Ne model

};


//! \brief GCS or hollow croissant model, based on CModel54, but with sinusoidal variation of Ne along loop axis direction.
class CModel69 : public CModelBase
{
    public:
      float Density(const Cvec &v); 
        void initParam(float* pparam);
        void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
    
    protected:
   float rb,alpha,rf,ratio,k2,mk2,nemin,thick,T_modulation,A_modulation,rs,rs2,rc,rc1,rc2,skinsigmain,skinsigmafr;
};




//! \brief GCS or hollow croissant model, based on CModel54, but with sinusoidal variation of Ne along loop axis direction, centered on circular part loop axis center.
class CModel70 : public CModelBase
{
    public:
      float Density(const Cvec &v); 
        void initParam(float* pparam);
        void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
    
    protected:
   float rb,alpha,rf,ratio,k2,mk2,nemin,thick,T_modulation,A_modulation,rs,rs2,rc,rc1,rc2,skinsigmain,skinsigmafr;
};



#endif
