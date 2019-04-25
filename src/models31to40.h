#ifndef MODELS31TO40_H
#define MODELS31TO40_H

#include "CModelBase.h"

//! Graduated Curved Cylindrical Shell Model
//!
//!       Flux rope model 8/19/04.
//!
//!		The top is now a hemisphere attached to the cone legs.
//!
//! The size of the circular cross section of the shell increases uniformly
//!   from the base to the top.
//!
//!  Electron density calculation: density is gradient from the base to the top.
class CModel31 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Streamer belt simulation with source surface field map: based on density 30
//!
//! d1 and d2 parameters passed in argument: they are constant (do no depend on r). Useful for circular profile automatic fit (well at least attempt of automatic fit !)
//!
//! User resolution:
//! - pparam[0] : longitude size in pix
//! - pparam[1] : latitude size in pix
//! - pparam[2] : d1 parameter
//! - pparam[3] : d2 parameter
//! - pparam[4] : neutral sheet map (lon,lat)
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
class CModel32 : public CModelBase
{
	public:
	  float Density(const Cvec &v);
		void initParam(float* pparam);
};
//! Tube shell model
class CModel33 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Test of a simple density rope
//!
//! raytracewl,sbt1,imsize=[512,512],losrange=[-15,15],losnbp=100,/c3,modelid=34
class CModel34 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Density and temperature cubes from S.Gibson and B.Low CME progam model.
//!
//! See "A time-dependent three-dimementional magnetohydrodynamic model of the coronal mass ejection", S.E.Gibson, B.C.Low, The Astrophiscal Journal, 493:460-473, 1998 January 20.
//!
//! - Use gibsoncmewrapper.pro to generate the density and temperature cube
//! - Use load_gibsondens.pro to load and build the modparam model input 
//! vector parameter
class CModel35 : public CModelBase
{
	public:
	  float Density(const Cvec &v,float &temperature); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Streamer belt simulation with source surface field map
//!
//! User resolution:
//! - pparam[0] : longitude size in pix
//! - pparam[1] : latitude size in pix
//! - pparam[2] : neutral sheet map (lon,lat)
//!
//! use Guhathakurta model: see CModel01
class CModel36 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Streamer belt simulation with source surface field map: distance from neutral line is given by the magnetic field value
//!
//! distance from neutral line is given by 
//! the magnetic field value
//! - pparam[0] : streamer belt thickness coefficient
//! User resolution:
//! - pparam[1] : longitude size in pix
//! - pparam[2] : latitude size in pix
//! - pparam[3] : neutral sheet map (lon,lat)
//!
//! Use the Hayes model:
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
class CModel37 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Slab model, same as density 23 but with modulation of the slab FO
//!
//! Streamer slab model
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!    - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
//! Ne parameters passing for pminimizerXX.pro (XX=  )
//! - pparam[0..3] : coeff of the Ne polynomial
//! - pparam[4] : slab half angle
//! - pparam[5] : subscript max of the modulation profile (size - 1)
//! - pparam[6] : first element of the modulation profile
class CModel38 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Test of a simple helix
//!
//! raytracewl,sbt1,imsize=[512,512],losrange=[-15,15],losnbp=100,/c3,modelid=34
class CModel39 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Streamer belt simulation with source surface field map
//!
//! From pfss simulation, the last onion layer (generally 2.5)
//! is used. The electron density is scaled on that value. The mag field 
//! is supposed to be radial after
//!
//! User resolution:
//! - pparam[0] : longitude size in pix
//! - pparam[1] : latitude size in pix
//! - pparam[2] : height of the layer in Rsun
//! - pparam[3] : Ne layer map (lon,lat)
//! use the Saito model scaled on the map value
class CModel40 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};

#endif

