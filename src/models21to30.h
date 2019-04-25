#ifndef MODELS21TO30_H
#define MODELS21TO30_H

#include "CModelBase.h"
#include "rtmiscfunc.h"


//! Sphericaly symmetric Ne, Saito equatorial density model
class CModel21 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Same as density 14 but with a modified value of the angular sector
//!
//! Streamer slab model
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
class CModel22 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Same as density 18 but with a modified value of the angular sector
//!
//! Streamer slab model
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
//!
//! Ne parameters passing for pminimizerXX.pro (XX=02)
class CModel23 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Cylinder model
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-abs(d/d0))
class CModel24 : public CModelBase
{
	public:
	  	float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Datacube passed in pparam pointer
//!
//! the center of a voxel has [.5,.5,.5] coordinates
//!
//!    element : description
//! -   0       : x size (sx)
//! -   1       : y size (sy)
//! -   2       : z size (sz)
//! -   3       : xc Sun center in pix
//! -   4       : yc Sun center in pix
//! -   5       : zc Sun center in pix
//!               -  note: (0,0,0) is the corner of the first voxel
//! -   6       : voxel size in rsun, same for the 3 directions of space
//! -   7       : data cube in lexicographical order (x,y,z)
class CModel25 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
float voxsize;
long sx,sy;

};
//! Datacube passed in pparam pointer. Same as density 25 but with trilinear interpolation
//!
//! Same as density 25 but with trilinear interpolation
//!
//!    element : description
//!  -  0       : x size (sx)
//!  -  1       : y size (sy)
//!  -  2       : z size (sz)
//!  -  3       : xc Sun center in pix
//!  -  4       : yc Sun center in pix
//!  -  5       : zc Sun center in pix
//!           -   note: (0,0,0) is the corner of the first voxel
//!  -  6       : voxel size in rsun, same for the 3 directions of space
//!  -  7       : data cube in lexicographical order (x,y,z)
//!
//! Still has some issues with edge effect...
class CModel26 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
  protected:
float voxsize;
long sx,sy,sxsy;

};
//! Spherical Shell Ne density by Dr. Russell Howard
class CModel27 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Cylindrical Shell Ne density By Dr. Russell Howward
class CModel28 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Bow Shock Model by Angelos Vourlidas
class CModel29 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Streamer belt simulation with source surface field map
//!
//! User resolution:
//! - pparam[0] : longitude size in pix
//! - pparam[1] : latitude size in pix
//! - pparam[2] : neutral sheet map (lon,lat)
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!    - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
class CModel30 : public CModelBase
{
	public:
		float Density(const Cvec &v); 
		void initParam(float *pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};

#endif

