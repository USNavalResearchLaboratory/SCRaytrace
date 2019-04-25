#ifndef MODELS11TO20_H
#define MODELS11TO20_H

#include <vector>
#include "CModelBase.h"

//! Streamer belt simulation with source surface field map.
//!
//! Fixed resolution: longitude 360 pix, latitude 181 pix.
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!    - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
//!
class CModel11 : public CModelBase
{
	public:
	  float Density(const Cvec &v);
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
	protected:
		float *pnsheetmap;
};
//! Simple neutral sheet: sin(theta)
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-abs(d/d0))
//!
class CModel12 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Sphericaly symmetric Ne
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r))
//!    - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
class CModel13 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Sphericaly symmetric Ne
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
class CModel14 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Sphericaly symmetric Ne
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
class CModel15 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void initParam(float* pparam);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
	protected:
		float *pmodul;
};
//! Cylinder model
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-abs(d/d0))
class CModel16 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Cylinder model with parameters passing
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-abs(d/d0))
class CModel17 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Sphericaly symmetric Ne
//!
//! Use the Hayes model
//! - Ne(r)=Sum_k alpha_k * r^(-k) , k=1,2,3,4
//! - Orthoradial model : exp(-d^2/w0(r)) and exp(-abs(d/u0(r)))
//!   - with w0(r) = Sum_k beta_k * r^(-k)), k=1,2,3,4
//!
//! Ne parameters passing for pminimizerXX.pro (XX=02)
class CModel18 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Sphericaly symmetric Ne, Saito coronal hole density model: Saito, Poland, Munro, SolPhy 55 (1977) pp. 121-134
class CModel19 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};
//! Sphericaly symmetric Ne, Saito polar density model
class CModel20 : public CModelBase
{
	public:
	  float Density(const Cvec &v); 
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);
};

#endif
