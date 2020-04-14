/*! \file CModel46.h
 * \brief Model using a NURBS skeleton.
 *
 */

#ifndef _CMODEL46_H_
#define _CMODEL46_H_

#include "CModelBase.h"
#include "CControlPoint.h"

//! Model using a NURBS skeleton
class CModel46 : public CModelBase
{
	public:
		// member function declarations...
		void initParam(float* pparam);
		float Density(const Cvec &v,float* pparam,float& temperature);
		//int getNbparam() {return 0;};
		void dumpDefaultParamForIDL(std::vector<moddefparam> & vp,int & flagcase);

	protected:
		std::vector<float> KnotVec;
		std::vector<CControlPoint> cp;
		int Ncontrol,Nbvertices,Korder;
		std::vector<Cvec> curve;  //!< the points of the nurbs curve 
		std::vector<float> dist;
		std::vector<float>::iterator itermin;
};


#endif	//_CMODEL46_H_


