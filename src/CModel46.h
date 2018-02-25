//
// File: CModel46.h
// $Id: CModel46.h,v 1.3 2009/02/09 20:50:59 thernis Exp $

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



/*
* $Log: CModel46.h,v $
* Revision 1.3  2009/02/09 20:50:59  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.2  2007/05/14 17:19:39  thernis
* Add CVS id and log in all files
*
*/
