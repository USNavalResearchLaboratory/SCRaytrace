//
// File: CModel48.h

#ifndef _CMODEL48_H_
#define _CMODEL48_H_

#include "CModelBase.h"

//! Testing model for the adaptative integration
class CModel48 : public CModelBase
{
	public:
		float Density(const Cvec &v,float* pparam,float& temperature);
		void dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase);

	protected:
	
};


#endif	//_CMODEL48_H_


