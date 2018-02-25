//
// File: CModel48.h
// $Id: CModel48.h,v 1.3 2009/02/09 20:51:00 thernis Exp $

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



/*
* $Log: CModel48.h,v $
* Revision 1.3  2009/02/09 20:51:00  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.2  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
