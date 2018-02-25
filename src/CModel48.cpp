//
// File: CModel48.cc
// $Id: CModel48.cpp,v 1.3 2009/02/09 20:50:59 thernis Exp $

#include "CModel48.h"
#include <iostream>
#include <vector>
#include <cmath>




// Testing model for the adaptative integration
float CModel48::Density(const Cvec &v,float* pmparam,float& temperature)
{	
	
	if (v.v[0] < 0) return 0.;
	if (v.v[1] < 0)	return 0.;	
	
	float a=v.v[2]*v.v[2]/1.;
	if (a > 8E1) return 0.; else return 1000.*exp(-a);
	
}

void CModel48::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=6;
	vp.push_back(moddefparam("","Testing model for the adaptative integration.","",""));	
	
	return;
}



/*
* $Log: CModel48.cpp,v $
* Revision 1.3  2009/02/09 20:50:59  thernis
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
