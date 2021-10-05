//
// File: CModel48.cc


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

