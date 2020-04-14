
#include "CModel46.h"
#include <vector>
#include <algorithm>
#include "Cvec.h"
#include "Cnurbs.h"

void CModel46::initParam(float* pparam)
{
	// define a nurbs skeleton
	
	// -- define control points	and weight
	Ncontrol=7;
	cp.resize(Ncontrol);
	cp[0]=CControlPoint(4,0,0,1);
	cp[1]=CControlPoint(8,0,1,0.5);
	cp[2]=CControlPoint(6,4,2,1);
	cp[3]=CControlPoint(4,8,3,0.5);
	cp[4]=CControlPoint(2,4,2,1);
	cp[5]=CControlPoint(0,0,1,0.5);
	cp[6]=CControlPoint(4,0,0,1);

	// -- define knot vector
	//KnotVec.resize(Ncontrol+Korder);
	Korder=3;
	KnotVec.push_back(0);
	KnotVec.push_back(0);
	KnotVec.push_back(0);
	KnotVec.push_back(1);
	KnotVec.push_back(1);
	KnotVec.push_back(2);
	KnotVec.push_back(2);
	KnotVec.push_back(3);
	KnotVec.push_back(3);
	KnotVec.push_back(3);

	Nbvertices=50;

	// -- calculate nurbs curve
	Cnurbs *pnurbs=new Cnurbs(Ncontrol,KnotVec,Nbvertices,cp);
	
	// -- get the curve
	curve=pnurbs->getCurve();

	//pnurbs->printCurve();
	
	//std::cout << "Size Curve : " << curve.size() << std::endl;


	// -- don't need the nurbs object anymore: free memory
	delete pnurbs;
	
	// -- init dist vector
	dist.resize(Nbvertices);
	
}

// Model using nurbs skeleton.
float CModel46::Density(const Cvec &v,float* pparam,float& temperature)
{
	// -- compute distance between curve and requested point
	for (int i=0;i<Nbvertices;i++)
		dist[i]=(curve[i]-v).norm();
	
	// -- find the minimum distance
	itermin=std::min_element(dist.begin(),dist.end());

	//float mindist=*itermin;
	
	//unsigned int posmin=distance(dist.begin(),itermin);
	
	// -- compute density depending on distance to skeleton
	float neout=0;
	
	if (*itermin < .5) neout=1e4;

	return neout;
}

void CModel46::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=3;
	vp.push_back(moddefparam("","Model using nurbs skeleton.","",""));	
	
	return;
}

