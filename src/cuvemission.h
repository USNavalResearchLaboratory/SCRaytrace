//
// C++ Interface: cuvemission
//
// $Id: cuvemission.h,v 1.1 2008-12-12 19:54:04 thernis Exp $
//
#ifndef CUVEMISSION_H
#define CUVEMISSION_H



using namespace std;


//! UV emission line
class CUVEmission{
public:
	CUVEmission();
	bool IsGood();
  float getyisel(const unsigned int &i,const unsigned int &j);
  float calcEmissivity(const unsigned int &kline,const float &te);

private:
  static const unsigned int NBSAMP=51;
  static const unsigned int NBLINES=4;
	float yisel[NBSAMP*NBLINES];
	float emis[NBSAMP*NBLINES];
	float ti[NBSAMP];
  float teg[NBSAMP];

	bool isgood;

};

#endif

// $Log: cuvemission.h,v $
// Revision 1.1  2008-12-12 19:54:04  thernis
// Implement UV emission raytracing
//
