//
// File: Cnurbs.cc

#include "Cnurbs.h"
#include <iostream>
#include "Cvec.h"


Cnurbs::~Cnurbs()
{
	delete3Darray();
}

Cnurbs::Cnurbs(unsigned int Ncontrol0,std::vector<float> KnotVec0,unsigned int Nbvertices0,std::vector<CControlPoint> cp0)
{
	Ncontrol=Ncontrol0;
	KnotVec=KnotVec0;
	KnotSize=KnotVec.size();
	Korder=KnotSize-Ncontrol;
	Nbvertices=Nbvertices0;
	
	
	// -- normalize Knot vector
	for (unsigned i=0;i<KnotVec.size();i++)
		KnotVec[i]=KnotVec[i]/KnotVec.back();
	
	// -- normalize vertices vector
	for (unsigned int i=0;i<Nbvertices;i++) t.push_back(VMIN+VDIFF*(((float) i)/((float) Nbvertices-1)));
	
//	std::cout << "t : ";
//	for (unsigned i=0;i<t.size();i++)
//		std::cout << t[i] << ",";
//	std::cout << std::endl;
	
	cp=cp0;

	allocate3Darray();
	
	calcBasis();
	calcCurve();

	
}

void Cnurbs::printBasis(unsigned int Ko)
{
	unsigned int j;
	if (Ko > Korder) 
	{
		// I should throw an exception here instead
		std::cout << "nurbs order too big " << std::endl;
		return;
	}
	
	for (unsigned int i=0;i<KnotSize;i++)
	{
		std::cout << "i,k : " << i << "," << Ko << std::endl;
		for (j=0;j< (Nbvertices-1);j++)
			std::cout << Nik[j][i][Ko] << ",";
		std::cout << Nik[j][i][Ko] << "END" << std::endl;
	}
}

void Cnurbs::printCurve()
{
	for (unsigned i=0;i<curve.size();i++)
	{
		std::cout << "point : " << i << ", xyz : " << curve[i] << std::endl;
	}	
}

std::vector<Cvec> Cnurbs::getCurve()
{
	return curve;
}


void Cnurbs::calcBasis()
{
	for (unsigned int i=1;i<KnotSize;i++)
	{
		for (unsigned int j=0;j<Nbvertices;j++)
		{
			if (t[j] >= KnotVec[i-1] && t[j] < KnotVec[i+1-1]) 
				Nik[j][i-1][0]=1;
			else
				Nik[j][i-1][0]=0;
		}
	}

	for (unsigned int kk=2;kk<=Korder;kk++)
		for (unsigned int i=1;i<=(KnotSize-kk);i++)
		{
	        float kipkmki=KnotVec[i+kk-1 -1]-KnotVec[i -1];
    	    float kipkmkip1=KnotVec[i+kk -1]-KnotVec[i+1 -1];
			for (unsigned int j=0;j<Nbvertices;j++)
			{
				Nik[j][i -1][kk -1]=((kipkmki != 0) ? 
					(t[j]-KnotVec[i -1])*Nik[j][i -1][kk-1 -1]/kipkmki : 0.)+
					((kipkmkip1 != 0) ? 
					(KnotVec[i+kk -1]-t[j])*Nik[j][i+1 -1][kk-1 -1]/kipkmkip1 : 0.);
			}
		}
}


int Cnurbs::allocate3Darray()
{
//	Nik float[Nbvertices][KnotSize][Korder];
	Nik = new float**[Nbvertices];
	for (unsigned int i = 0; i < Nbvertices; i++)
	{
        Nik[i] = new float*[KnotSize];
		for (unsigned int j=0;j < KnotSize;j++)
			Nik[i][j]=new float[Korder];
	}


	return 1;
}

int Cnurbs::delete3Darray()
{
	for (unsigned int i = 0; i < Nbvertices ;i++)
	{
		for (unsigned int j=0;j < KnotSize;j++)
			delete[] Nik[i][j];
        delete[] Nik[i];
		}
    delete[] Nik;
	return 1;
}



void Cnurbs::calcCurve()
{
	//std::vector<Cvec> curve;  // the points of the nurbs curve 

	Cvec pup;
	float plo;

	for (unsigned int j=0;j<Nbvertices;j++)
	{
		plo=0;pup=Cvec(0,0,0);
		for (unsigned int i=1;i<=Ncontrol;i++)
		{	
			pup+=cp[i -1].p*cp[i -1].w * Nik[j][i -1][Korder -1];
			plo+=cp[i -1].w*Nik[j][i -1][Korder -1];
		}
		curve.push_back(pup/plo);
	}

}



