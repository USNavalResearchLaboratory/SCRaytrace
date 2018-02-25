//
// C++ Implementation: cuvemission
//
// $Id: cuvemission.cpp,v 1.2 2010-09-17 15:19:37 thernis Exp $
//
#include "cuvemission.h"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <rpc/rpc.h>
#include <string>
#include <sstream>
#include <cmath>

#include "config.h"
#include "rtmiscfunc.h"

const unsigned int CUVEmission::NBSAMP;
const unsigned int CUVEmission::NBLINES;

CUVEmission::CUVEmission()
{

	isgood=false;
	char *ppathtoxdr;
	string pathtoxdr;
  // ---- get the path where emission files are located
  //ppathtoxdr=getenv("RT_UVEMISSIONPATH");
  pathtoxdr.assign(SCRAYTRACE_DATA_DIR);
  
//   if (ppathtoxdr!=NULL) pathtoxdr.assign(ppathtoxdr);
// 
// 	if (pathtoxdr.empty()) {
// 		cout << "The RT_UVEMISSIONPATH environment variable is not defined\nI just use the current directory" << endl;
// 		pathtoxdr="~/work/cpp/scraytrace/data";
// 	}

//pathtoxdr="/home/arnaud/work/cpp";

  // ---- load the emission vs temperature tables
  // -- yisel
  FILE *pFile;
  string filename=pathtoxdr+"/yisel.xdr";
  pFile = fopen (filename.c_str() , "r");
	if (pFile == NULL) {
		cout << "Error opening file "<< filename << endl;
		exit(1);
  }

	XDR xdrs;
  xdrstdio_create(&xdrs, pFile, XDR_DECODE);

	float a;

   for (unsigned int j=0;j<NBLINES;j++) for (unsigned int i=0;i<NBSAMP;i++) {
		if (!xdr_float(&xdrs, &a)) {
			printf("Error parsing yisel.xdr file!\n");
			fclose(pFile);
			exit(1);
  	}
		yisel[i+j*NBSAMP]=log10(a+1e-30);
	}

  xdr_destroy(&xdrs);
	fclose (pFile);


  // -- emis
  filename.clear();
  filename=pathtoxdr+"/emis.xdr";
  pFile = fopen (filename.c_str() , "r");
	if (pFile == NULL) {
		cout << "Error opening file " << filename << endl;
		exit(1);
  }

  xdrstdio_create(&xdrs, pFile, XDR_DECODE);

  for (unsigned int j=0;j<NBLINES;j++) for (unsigned int i=0;i<NBSAMP;i++) {
		if (!xdr_float(&xdrs, &a)) {
			printf("Error parsing emis.xdr file!\n");
			fclose(pFile);
			exit(1);
  	}
		emis[i+j*NBSAMP]=log10(a);
	}

  xdr_destroy(&xdrs);
	fclose (pFile);


  // -- ti
  filename.clear();
  filename=pathtoxdr+"/ti.xdr";
  pFile = fopen (filename.c_str() , "r");
	if (pFile == NULL) {
		cout << "Error opening file "<< filename << endl;
		exit(1);
  }

  xdrstdio_create(&xdrs, pFile, XDR_DECODE);

  for (unsigned int i=0;i<NBSAMP;i++) {
		if (!xdr_float(&xdrs, &a)) {
			printf("Error parsing ti.xdr file!\n");
			fclose(pFile);
			exit(1);
  	}
		ti[i]=a;
	}

  xdr_destroy(&xdrs);
	fclose (pFile);


  // -- teg
  filename.clear();
  filename=pathtoxdr+"/teg.xdr";
  pFile = fopen (filename.c_str() , "r");
	if (pFile == NULL) {
		cout << "Error opening file "<< filename << endl;
		exit(1);
  }

  xdrstdio_create(&xdrs, pFile, XDR_DECODE);

  for (unsigned int i=0;i<NBSAMP;i++) {
		if (!xdr_float(&xdrs, &a)) {
			printf("Error parsing teg.xdr file!\n");
			fclose(pFile);
			exit(1);
  	}
		teg[i]=a;
	}

  xdr_destroy(&xdrs);

	// -- set isgood flag to true since it looks like everything was alright
  isgood=true;

}


//! Return True if the tables have been loaded properly
bool CUVEmission::IsGood()
{
	return isgood;
}


//! Return yisel
float CUVEmission::getyisel(const unsigned int &i,const unsigned int &j)
{
return yisel[i+j*NBSAMP];
}

//! Compute the emissivity function of the line and the temperature
float CUVEmission::calcEmissivity(const unsigned int &kline,const float &te)
{
if (kline > 4) return 0.;
unsigned int klineoffset=kline-1;
float logte=log10(te);
float yif=nearestneighbor1dinterp(logte,ti,NBSAMP,&yisel[klineoffset*NBSAMP]);
float emiss=nearestneighbor1dinterp(logte,teg,NBSAMP,&emis[klineoffset*NBSAMP]);

return(pow(float(10.),yif)*pow(float(10.),emiss));
}


// $Log: cuvemission.cpp,v $
// Revision 1.2  2010-09-17 15:19:37  thernis
// Add declaration of static constants.
//
// Revision 1.1  2008/12/12 19:54:03  thernis
// Implement UV emission raytracing
//
