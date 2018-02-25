//
// C++ Interface: model67test
//
// $Id: model67test.h,v 1.1 2009/02/09 20:46:29 thernis Exp $
//
#ifndef MODEL67TEST_H
#define MODEL67TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "models61to70.h"

#define SR 3
#define SLON 4
#define SLAT 5
#define SC 2

//! Unit tests of CModel67
class model67Test : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (model67Test);
  CPPUNIT_TEST (testmodel67);
  CPPUNIT_TEST_SUITE_END ();

public:
  void setUp (void);
	void tearDown (void);

protected:
	void testmodel67 (void);

private:

  CModel67 *a,*b;
  float *cube;
  float *cubeb;
	float *cucube;

  static const unsigned int sr;
  static const unsigned int slon;
  static const unsigned int slat;

  static const float r[SR];
  static const float lon[SLON];
  static const float lat[SLAT];

  static const unsigned int CubeTotalSize;

  static const unsigned int offsetnele;
  static const unsigned int offsettemp;

	static const unsigned int sc;


};

const unsigned int model67Test::sr=SR;
const unsigned int model67Test::slon=SLON;
const unsigned int model67Test::slat=SLAT;

const float model67Test::r[SR]={1,2,3};
const float model67Test::lon[SLON]={0,PI/2,PI,3*PI/2};
const float model67Test::lat[SLAT]={-PI/2,-PI/4,0.,PI/4,PI/2};

const unsigned int model67Test::CubeTotalSize=3+SR+SLON+SLAT+SR*SLON*SLAT*2;

const unsigned int model67Test::offsetnele=3+SR+SLON+SLAT;
const unsigned int model67Test::offsettemp=3+SR+SLON+SLAT+SR*SLON*SLAT;


const unsigned int sc=SC;



#endif

// $Log: model67test.h,v $
// Revision 1.1  2009/02/09 20:46:29  thernis
// - Clean up the code
// - Change CModel::Density prototype
// - Update documentation
// - Implement multi-threading using boost thread
// - Add new models
//
