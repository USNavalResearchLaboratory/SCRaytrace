/*! \file model26test.h 
 * \brief Test suite for model 26.
 *
 */

#ifndef MODEL26TEST_H
#define MODEL26TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "models21to30.h"

#define SX 3
#define SY 3
#define SZ 3
#define CNTRX 1.5
#define CNTRY 1.5
#define CNTRZ 1.5
#define VOXSIZERSUN 1.


//! Unit tests of CModel26
class model26Test : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (model26Test);
  CPPUNIT_TEST (testmodel26);
  CPPUNIT_TEST_SUITE_END ();

public:
  void setUp (void);
	void tearDown (void);

protected:
	void testmodel26 (void);

private:

  CModel26 *a;

  float *pparam;

  static const unsigned int sx;
  static const unsigned int sy;
  static const unsigned int sz;
  static const float cntrx;
  static const float cntry;
  static const float cntrz;
  static const float voxsizersun;

	static const float pcube[SX*SY*SZ];

};

const unsigned int model26Test::sx=SX;
const unsigned int model26Test::sy=SY;
const unsigned int model26Test::sz=SZ;
const float model26Test::cntrx=CNTRX;
const float model26Test::cntry=CNTRY;
const float model26Test::cntrz=CNTRZ;
const float model26Test::voxsizersun=VOXSIZERSUN;


const float model26Test::pcube[SX*SY*SZ]={   0,0,0 , 0,0,0 , 0,0,0,
                                             1,1,1 , 1,3,1 , 1,1,10,
                                             2,2,2 , 2,2,2 , 2,2,20};




#endif
