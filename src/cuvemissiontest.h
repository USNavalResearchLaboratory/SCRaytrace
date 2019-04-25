/*! \file cuvemissiontest.cpp
 * \brief Test of the UV Emission physics
 *
 *  
 */


#ifndef CUVEMISSIONTEST_H
#define CUVEMISSIONTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cmath>
#include "cuvemission.h"

//! Unit tests of CUVEmission
class CUVEmissionTest : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (CUVEmissionTest);
  CPPUNIT_TEST (testCUVEmission);
  CPPUNIT_TEST_SUITE_END ();

public:
	void setUp (void);
	void tearDown (void);

protected:
	void testCUVEmission (void);

private:
	CUVEmission *a;

};

#endif

