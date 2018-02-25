//
// C++ Interface: cuvemissiontest
//
// $Id: cuvemissiontest.h,v 1.1 2008-12-12 19:54:06 thernis Exp $
//
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

// $Log: cuvemissiontest.h,v $
// Revision 1.1  2008-12-12 19:54:06  thernis
// Implement UV emission raytracing
//
