// $Id: suntest.h,v 1.1 2009-02-09 20:46:46 thernis Exp $
#ifndef SUNTEST_H
#define SUNTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "sun.h"

//! Unit tests for the Sun class
class SunTest : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (SunTest);
  CPPUNIT_TEST (test);
  CPPUNIT_TEST (testgetThomsonCoeff);
  CPPUNIT_TEST (testgetPosition);
  CPPUNIT_TEST_SUITE_END ();

  public:
    void setUp (void);
    void tearDown (void);

  protected:
    void test (void);
    void testgetThomsonCoeff (void);
    void testgetPosition (void);

  private:
    Sun *a, *b, *c;

};

#endif
