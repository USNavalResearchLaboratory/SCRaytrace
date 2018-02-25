// $Id: cvectest.h,v 1.1 2009/02/09 20:46:27 thernis Exp $

#ifndef CVECTEST_H
#define CVECTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "Cvec.h"

using namespace std;
//! Unit tests for the Cvec class
class CvecTest : public CPPUNIT_NS :: TestFixture
{
  CPPUNIT_TEST_SUITE (CvecTest);
  CPPUNIT_TEST (Test);
  CPPUNIT_TEST_SUITE_END ();

  public:
    void setUp (void);
    void tearDown (void);

  protected:
    void Test (void);

  private:
    Cvec *a, *b, *c;
};


#endif
