/*! \file cvectest.h 
 * \brief Test of the Cvec class.
 *
 */

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
