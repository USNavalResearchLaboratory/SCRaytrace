/*! \file Cbasistest.h 
 * \brief Test suite for Cbasis.
 *
 */

#ifndef CBASISTEST_H
#define CBASISTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "Cbasis.h"

//! Unit test for the Cbasis class
class CbasisTest : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (CbasisTest);
  CPPUNIT_TEST (testCbasis);
  CPPUNIT_TEST_SUITE_END ();

  public:
    void setUp (void);
    void tearDown (void);

  protected:
    void testCbasis (void);

  private:
    Cbasis *a, *b, *c;
   
};


#endif
