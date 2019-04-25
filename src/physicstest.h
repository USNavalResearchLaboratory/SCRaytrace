/** \file physicstest.h
 * \brief Implement test suite for various physics.
 */

#ifndef PHYSICSTEST_H
#define PHYSICSTEST_H


#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include "physicsbase.h"
#include "physicsthomson.h"
#include "physicsuv.h"
#include "physicsvsf.h"
#include "physicsvsfvarydist.h"


//! Unit tests for the Physics classes
class PhysicsTest : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (PhysicsTest);
  CPPUNIT_TEST (generalTests);
  CPPUNIT_TEST_SUITE_END ();

  public:
    void setUp (void);
    void tearDown (void);

  protected:
    void generalTests (void);

  private:
   PhysicsBase base;
   PhysicsBase *pthom,*puv;

};

#endif
