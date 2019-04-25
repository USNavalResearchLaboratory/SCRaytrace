/** \file cameratest.h
 * \brief Test suite for the camera class.
 */
 
#ifndef CAMERATEST_H
#define CAMERATEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "camera.h"

//! Unit test for the Camera class
class CameraTest : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (CameraTest);
  CPPUNIT_TEST (testCCD);
  CPPUNIT_TEST (testCamera);
  CPPUNIT_TEST_SUITE_END ();

  public:
    void setUp (void);
    void tearDown (void);

  protected:
    void testCCD (void);
    void testCamera (void);

  private:
    Camera *a, *b, *c;
    CCD *c1,*c2;

};


#endif
