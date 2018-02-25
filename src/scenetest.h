
#ifndef SCENETEST_H
#define SCENETEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "scene.h"

//! Unit tests for the Scene class
class SceneTest : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (SceneTest);
  //CPPUNIT_TEST (testScene);
  //CPPUNIT_TEST (testSceneWithPhysics);
  CPPUNIT_TEST (testSceneWTF);
  CPPUNIT_TEST_SUITE_END ();

public:
	void setUp (void);
	void tearDown (void);

protected:
	void testScene (void);
    void testSceneWithPhysics (void);
    void testSceneWTF (void);
  
    Scene* pscene;
    float *btot,*bpol,*netot;
    unsigned int sx,sy;


};

#endif
