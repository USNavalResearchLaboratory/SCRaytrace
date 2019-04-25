
#ifndef MODEL77TEST_H
#define MODEL77TEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "models71to80.h"


//! Unit tests of CModel77
class model77Test : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (model77Test);
  CPPUNIT_TEST (testmodel77);
  CPPUNIT_TEST_SUITE_END ();

public:
  void setUp (void);
    void tearDown (void);

protected:
    void testmodel77 (void);

private:

  CModel77 *silicate, *carbon;

  float pparam[2];


};




#endif
