
#ifndef STDRTMISCFUNCTEST_H
#define STDRTMISCFUNCTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#define NBSAMP 3

//! Unit tests of the functions in rtmiscfunc.h
class rtmiscfunctest : public CPPUNIT_NS::TestFixture
{

  CPPUNIT_TEST_SUITE (rtmiscfunctest);
  CPPUNIT_TEST (testrtmiscfunc);
  CPPUNIT_TEST_SUITE_END ();

public:
  void setUp (void);
	void tearDown (void);
protected:
	void testrtmiscfunc (void);
private:
  static const float xarray[NBSAMP];
  static const float yarray[NBSAMP];
  static const unsigned int nbsamp;

};
const float rtmiscfunctest::xarray[NBSAMP]={1,2,3};
const float rtmiscfunctest::yarray[NBSAMP]={-2,-1,0};
const unsigned int rtmiscfunctest::nbsamp=NBSAMP;

#endif
