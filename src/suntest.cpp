
#include "suntest.h"
#include "Cvec.h"
#include "Cbasis.h"


CPPUNIT_TEST_SUITE_REGISTRATION (SunTest);

void SunTest :: setUp (void)
{
    // set up test environment (initializing objects)
  a = new Sun;
  b = new Sun;
  c = new Sun;
}

void SunTest :: tearDown (void)
{
    // finally delete objects
  delete a; delete b; delete c;
}

void SunTest :: test (void)
{
  CPPUNIT_ASSERT_EQUAL (a->getLimbDarkening(), float(0.7));
  a->setLimbDarkening(0.8);
  CPPUNIT_ASSERT_EQUAL (a->getLimbDarkening(), float(0.8));
  CPPUNIT_ASSERT_EQUAL (a->getRadius(), float(696000e3));
}

void SunTest :: testgetThomsonCoeff (void)
{
  float btotc,bpolc;
  b->getThomsonCoeff(3.,3.,btotc,bpolc);

  CPPUNIT_ASSERT_DOUBLES_EQUAL (btotc, float(0.0895308),1e6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL (bpolc, float(0.0809165),1e6);
}

void SunTest :: testgetPosition (void)
{
  Cbasis *base;
  base=c->getPosition();

  CPPUNIT_ASSERT_EQUAL (Cvec(0.,0.,0.),base->o);
  
}
