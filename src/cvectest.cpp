// $Id: cvectest.cpp,v 1.1 2009/02/09 20:46:26 thernis Exp $

#include "cvectest.h"
CPPUNIT_TEST_SUITE_REGISTRATION (CvecTest);

void CvecTest :: setUp (void)
{
  a = new Cvec (1, 2, 3);
  b = new Cvec (2, 3, 4);
  c = new Cvec (1, 2, 2);
}

void CvecTest :: tearDown (void)
{
  delete a; delete b;
}

void CvecTest :: Test (void)
{
  CPPUNIT_ASSERT_EQUAL (*a + *b, Cvec(3,5,7));
  CPPUNIT_ASSERT_EQUAL (c->norm(), float(3.));

}
