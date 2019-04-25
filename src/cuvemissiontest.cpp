/*! \file cuvemissiontest.cpp
 * \brief Test of the UV Emission physics
 *
 *  
 */




#include "cuvemissiontest.h"

CPPUNIT_TEST_SUITE_REGISTRATION (CUVEmissionTest);


void CUVEmissionTest :: setUp (void)
{
	// set up test environment (initializing objects)
	a=new CUVEmission;

}

void CUVEmissionTest :: tearDown (void)
{
	// finally delete objects
	delete a; 
}


void CUVEmissionTest :: testCUVEmission (void)
{
CPPUNIT_ASSERT(a->IsGood());
CPPUNIT_ASSERT_DOUBLES_EQUAL(log10(0.997800+1e-30),(double)a->getyisel(0,0),1e-8);
CPPUNIT_ASSERT_DOUBLES_EQUAL(log10(0.0116900+1e-30),(double)a->getyisel(30,2),1e-4);

//CPPUNIT_ASSERT_DOUBLES_EQUAL(-7.26769,(double)a->calcEmissivity(1,1e5),1e-3); // test on yif
//CPPUNIT_ASSERT_DOUBLES_EQUAL(-12.3636,(double)a->calcEmissivity(1,1e5),1e-3); // test on emiss
CPPUNIT_ASSERT_DOUBLES_EQUAL(2.33728e-20,(double)a->calcEmissivity(2,1e5),1e-24);


}


