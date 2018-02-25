
#include <iostream>
#include "Cbasistest.h"
#include "constant.h"

CPPUNIT_TEST_SUITE_REGISTRATION (CbasisTest);

void CbasisTest :: setUp (void)
{
    // set up test environment (initializing objects)
  a = new Cbasis;
  b = new Cbasis;
  c = new Cbasis;
}

void CbasisTest :: tearDown (void)
{
    // finally delete objects
  delete a; delete b; delete c;
}


void CbasisTest :: testCbasis (void)
{

  CPPUNIT_ASSERT_EQUAL (Cvec(0,0,0),a->o);
 // CPPUNIT_ASSERT_EQUAL (Cvec(0,0,0),a->translation);
  CPPUNIT_ASSERT_EQUAL (Cmat(1,0,0,0,1,0,0,0,1),a->u);

  //Cbasis d(Cvec(1,2,3),0,0,0,0,0,0,Cvec(5,4,3));
	//CPPUNIT_ASSERT_EQUAL (Cvec(1,2,3),d.o);
	//CPPUNIT_ASSERT_EQUAL (Cvec(5,4,3),d.translation);


  Cbasis e1(Cvec(0,0,0),0,0,0,0,0,0);
	Cbasis e2(Cvec(1,1,0),0,0,0,0,0,0);
	Cbasis e3(Cvec(1,1,0),0,0,PI/2,0,0,0);
	//Cbasis e4(Cvec(1,1,0),0,0,0,0,0,0,Cvec(0,0,1));

  Cvec v1_1(1,0,0);
  Cvec v1_2=ChangeBase(v1_1,e1,e2);
  CPPUNIT_ASSERT_EQUAL (Cvec(0,-1,0),v1_2);

  Cvec v1_3=ChangeBase(v1_1,e1,e3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL (-1,v1_3[0],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,v1_3[1],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,v1_3[2],1e-5);

  Cvec cntr(1,2,3);
  e1.setCenter(cntr);
  CPPUNIT_ASSERT_EQUAL (cntr,e1.o);


  e1.setRotationPerAxis(0.,1,PI/2.,2,0.,3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,e1.u[0][0],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,e1.u[0][1],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1,e1.u[0][2],1e-5);

  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,e1.u[1][0],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1,e1.u[1][1],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,e1.u[1][2],1e-5);

  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,e1.ui[0][0],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL ( 0,e1.ui[0][1],1e-5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL (-1,e1.ui[0][2],1e-5);



	//Cvec v2_4(1,0,0);
  //Cvec v2_1=ChangeBase(v2_4,e4,e1);
  //CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2,v2_1[0],1e-5);
  //CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1,v2_1[1],1e-5);
  //CPPUNIT_ASSERT_DOUBLES_EQUAL ( 1,v2_1[2],1e-5);




// 	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosabs[0],vlosabs0[0],1e-5);
// 	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosabs[1],vlosabs0[1],1e-5);
// 	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosabs[2],vlosabs0[2],1e-5);


}

