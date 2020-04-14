
#include <iostream>
#include "Cbasis.h"
#include "Cvec.h"
#include "constant.h"

#define BOOST_TEST_MODULE CbasisTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


struct CbasisTest {

  Cbasis *a, *b, *c;

  CbasisTest() { 
    a = new Cbasis;
    b = new Cbasis;
    c = new Cbasis;

}
  ~CbasisTest() { 
    delete a; 
    delete b; 
    delete c;

}
};

BOOST_FIXTURE_TEST_SUITE(s, CbasisTest)

  BOOST_AUTO_TEST_CASE(test_Cbasis)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running test_Cbasis");

   BOOST_REQUIRE_EQUAL(Cvec(0,0,0), a->o);
   BOOST_REQUIRE_EQUAL(Cmat(1, 0, 0, 0, 1, 0, 0, 0, 1), a->u);
   

  Cbasis e1(Cvec(0,0,0),0,0,0,0,0,0);
  Cbasis e2(Cvec(1,1,0),0,0,0,0,0,0);
  Cbasis e3(Cvec(1,1,0),0,0,PI/2,0,0,0);

  Cvec v1_1(1,0,0);
  Cvec v1_2 = ChangeBase(v1_1, e1, e2);
  BOOST_REQUIRE_EQUAL(Cvec(0,-1,0), v1_2);

  Cvec v1_3 = ChangeBase(v1_1, e1, e3);
  
  BOOST_TEST(-1. == v1_3[0], tt::tolerance(0.0001));
  BOOST_TEST( 0. == v1_3[1], tt::tolerance(0.0001));
  BOOST_TEST( 0. == v1_3[2], tt::tolerance(0.0001));

  Cvec cntr(1,2,3);
  e1.setCenter(cntr);
  BOOST_REQUIRE_EQUAL(cntr, e1.o);

  e1.setRotationPerAxis(0., 1, PI/2., 2, 0., 3);
  
  BOOST_TEST( 0. == e1.u[0][0], tt::tolerance(0.0001));
  BOOST_TEST( 0. == e1.u[0][1], tt::tolerance(0.0001));
  BOOST_TEST( 1. == e1.u[0][2], tt::tolerance(0.0001));

  BOOST_TEST( 0. == e1.u[1][0], tt::tolerance(0.0001));
  BOOST_TEST( 1. == e1.u[1][1], tt::tolerance(0.0001));
  BOOST_TEST( 0. == e1.u[1][2], tt::tolerance(0.0001));

  BOOST_TEST( 0. == e1.ui[0][0], tt::tolerance(0.0001));
  BOOST_TEST( 0. == e1.ui[0][1], tt::tolerance(0.0001));
  BOOST_TEST(-1. == e1.ui[0][2], tt::tolerance(0.0001));

}


BOOST_AUTO_TEST_SUITE_END()





