
#include "Cvec.h"


#define BOOST_TEST_MODULE CvecTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


struct CvecTest {
    Cvec *a, *b, *c;

CvecTest() { 
  a = new Cvec (1, 2, 3);
  b = new Cvec (2, 3, 4);
  c = new Cvec (1, 2, 2);

}

  ~CvecTest() { 
 delete a; 
 delete b;}

};

                  
BOOST_FIXTURE_TEST_SUITE(s, CvecTest)

  BOOST_AUTO_TEST_CASE(test_Cvec)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running CvecTest");
     
   
     BOOST_TEST((*a + *b) == Cvec(3,5,7));
  BOOST_TEST(c->norm() == 3.);

     
}

BOOST_AUTO_TEST_SUITE_END()

