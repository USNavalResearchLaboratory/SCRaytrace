
#include "constant.h"
#include "Cvec.h"
#include "models71to80.h"


#define BOOST_TEST_MODULE model77Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


struct model77Test {

  CModel77 *silicate, *carbon;

  float pparam[2];

model77Test() { 
    silicate = new CModel77;
    carbon = new CModel77;
    
    pparam[0] = 1.;
    pparam[1] = 0.; 
    
    silicate->initParam(pparam);

    
    pparam[1] = 1.;
    carbon->initParam(pparam);

}

  ~model77Test() { 
    delete silicate;
    delete carbon;
}

};

                  
BOOST_FIXTURE_TEST_SUITE(s, model77Test)

  BOOST_AUTO_TEST_CASE(test_model77)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running model77Test");
     
  silicate->checkData();
  float ef;
  ef = silicate->getEnhanceFactor(5.);
  
  std::cout << "ef = silicate->getEnhanceFactor(5.) : " << ef << std::endl;
  BOOST_TEST(1.02351 == ef, tt::tolerance(0.001));

  carbon->checkData();
  ef = carbon->getEnhanceFactor(4.00108);
  std::cout << "ef = carbon->getEnhanceFactor(4.00108) : " << ef << std::endl;
  BOOST_TEST(1.28268 == ef, tt::tolerance(0.001));

   
}

BOOST_AUTO_TEST_SUITE_END()

