
#include "constant.h"
#include "Cvec.h"
#include "model77test.h"

CPPUNIT_TEST_SUITE_REGISTRATION (model77Test);


void model77Test :: setUp (void)
{
    // ---- create a model instance
    silicate = new CModel77;
    carbon = new CModel77;
    
    pparam[0] = 1.;
    pparam[1] = 0.; 
    
    silicate->initParam(pparam);

    
    pparam[1] = 1.;
    carbon->initParam(pparam);
    
    
}

void model77Test :: tearDown (void)
{
    // finally delete objects
    delete silicate;
    delete carbon;

}

void model77Test :: testmodel77 (void)
{
  
  silicate->checkData();
  float ef;
  ef = silicate->getEnhanceFactor(5.);
  
  std::cout << "ef = silicate->getEnhanceFactor(5.) : " << ef << std::endl;
  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)1.02351,(double)ef,1e-3);

  carbon->checkData();
  ef = carbon->getEnhanceFactor(4.00108);
  std::cout << "ef = carbon->getEnhanceFactor(4.00108) : " << ef << std::endl;
  CPPUNIT_ASSERT_DOUBLES_EQUAL((double)1.28268,(double)ef,1e-3);

  
  
/*
float dens=-1.;
Cvec pos(0.,0.,0.);
dens=a->Density(pos);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)3.,(double)dens,1e-3);


pos=Cvec(0.,0.,1.);
dens=a->Density(pos);
cout << "dens : " << dens << endl;
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)2.,(double)dens,1e-3);
*/

}
