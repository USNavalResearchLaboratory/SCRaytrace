
#include <iostream>
#include <string>
#include "Cvec.h"

#include "scene.h"
#include "CModelBase.h"

#include "physicsbase.h"
#include "physicsthomson.h"
#include "physicsuv.h"
#include "physicsvsf.h"
#include "physicsvsfvarydist.h"


#define BOOST_TEST_MODULE PhysicsTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


struct PhysicsTest {
   PhysicsBase base;

PhysicsTest() { 

}

  ~PhysicsTest() { 

}

};

                  
BOOST_FIXTURE_TEST_SUITE(s, PhysicsTest)

  BOOST_AUTO_TEST_CASE(test_PhysicsTest)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running PhysicsTest");
   
    string s1;
    bool flagok;
    float bto,bpo,neo;
    Scene *pscene;
    pscene = new Scene;
    
    base.setParentScene(pscene);
    
    unsigned int sx=256,sy=256;
    pscene->camera.setDetector(Detector(sx,sy));
    pscene->camera.setProjType(ARC);
    pscene->camera.setFovpix(0.01);
    pscene->camera.setCrpix(63.5,63.5);

    // -- LOS integration parameter definition
    
    int losnbp=20;
    float losrange[2]={190,230};
    pscene->los.setLOS(losnbp,losrange[0],losrange[1]);

    // -- density model
    
    int modelid=14;
    float *pmodparam;
    pmodparam=NULL;
    pscene->setDensityModel(14,pmodparam);
    pscene->setobs(Cbasis());
    
    s1=base.getPhysics();
    BOOST_TEST(string("Physics Base").compare(s1) == 0);
    
    flagok=base.computeRadiation(Cvec(0,0,0),1,1,bto,bpo,neo);
    BOOST_TEST(1 == flagok);


   PhysicsBase *puv;
   puv = physicsSelect(UV);
   if (puv) {
      cout<< "UV physics selected" << endl;
   }
    puv->printParam();
    float phyparam=1;
    puv->setParam(&phyparam);
    puv->printParam();

   if (puv) {delete puv;}

   if (pscene) delete pscene;
    
    // -- test physicsVSF constructor
    PhysicsVSF vsf;

    // -- test physicsVSFVaryDist constructor
    PhysicsVSFVaryDist vsfvd;
  
}

BOOST_AUTO_TEST_SUITE_END()

