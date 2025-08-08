
#include <string>
#include "scene.h"




#define BOOST_TEST_MODULE SceneTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


#define NBSAMP 3

struct SceneTest {

    Scene* pscene;
    float *btot,*bpol,*netot;
    unsigned int sx,sy;

    
SceneTest() { 
    pscene=new Scene;
    sx=256;sy=256;
    btot=new float[sx*sy];
    bpol=new float[sx*sy];
    netot=new float[sx*sy];

}

  ~SceneTest() { 
  delete pscene;
  delete[] btot;
  delete[] bpol;
  delete[] netot;

}

};
                 


BOOST_FIXTURE_TEST_SUITE(s, SceneTest)

  BOOST_AUTO_TEST_CASE(test_Cvec)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running SceneTest");
     
   
    string st;
    st=pscene->getPhysics();
    std::cout << st << std::endl;
    BOOST_TEST(string("Thomson Scattering").compare(st) == 0);
    
//     pscene->setPhysics(UV);
//     st=pscene->getPhysics();
//     std::cout << st << std::endl;
//     BOOST_TEST(string("UV emission").compare(st) == 0);
    
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
    
    float phyparam=1;
    pscene->setPhysicsParam(&phyparam);
    
    pscene->computeImagebyRay(btot,bpol,netot,2);
    pscene->computeImagebyChunk(btot,bpol,netot, 8, 8);
    
//     BOOST_TEST(string("UV emission").compare(pscene->getPhysics()) == 0);


}

BOOST_AUTO_TEST_SUITE_END()

















