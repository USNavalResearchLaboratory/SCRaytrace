
#include <string>
#include "scenetest.h"
#include "scene.h"


CPPUNIT_TEST_SUITE_REGISTRATION (SceneTest);


void SceneTest :: setUp (void)
{
    pscene=new Scene;
    sx=256;sy=256;
    btot=new float[sx*sy];
    bpol=new float[sx*sy];
    netot=new float[sx*sy];

}

void SceneTest :: tearDown (void)
{
  delete pscene;
  delete [] btot,bpol,netot;
}

void SceneTest :: testScene (void)
{


  unsigned int sx=256,sy=256;
  Scene scene;
  scene.camera.setCCD(CCD(sx,sy));
  scene.camera.setProjType(ARC);
  scene.camera.setFovpix(0.01);
  scene.camera.setCrpix(63.5,63.5);

  // -- LOS integration parameter definition
  int losnbp=20;
  float losrange[2]={190,230};
  scene.los.setLOS(losnbp,losrange[0],losrange[1]);

  // -- output images
  float *btot=new float[sx*sy];
  float *bpol=new float[sx*sy];
  float *netot=new float[sx*sy];

  // -- density model
  int modelid=14;
  float *pmodparam;
  scene.setDensityModel(14,pmodparam);
  
  // -- position of camera and density model
  scene.setobs(Cbasis());	

  // ---- run raytracing
  scene.computeImagebyRay(btot,bpol,netot,2);


  //std::cout << scene.getPhysics() << std::endl;
  //CPPUNIT_ASSERT(string("Thomson Scattering").compare(scene.getPhysics()) ==0);

    delete [] btot,bpol,netot;

  std::cout << "Yo ! " << std::endl;

}





void SceneTest :: testSceneWTF (void)
{
    string st;
    st=pscene->getPhysics();
    std::cout << st << std::endl;
    CPPUNIT_ASSERT(string("Thomson Scattering").compare(st) ==0);
    
    pscene->setPhysics(UV);
    st=pscene->getPhysics();
    std::cout << st << std::endl;
    CPPUNIT_ASSERT(string("UV emission").compare(st) ==0);
    
    pscene->camera.setCCD(CCD(sx,sy));
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
    pscene->computeImagebyChunk(btot,bpol,netot,8,4);
    
    CPPUNIT_ASSERT(string("UV emission").compare(pscene->getPhysics()) ==0);

    
}



