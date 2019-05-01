
#include <iostream>
#include "Cvec.h"
#include "camera.h"


#define BOOST_TEST_MODULE cameratest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


struct MyGlobalFixture {
  
  MyGlobalFixture() {
  }
  void setup() {
  }
  void teardown() {
  }
  ~MyGlobalFixture() {
  }
};



BOOST_TEST_GLOBAL_FIXTURE( MyGlobalFixture );

BOOST_AUTO_TEST_CASE(test_Detector)
{
  BOOST_TEST_MESSAGE("running test_Detector");
  
  Detector c1;
  Detector c2;
  
  BOOST_REQUIRE_EQUAL((unsigned int)0, c1.getSizePixX());
  BOOST_REQUIRE_EQUAL((unsigned int)0, c1.getSizePixX());
  BOOST_REQUIRE_CLOSE((float)1., c1.getSizemmX(), 0.001);
  BOOST_REQUIRE_CLOSE((float)1., c1.getSizemmY(), 0.001);

  Detector c3(1024, 1024), c4;
  BOOST_REQUIRE_EQUAL((unsigned int)1024, c3.getSizePixX());
  BOOST_REQUIRE_EQUAL((unsigned int)1024, c3.getSizePixY());

  c3.setSizemm(10., 10.);
  BOOST_REQUIRE_CLOSE((float)10., c3.getSizemmX(), 0.001);
  BOOST_REQUIRE_CLOSE((float)10., c3.getSizemmY(), 0.001);
  
  c4 = c3;
  bool foo;
  foo = (c4==c3);
  BOOST_TEST(foo);
  BOOST_REQUIRE_EQUAL((unsigned int)1024, c4.getSizePixX());
  BOOST_REQUIRE_EQUAL((unsigned int)1024, c4.getSizePixY());
  BOOST_REQUIRE_CLOSE((float)10., c4.getSizemmX(), 0.001);
  BOOST_REQUIRE_CLOSE((float)10., c4.getSizemmY(), 0.001);

  c2.setSizePix(2048, 2048);
  BOOST_REQUIRE_EQUAL((unsigned int)2048, c2.getSizePixX());
  BOOST_REQUIRE_EQUAL((unsigned int)2048, c2.getSizePixY());

}

BOOST_AUTO_TEST_CASE(test_camera)
{
  BOOST_TEST_MESSAGE("running test_camera");
  
  Camera a;
  
  BOOST_REQUIRE_EQUAL((ProjType)ARC, a.getProjType());
  BOOST_REQUIRE_CLOSE((float)0., a.getFovpix(), 0.001);
  bool foo;
  foo = Detector()==a.getDetector();
  BOOST_TEST(foo);
  
  Camera b;
  
  b.setFovpix(10.);
  BOOST_REQUIRE_CLOSE((float)10., b.getFovpix(), 0.001);
  
  b.setProjType(TAN);
  BOOST_REQUIRE_EQUAL((ProjType)TAN, b.getProjType());

  Camera cam(0.01, ARC, Detector(512,512), 256., 255., 1.5);
  BOOST_REQUIRE_CLOSE((float)0.01, cam.getFovpix(), 0.001);
  BOOST_REQUIRE_EQUAL((ProjType)ARC, cam.getProjType());
  BOOST_REQUIRE_CLOSE((float)256, cam.getCrpix1(), 0.001);
  BOOST_REQUIRE_CLOSE((float)255, cam.getCrpix2(), 0.001);
  BOOST_REQUIRE_CLOSE((float)1.5, cam.getPv2_1(), 0.001);
  BOOST_REQUIRE_CLOSE((float)1., cam.getPc(0), 0.001);
  BOOST_REQUIRE_CLOSE((float)0, cam.getPc(1), 0.001);
  float pc[4]={1,2,3,4};
  cam.setPc(pc);
  BOOST_REQUIRE_CLOSE((float)3., cam.getPc(2), 0.001);
  BOOST_REQUIRE_CLOSE((float)4., cam.getPc(3), 0.001);

  float crpix[2] = {256.,255.};
  float fovpix = 0.01;
  float *fp;
  fp = &fovpix;

  float *pc2 = new float[4];
  pc[0]=1; pc[1]=0; pc[2]=0; pc[3]=1;

  Camera cam2(0.01, ARC, Detector(512,512), 256., 255., 0.);
  Cvec vlosobs0 = cam2.ij2los(256., 255.);
  BOOST_REQUIRE_EQUAL(Cvec(0,0,1), vlosobs0);

  delete[] pc2;
}

