
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
  
  BOOST_REQUIRE_EQUAL((unsigned int)128, c1.getSizePixX());
  BOOST_REQUIRE_EQUAL((unsigned int)128, c1.getSizePixY());
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
  BOOST_REQUIRE_CLOSE((float)0.1, a.getFovpix(), 0.001);
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

  static const int NPV=9;
  int pv_i[NPV] = {1,1,1,2,2,2,2,2,2};
  int pv_m[NPV] = {1,2,3,0,1,2,3,4,5};
  float pv[NPV] = {0.0, 90., 180.,
                      -1.98789997796E-08,
                        1.00501000881,
                        0.0729582980275,
                        0.275292009115,
                        -0.701880991459,
                        1.97518002987};

  Detector detB(2048, 1920, 20.48, 19.20);
  float fovpix_deg = 0.0211525;

  // Camera camB(fovpix_deg / RADEG,
  //             ZPN,
  //             detB,
  //             crpix[0], crpix[1],
  //             NPV,
  //             pv, pv_i, pv_m);

  // ---- Check between old computation of the los and the new which uses libwcs
  int pv_iC[3] = {1,1,1};
  int pv_mC[3] = {1,2,3};
  float pvC[3] = {0.0, 90., 180.};

  Detector detC(2048, 1920, 20.48, 19.20);
  // -- Test ARC projection
  BOOST_TEST_MESSAGE("Testing ARC projection");
  Camera camC(fovpix_deg / RADEG,
                ARC,
                detC,
                crpix[0], crpix[1],
                3,
                pvC, pv_iC, pv_mC);

  Cvec v_wcs, v_rt;
  v_wcs = camC.ij2loswcs(1024, 1024);
  v_rt = camC.ij2losold(1024, 1024);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  v_wcs = camC.ij2loswcs(991.547, 1015.11);
  v_rt = camC.ij2losold(991.547, 1015.11);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  v_wcs = camC.ij2loswcs(991.547, 1015.11 + 300);
  v_rt = camC.ij2losold(991.547, 1015.11 + 300);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  v_wcs = camC.ij2loswcs(991.547, 1015.11 - 300);
  v_rt = camC.ij2losold(991.547, 1015.11 - 300);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  v_wcs = camC.ij2loswcs(991.547 + 300, 1015.11);
  v_rt = camC.ij2losold(991.547 + 300, 1015.11);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  v_wcs = camC.ij2loswcs(991.547 - 300, 1015.11);
  v_rt = camC.ij2losold(991.547 - 300, 1015.11);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  // -- Test TAN projection
  BOOST_TEST_MESSAGE("Testing TAN projection");
  Camera camD(fovpix_deg / RADEG,
                TAN,
                detC,
                crpix[0], crpix[1],
                3,
                pvC, pv_iC, pv_mC);

  v_wcs = camD.ij2loswcs(1024, 1024);
  v_rt = camD.ij2losold(1024, 1024);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  // -- Test SIN projection
  BOOST_TEST_MESSAGE("Testing SIN projection");
  Camera camE(fovpix_deg / RADEG,
                SIN,
                detC,
                crpix[0], crpix[1],
                3,
                pvC, pv_iC, pv_mC);

  v_wcs = camE.ij2loswcs(1024, 1024);
  v_rt = camE.ij2losold(1024, 1024);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.001);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.001);

  // -- Test AZP projection
  BOOST_TEST_MESSAGE("Testing AZP projection");
  Camera camF(0.000369181, // fovpix, in rad/pix
                AZP,
                detC,
                991, 1015, // crpix
                0.3);  // pv2_1

  v_wcs = camF.ij2loswcs(1024, 1024);
  v_rt = camF.ij2losold(1024, 1024);
  BOOST_REQUIRE_CLOSE(v_wcs[0], v_rt[0], 0.01);
  BOOST_REQUIRE_CLOSE(v_wcs[1], v_rt[1], 0.01);
  BOOST_REQUIRE_CLOSE(v_wcs[2], v_rt[2], 0.01);

}




// BOOST_AUTO_TEST_CASE(test_cam_wcs_proj)
// {
//   BOOST_TEST_MESSAGE("running test_cam_wcs_proj");
 

//   float pc[4]={1,2,3,4};
//   pc[0]=1; pc[1]=0; pc[2]=0; pc[3]=1;

//   Camera cam2(0.01, ARC, Detector(512,512), 256., 255., 0.);
//   Cvec vlosobs0 = cam2.ij2los(256., 255.);
//   Cvec vlosobswcs = cam2.ij2loswcs(256., 255.);
//   BOOST_REQUIRE_EQUAL(Cvec(0,0,1), vlosobs0);
  
// //TODO: test error throw NPV_MAX


// }



