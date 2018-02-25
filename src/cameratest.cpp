
#include <iostream>
#include "cameratest.h"
#include "Cvec.h"
#include "raytrace.h"


CPPUNIT_TEST_SUITE_REGISTRATION (CameraTest);

void CameraTest :: setUp (void)
{
    // set up test environment (initializing objects)
  a = new Camera;
  b = new Camera;
  c = new Camera;
  c1 = new CCD;
  c2 = new CCD;
}

void CameraTest :: tearDown (void)
{
    // finally delete objects
  delete a; delete b; delete c;
  delete c1; delete c2;
}

void CameraTest :: testCCD (void)
{
  CPPUNIT_ASSERT_EQUAL ((unsigned int)0,c1->getSizePixX());
  CPPUNIT_ASSERT_EQUAL ((unsigned int)0,c1->getSizePixY());
  CPPUNIT_ASSERT_EQUAL ((float)1.,c1->getSizemmX());
  CPPUNIT_ASSERT_EQUAL ((float)1.,c1->getSizemmY());

  CCD c3(1024,1024),c4;
  CPPUNIT_ASSERT_EQUAL ((unsigned int)1024,c3.getSizePixX());
  CPPUNIT_ASSERT_EQUAL ((unsigned int)1024,c3.getSizePixY());

  c3.setSizemm(10.,10.);
  CPPUNIT_ASSERT_EQUAL ((float)10.,c3.getSizemmX());
  CPPUNIT_ASSERT_EQUAL ((float)10.,c3.getSizemmY());

  
  c4=c3;
  CPPUNIT_ASSERT(c4==c3);
  CPPUNIT_ASSERT_EQUAL ((unsigned int)1024,c4.getSizePixX());
  CPPUNIT_ASSERT_EQUAL ((unsigned int)1024,c4.getSizePixY());
  CPPUNIT_ASSERT_EQUAL ((float)10.,c4.getSizemmX());
  CPPUNIT_ASSERT_EQUAL ((float)10.,c4.getSizemmY());

  c2->setSizePix(2048,2048);
  CPPUNIT_ASSERT_EQUAL ((unsigned int)2048,c2->getSizePixX());
  CPPUNIT_ASSERT_EQUAL ((unsigned int)2048,c2->getSizePixY());

}

void CameraTest :: testCamera (void)
{
  CPPUNIT_ASSERT_EQUAL ((ProjType)ARC,a->getProjType());
  CPPUNIT_ASSERT_EQUAL ((float)0.,a->getFovpix());
  CPPUNIT_ASSERT (CCD()==a->getCCD());

  b->setFovpix(10.);
  CPPUNIT_ASSERT_EQUAL((float)10.,b->getFovpix());
  
  b->setProjType(TAN);
  CPPUNIT_ASSERT_EQUAL((ProjType)TAN,b->getProjType());

  Camera cam(0.01,ARC,CCD(512,512),256.,255.,1.5);
  CPPUNIT_ASSERT_EQUAL((float)0.01,cam.getFovpix());
  CPPUNIT_ASSERT_EQUAL((ProjType)ARC,cam.getProjType());
  CPPUNIT_ASSERT_EQUAL((float)256,cam.getCrpix1());
  CPPUNIT_ASSERT_EQUAL((float)255,cam.getCrpix2());
  CPPUNIT_ASSERT_EQUAL((float)1.5,cam.getPv2_1());
  CPPUNIT_ASSERT_EQUAL((float)1.,cam.getPc(0));
  CPPUNIT_ASSERT_EQUAL((float)0,cam.getPc(1));
  float pc[4]={1,2,3,4};
  cam.setPc(pc);
  CPPUNIT_ASSERT_EQUAL((float)3.,cam.getPc(2));
  CPPUNIT_ASSERT_EQUAL((float)4.,cam.getPc(3));


  float crpix[2]={256.,255.};
	float fovpix=0.01;
	float *fp;
	fp=&fovpix;

//	float *pc2=new float[4];
	pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;

  Camera cam2(0.01,ARC,CCD(512,512),256.,255.,0.);
	Cvec vlosobs0=cam2.ij2los(256.,255.);
	CPPUNIT_ASSERT_EQUAL(Cvec(0,0,1),vlosobs0);

	float iii=270.,jjj=240.;
	vlosobs0=cam2.ij2los(iii,jjj);
	std::cout << "vlosobs0 : " << vlosobs0 << std::endl;

  Cvec vlosabs,vlosobs;
  float rho,rrr,alpha;
  //vlosobs=xxyy2losobs(iii,jjj,rrr,alpha,crpix,fovpix,pc);

	std::cout << "*fp : " << *fp << std::endl;	

	Cbasis obs=Cbasis(Cvec(0,0,-214),0.01,0,0);
 	Cbasis abs=Cbasis(Cvec(0,0,0),0,0,0);
  getimpactandlos(iii,jjj,obs,abs,1,0.,vlosabs,rho,rrr,alpha,vlosobs,crpix,fp,pc);

	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosobs[0],vlosobs0[0],1e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosobs[1],vlosobs0[1],1e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosobs[2],vlosobs0[2],1e-5);

	std::cout << "vlosabs : " << vlosabs << std::endl;
	std::cout << "vlosobs : " << vlosobs << std::endl;
	std::cout << "rho : " << rho << std::endl;
	std::cout << "rrr : " << rrr << std::endl;

	Cvec vlosabs0=obs.ui * vlosobs0;
	std::cout << "vlosabs0 : " << vlosabs0 << std::endl;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosabs[0],vlosabs0[0],1e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosabs[1],vlosabs0[1],1e-5);
	CPPUNIT_ASSERT_DOUBLES_EQUAL(vlosabs[2],vlosabs0[2],1e-5);

  float rho0=psldist(obs.o,vlosabs0,abs.o);
	std::cout << "rho0 : " << rho0 << std::endl;
	CPPUNIT_ASSERT_DOUBLES_EQUAL(rho,rho0,1e-5);


}

