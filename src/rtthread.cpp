
#include <iostream>
#include <cstdlib>
#include "rtthread.h"
#include "camera.h"
#include "Cbasis.h"
#include "scene.h"
#include "physicsbase.h"



extern "C" int rtthread(int sx,
                        int sy,
                        float fovpix,
                        float *obspos,
                        float *obsang,
                        float *nepos,
                        float *neang,
                        int losnbp,
                        float *losrange,
                        int modelid,
                        float *btot,
                        float *bpol,
                        float *netot,
                        float *pmodparam,
                        float *crpix,
                        int quiet,
                        int neonly,
                        float *hlonlat,
                        float occrad,
                        float limbdark,
                        float *obslonlat,
                        int obslonlatflag,
                        unsigned int projtypecode,
                        float pv2_1,
                        float *pc,
                        int frontinteg,
                        unsigned int nbthreads,
                        unsigned int nbchunk,
                        float *nerotcntr,
                        float *nerotang,
                        float *netranslation,
                        int *nerotaxis,
                        int physics,
                        float *phyparam,
                        float fracmax,
                        int runDumpInteg, 
                        float *pIntegrand,
                        int npv,
                        float *pv,
                        int *pv_i,
                        int *pv_m)
{
  if (quiet!=1) {
  std::cout << "In rtthread..." << std::endl;
	printvar(sx);
	printvar(sy);
	printvar(fovpix);
	printvarnoend(obspos[0]);printvarnoend(obspos[1]);printvar(obspos[2]);
	printvarnoend(obsang[0]);printvarnoend(obsang[1]);printvar(obsang[2]);
	printvarnoend(nepos[0]);printvarnoend(nepos[1]);printvar(nepos[2]);
	printvarnoend(neang[0]);printvarnoend(neang[1]);printvar(neang[2]);
	printvarnoend(nerotcntr[0]);printvarnoend(nerotcntr[1]);printvar(nerotcntr[2]);
	printvarnoend(nerotang[0]);printvarnoend(nerotang[1]);printvar(nerotang[2]);
	printvarnoend(nerotaxis[0]);printvarnoend(nerotaxis[1]);printvar(nerotaxis[2]);
	printvarnoend(netranslation[0]);printvarnoend(netranslation[1]);printvar(netranslation[2]);
	printvar(losnbp);
	printvarnoend(losrange[0]);printvar(losrange[1]);
	printvar(modelid);
    printvarnoend(crpix[0]);printvar(crpix[1]);
    printvar(quiet);
	printvar(neonly);
	printvarnoend(hlonlat[0]);printvarnoend(hlonlat[1]);printvar(hlonlat[2]);
	printvar(occrad);
// 	printvar(limbdark);
	printvarnoend(obslonlat[0]);printvarnoend(obslonlat[1]);printvar(obslonlat[2]);
	printvar(obslonlatflag);
	printvar(projtypecode);
	printvar(pv2_1);
	printvarnoend(pc[0]);printvar(pc[1]);
	printvarnoend(pc[2]);printvar(pc[3]);
	printvar(frontinteg);
	printvar(nbthreads);
	printvar(nbchunk);
    printvar(physics);
    printvar(fracmax);
    printvar(runDumpInteg);
    cout << "Compilation date : " << __DATE__ << " " << __TIME__ << endl; 
  printvar(npv);
  int i;
  for (i=0; i<npv; i++){
    printvarnoend(pv_i[i]);printvarnoend(pv_m[i]);printvar(pv[i]);
  }
}
  // ---- setup the scene
  // -- Camera parameters
  Scene scene;
  scene.camera.setDetector(Detector(sx,sy));
  scene.camera.setProjType(projtypecode);
  scene.camera.setFovpix(fovpix);
  scene.camera.setCrpix(crpix[0],crpix[1]);
  scene.camera.setPv2_1(pv2_1);
  scene.camera.setPc(pc);
  scene.camera.setPv(npv, pv, pv_i, pv_m);
//   scene.csun.setLimbDarkening(limbdark);
  scene.setNeonly(neonly);
  scene.setFrontInteg(frontinteg);
  scene.setQuiet(quiet);
  scene.setFracMax(fracmax);
  scene.setDumpIntegrand((runDumpInteg == 1),pIntegrand);
  
  // -- LOS integration parameter definition
  scene.los.setLOS(losnbp,losrange[0],losrange[1]);

  // -- density model
  scene.setDensityModel(modelid,pmodparam);
  // -- physics
  scene.setPhysics((PhysicsType)physics);
  scene.setPhysicsParam(phyparam);
  std::cout<<scene.getPhysics()<<std::endl;
  scene.printPhysicsParam();

  // -- position of camera and density model
  if (obslonlatflag == 0) 
    scene.setobs(Cbasis(Cvec(obspos[0],obspos[1],obspos[2]),obsang[0],obsang[1],obsang[2]));
  else 
    scene.setobs(Cbasis(obslonlat[0],obslonlat[1],obslonlat[2],obsang[0],obsang[1],obsang[2]));
   // ---- Ne position
  scene.modelposition.setBasis(Cbasis(Cvec(nepos[0],nepos[1],nepos[2]),neang[0],neang[1],neang[2],hlonlat[0],hlonlat[1],hlonlat[2]));
  scene.modelposition.rotation.setCenter(Cvec(nerotcntr[0],nerotcntr[1],nerotcntr[2]));
  scene.modelposition.rotation.setRotationPerAxis(nerotang[0],nerotaxis[0],nerotang[1],nerotaxis[1],nerotang[2],nerotaxis[2]);

  scene.modelposition.setTranslation(Cvec(netranslation[0],netranslation[1],netranslation[2]));


  // ---- run raytracing
  if (nbchunk < 2) 
    scene.computeImagebyRay(btot,bpol,netot,nbthreads); 
  else
    scene.computeImagebyChunk(btot,bpol,netot,nbthreads,nbchunk);

  if (quiet != 1) std::cout << "Done ! " << std::endl;

  return EXIT_SUCCESS;

}



int rtthreadtest()
{
  // ---- setup the scene
  // -- Camera parameters
  unsigned int sx=256,sy=256;
  Scene scene;
  scene.camera.setDetector(Detector(sx,sy));
  scene.camera.setProjType(ARC);
  scene.camera.setFovpix(0.01);
  scene.camera.setCrpix(63.5,63.5);

  // -- LOS integration parameter definition
  int losnbp=200;
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
//  scene.setnps(Cbasis());	


  // ---- run raytracing
  scene.computeImagebyRay(btot,bpol,netot,2);

    delete[] btot;
    delete[] bpol;
    delete[] netot;

  std::cout << "Yo ! " << std::endl;

  return EXIT_SUCCESS;
}


extern "C" void footest()
{
    std::cout << "In footest ! " << std::endl;
    
}


extern "C" void testPassPython(float* array, int N, float f)
{
    std::cout << "In testPassPython ! " << f << std::endl;
    for (int i=0; i<N; i++) 
        std::cout << i << " " << array[i] << std::endl;
    
}



extern "C" int rtthreadtestExt()
{
 
    int foo;
    
    foo = rtthreadtest();
    
    // Try calling an extern C function
    footest();
    
    
    return EXIT_SUCCESS;
}



//! Raytracing with threads to speed up on multi-core computers
extern "C" int rtthreadidl(int argc, void **argv)
{
  int sx=*((int*) argv[0]);
  int sy=*((int*) argv[1]);
  float fovpix=*((float*) argv[2]);
  float *obspos=(float*) argv[3];
  float *obsang=(float*) argv[4];
  float *nepos=(float*) argv[5];
  float *neang=(float*) argv[6];
  int losnbp=*((int*) argv[7]);
  float *losrange=(float*) argv[8];
  int modelid=*((int*) argv[9]);
  float *btot=(float*) argv[10];
  float *bpol=(float*) argv[11];
  float *netot=(float*) argv[12];
  float *pmodparam=(float*) argv[13];
  float *crpix=(float*) argv[14];
  int quiet=*((int*) argv[15]);
  int neonly=*((int*) argv[16]);
  float *hlonlat=(float*) argv[17];
  float occrad=*((float*) argv[18]);
  float limbdark=*((float*) argv[19]);
  float *obslonlat=(float*) argv[20];
  int obslonlatflag=*((int*) argv[21]);
  unsigned int projtypecode=*((int*) argv[22]);
  float pv2_1=*((float*) argv[23]);
  float *pc=(float*) argv[24];
  int frontinteg=*((int*) argv[25]);
  unsigned int nbthreads=*((unsigned int*) argv[26]);
  unsigned int nbchunk=*((unsigned int*) argv[27]);
  float *nerotcntr=(float*) argv[28];
  float *nerotang=(float*) argv[29];
  float *netranslation=(float*) argv[30];
  int *nerotaxis=(int*) argv[31];
  int physics=*((int*) argv[32]);
  float *phyparam=(float*) argv[33];
  float fracmax=*((float*) argv[34]);
  int runDumpInteg = *((int*) argv[35]);
  float *pIntegrand=(float*) argv[36];
  int npv=*((int*) argv[37]);
  float *pv=(float*) argv[38];
  int *pv_i=(int*) argv[39];
  int *pv_m=(int*) argv[40];


	int out=rtthread(sx,sy,
    fovpix,obspos,obsang,nepos,neang,
    losnbp,losrange,modelid,btot,bpol,netot,
    pmodparam,crpix,quiet,neonly,hlonlat,
    occrad,limbdark,obslonlat,obslonlatflag,
    projtypecode,pv2_1,pc,frontinteg,nbthreads,
    nbchunk,nerotcntr,nerotang,netranslation,
    nerotaxis,physics,phyparam,fracmax,
    runDumpInteg,pIntegrand, npv, pv, pv_i, pv_m);

	return EXIT_SUCCESS;

}


