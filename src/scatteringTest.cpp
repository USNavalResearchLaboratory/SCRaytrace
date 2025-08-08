
#include <iostream>
#include <string>
#include "Cvec.h"
#include <fstream>
#include "config.h"

#include "scene.h"
#include "CModelBase.h"

#include "physicsbase.h"
#include "physicsthomson.h"
// #include "physicsuv.h"
#include "physicsvsf.h"
#include "physicsvsfvarydist.h"
#include "rtthread.h"
#include "readDataScatteringTest.h"


#define BOOST_TEST_MODULE ScatteringTest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

const char* fname1 = SCRAYTRACE_DATA_DIR "/ScatteringTest/Test1.dat";
const char* fname2 = SCRAYTRACE_DATA_DIR "/ScatteringTest/Test2.dat";


struct ScatteringTest {
   PhysicsBase base;
   PhysicsBase *pthom; //,*puv;


ScatteringTest() { 
//      puv = physicsSelect(UV);

}

  ~ScatteringTest() { 
//      delete puv;
}

};

                  
BOOST_FIXTURE_TEST_SUITE(s, ScatteringTest)

  BOOST_AUTO_TEST_CASE(test_ScatteringTest)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running ScatteringTest");

   int sx = 256;
   int sy = 256;

   float* realData = readDatFile(fname1);

   float fovpix = 0.000681769;
   float obspos[] = {1.2807016e-01, -2.1597208e+02, 1.6441868e+01};
   float obsang[] = {0., 0., -0.05235988};
   float nepos[] = {0., 0., 0.};
   float neang[] = {0., 0., 0.};
   int losnbp = 2000;
   float losrange[] = {115., 265.};
   int modelid = 73;
   float pmodparam[] = {0.00118675, 0., 0.};
   float crpix[] = {127.5, 127.5};
   int quiet = 1;
   int neonly = 0;
   float hlonlat[] = {0., 0., 0.};
   float occrad = 0.0;
   float limbdark = 0.0;
   float obslonlat[] = {4.6364059e00, 5.9128302e-04, 2.1659708e02};
   int obslonlatflag = 1;
   unsigned int projtypecode = 1;
   float pv2_1 = 0.0;
   float pc[] = {1., 0., 0., 1.};
   int frontinteg = 1;
   unsigned int nbthreads = 32;
   unsigned int nbchunks = 32;
   float nerotcntr[] = {0., 0., 0.};
   float nerotang[] = {0., 0., 0.};
   float netranslation[] = {0., 0., 0.};
   int nerotaxis[] = {3, 2, 1};
   int physics = 5;
   float phyparam[] = {0.58};
   float fracmax = 0.0;
   int runDumpInteg = 0;
   float pIntegrand[] = {0.};

   float* btot = (float*)calloc(sy*sx, sizeof(float*));
   float* bpol = (float*)calloc(sy*sx, sizeof(float*));
   float* netot = (float*)calloc(sy*sx, sizeof(float*));

    rtthread(sx,
            sy,
            fovpix,
            obspos,
            obsang,
            nepos,
            neang,
            losnbp,
            losrange,
            modelid,
            btot,
            bpol,
            netot,
            pmodparam,
            crpix,
            quiet,
            neonly,
            hlonlat,
            occrad,
            limbdark,
            obslonlat,
            obslonlatflag,
            projtypecode,
            pv2_1,
            pc,
            frontinteg,
            nbthreads,
            nbchunks,
            nerotcntr,
            nerotang,
            netranslation,
            nerotaxis,
            physics,
            phyparam,
            fracmax,
            runDumpInteg, 
            pIntegrand);
 
    for(int i = 0; i < sy*sx; i++){
        std::stringstream temp;
        temp << std::scientific << abs(realData[i] - btot[i]);
        std::string prctErr = temp.str();
        float diffMag = 1*pow(10, std::stof(prctErr.substr(prctErr.find('e') + 1, prctErr.length() - 1)) + 1);
        if(abs(realData[i] - btot[i]) > diffMag){
            BOOST_TEST(false);
        }
    }

    free(realData);
    free(btot);
    free(bpol);
    free(netot);

    BOOST_TEST(true);
     
}

BOOST_AUTO_TEST_SUITE_END()

