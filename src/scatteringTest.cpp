
#include <iostream>
#include <string>
#include "Cvec.h"
#include <fstream>

#include "scene.h"
#include "CModelBase.h"

#include "physicsbase.h"
#include "physicsthomson.h"
#include "physicsuv.h"
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
   PhysicsBase *pthom,*puv;


ScatteringTest() { 
     puv = physicsSelect(UV);

}

  ~ScatteringTest() { 
     delete puv;
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

   float fovpix = 0.00136353847812057;
   float obspos[] = {0., 0., -205.99867};
   float obsang[] = {-0., 0., 0.05042338};
//    float* obsang = [-3.14159265, 0., 0.05042338];
   float nepos[] = {0., 0., 0.};
   float neang[] = {0., 0., 0.};
   int losnbp = 2000;
   float losrange[] = {115., 265.};
   int modelid = 83;
   float pmodparam[] = {1.0, 0., 0., 1.6891547036652756, 450.18827707689974, 174.174766817585, 5.};
   float crpix[] = {127.5, 127.5};
   int quiet = 1;
   int neonly = 0;
   float hlonlat[] = {0., 0., 0.};
   float occrad = 0.0;
   float limbdark = 0.0;
   float obslonlat[] = {0., 0., 215.};
   int obslonlatflag = 0;
   unsigned int projtypecode = 4;
   float pv2_1 = 0.0;
   float pc[] = {1., 0., 0., 1.};
   int frontinteg = 1;
   unsigned int nbthreads = 32;
   unsigned int nbchunks = 32;
   float nerotcntr[] = {0., 0., 0.};
   float nerotang[] = {0., 0., 0.};
   float netranslation[] = {0., 0., 0.};
   int nerotaxis[] = {3, 2, 1};
   int physics = 6;
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
        if(!(realData[i] - realData[i]*0.079 <= btot[i] && realData[i] + realData[i]*0.079 >= btot[i])){ //Could it be that the python code is saving the data as double precision float but C is reading as single and that's what's causing the issue? (DEBUG)
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

