
#include <iostream>

#include "constant.h"
#include "Cvec.h"
#include "Cbasis.h"
#include "models21to30.h"

#define BOOST_TEST_MODULE model26Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>



// #include "model26test.h"

#define SX 3
#define SY 3
#define SZ 3
#define CNTRX 1.5
#define CNTRY 1.5
#define CNTRZ 1.5
#define VOXSIZERSUN 1.



struct model26Test {

  CModel26 *a;

  float *pparam;

  static const unsigned int sx;
  static const unsigned int sy;
  static const unsigned int sz;
  static const float cntrx;
  static const float cntry;
  static const float cntrz;
  static const float voxsizersun;

  static const float pcube[SX*SY*SZ];



model26Test() { 
 
    a = new CModel26;
 
    pparam=new float[7*sx*sy*sz];

    pparam[0]=sx;
    pparam[1]=sy;
    pparam[2]=sz;
    pparam[3]=cntrx;
    pparam[4]=cntry;
    pparam[5]=cntrz;
    pparam[6]=voxsizersun;

    for(unsigned int i=0;i<sx*sy*sz;i++) pparam[7+i]=pcube[i];

  a->initParam(pparam);

}

  ~model26Test() { 
    
    delete a;
    delete[] pparam;
}
 


};

const unsigned int model26Test::sx=SX;
const unsigned int model26Test::sy=SY;
const unsigned int model26Test::sz=SZ;
const float model26Test::cntrx=CNTRX;
const float model26Test::cntry=CNTRY;
const float model26Test::cntrz=CNTRZ;
const float model26Test::voxsizersun=VOXSIZERSUN;


const float model26Test::pcube[SX*SY*SZ]={   0,0,0 , 0,0,0 , 0,0,0,
                                             1,1,1 , 1,3,1 , 1,1,10,
                                             2,2,2 , 2,2,2 , 2,2,20};


                  
BOOST_FIXTURE_TEST_SUITE(s, model26Test)

  BOOST_AUTO_TEST_CASE(test_model26)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running model26Test");


    float dens=-1.;
    Cvec pos(0., 0., 0.);
    dens=a->Density(pos);
    BOOST_TEST(3. == dens, tt::tolerance(0.0001));

    pos=Cvec(0.,0.,1.);
    dens=a->Density(pos);
    BOOST_TEST(2. == dens, tt::tolerance(0.0001));

    pos=Cvec(0.,0.,-1.5);
    dens=a->Density(pos);
    BOOST_TEST(0. == dens, tt::tolerance(0.0001));

    pos=Cvec(1.49,1.49,1.49);
    dens=a->Density(pos);
    BOOST_TEST(20. == dens, tt::tolerance(0.0001));


    pos=Cvec(0.,0.,.5);
    dens=a->Density(pos);
    BOOST_TEST(2.5 == dens, tt::tolerance(0.0001));

}

BOOST_AUTO_TEST_SUITE_END()
