
#include "constant.h"
#include "Cvec.h"

#define BOOST_TEST_MODULE model67Test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

// -- The define private public is a hack to switch all the private variables of the model61to70.h to public. This is to allow unit testing of private members. Apparently, this will probably fail with a Windows compiler.
#define private public
#include "models61to70.h"


#define SR 3
#define SLON 4
#define SLAT 5
#define SC 2


struct model67Test {

model67Test() { 
    // ---- create a model instance
    a=new CModel67;
  b=new CModel67;

  // ---- create a density and temperature cube
  // -- allocate the cube
  cube=new float[CubeTotalSize];
  cubeb=new float[CubeTotalSize];

  // -- fill in the cube
  cube[0]=(float) sr;
  cube[1]=(float) slon;
  cube[2]=(float) slat;
  cubeb[0]=(float) sr;
  cubeb[1]=(float) slon;
  cubeb[2]=(float) slat;
  for (unsigned int i=0;i<sr;i++) cube[3+i]=r[i];
  for (unsigned int i=0;i<slon;i++) cube[3+sr+i]=lon[i];
  for (unsigned int i=0;i<slat;i++) cube[3+sr+slon+i]=lat[i];
  for (unsigned int i=0;i<sr;i++) cubeb[3+i]=r[i];
  for (unsigned int i=0;i<slon;i++) cubeb[3+sr+i]=lon[i];
  for (unsigned int i=0;i<slat;i++) cubeb[3+sr+slon+i]=lat[i];

  for (unsigned int k=0;k<(slat);k++) for (unsigned int j=0;j<(slon);j++) for (unsigned int i=0;i<(sr);i++)
{
    cube[offsetnele+i+j*sr+k*sr*slon]=r[i];
    cube[offsettemp+i+j*sr+k*sr*slon]=lon[j];
}

  for (unsigned int k=0;k<(slat);k++) for (unsigned int j=0;j<(slon);j++) for (unsigned int i=0;i<(sr);i++)
{
    cubeb[offsetnele+i+j*sr+k*sr*slon]=r[i];
    cubeb[offsettemp+i+j*sr+k*sr*slon]=lat[k];
}

  a->initParam(cube);
  b->initParam(cubeb);

// ---- simple cube to test trilininterp function
cucube=new float[2*2*2];
for (unsigned int i=0;i<8;i++) cucube[i]=0;
cucube[0]=1;

}

  ~model67Test() { 
    delete a,
    delete b;
    delete[] cube;
    delete[] cubeb;
    delete[] cucube;
}


  CModel67 *a,*b;
  float *cube;
  float *cubeb;
    float *cucube;

  static const unsigned int sr;
  static const unsigned int slon;
  static const unsigned int slat;

  static const float r[SR];
  static const float lon[SLON];
  static const float lat[SLAT];

  static const unsigned int CubeTotalSize;

  static const unsigned int offsetnele;
  static const unsigned int offsettemp;

    static const unsigned int sc;




};

const unsigned int model67Test::sr=SR;
const unsigned int model67Test::slon=SLON;
const unsigned int model67Test::slat=SLAT;

const float model67Test::r[SR]={1,2,3};
const float model67Test::lon[SLON]={0,PI/2,PI,3*PI/2};
const float model67Test::lat[SLAT]={-PI/2,-PI/4,0.,PI/4,PI/2};

const unsigned int model67Test::CubeTotalSize=3+SR+SLON+SLAT+SR*SLON*SLAT*2;

const unsigned int model67Test::offsetnele=3+SR+SLON+SLAT;
const unsigned int model67Test::offsettemp=3+SR+SLON+SLAT+SR*SLON*SLAT;

const unsigned int sc=SC;




                  
BOOST_FIXTURE_TEST_SUITE(s, model67Test)

  BOOST_AUTO_TEST_CASE(test_model67)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running model67Test");
   
   // ---- test trilininterp function
    float neout;
    neout=trilininterp(0.,0.,0.,0,0,0,cucube,2,2);
    BOOST_TEST(1. == neout, tt::tolerance(0.0001));

    neout=trilininterp(1.,0.,0.,0,0,0,cucube,2,2);
    BOOST_TEST(0. == neout, tt::tolerance(0.0001));

    neout=trilininterp(0.5,0.,0.,0,0,0,cucube,2,2);
    BOOST_TEST(0.5 == neout, tt::tolerance(0.0001));

    neout=trilininterp(0.,0.6,0.,0,0,0,cucube,2,2);
    BOOST_TEST(0.4 == neout, tt::tolerance(0.0001));

    neout=trilininterp(0.,0.,0.3,0,0,0,cucube,2,2);
    BOOST_TEST(0.7 == neout, tt::tolerance(0.0001));


    cucube[0]=0;
    cucube[7]=1;

    neout=trilininterp(0.,0.,0.,0,0,0,cucube,2,2);
    BOOST_TEST(0. == neout, tt::tolerance(0.0001));

    neout=trilininterp(1.,1.,1.,0,0,0,cucube,2,2);
    BOOST_TEST(1. == neout, tt::tolerance(0.0001));



    // ---- test initialization
    // -- object a
    BOOST_TEST(sr == a->sr);
    BOOST_TEST(slon == a->slon);
    BOOST_TEST(slat == a->slat);

    BOOST_TEST(r[0] == a->r[0], tt::tolerance(0.0001));
    BOOST_TEST(lon[0] == a->lon[0], tt::tolerance(0.0001));
    BOOST_TEST(lat[0] == a->lat[0], tt::tolerance(0.0001));

    BOOST_TEST(cube[offsetnele] == a->nele[0], tt::tolerance(0.0001));
    BOOST_TEST(cube[offsettemp] == a->temp[0], tt::tolerance(0.0001));

    BOOST_TEST(cube[offsetnele+2] == a->nele[2], tt::tolerance(0.0001));
    BOOST_TEST(cube[offsettemp+2] == a->temp[2], tt::tolerance(0.0001));

    // -- object b
    BOOST_TEST(cubeb[offsettemp] == b->temp[0], tt::tolerance(0.0001));
    BOOST_TEST(cubeb[offsettemp+2] == b->temp[2], tt::tolerance(0.0001));

    BOOST_TEST(cubeb[offsettemp+1+1*sr+0*sr*slon] == lat[0], tt::tolerance(0.0001));
    BOOST_TEST(cubeb[offsettemp+1+1*sr+1*sr*slon] == lat[1], tt::tolerance(0.0001));
    BOOST_TEST(cubeb[offsettemp+1+1*sr+2*sr*slon] == lat[2], tt::tolerance(0.0001));

    // ---- test radius
    float temp=-1.;
    float dens=-1.;
    Cvec pos(0.,0.,1.);
    dens=a->Density(pos,temp);
    BOOST_TEST(r[0] == dens, tt::tolerance(0.0001));
    pos=Cvec(0.,0.,2.);
    dens=a->Density(pos,temp);
    BOOST_TEST(r[1] == dens, tt::tolerance(0.0001));
    pos=Cvec(0.,0.,3.);
    dens=a->Density(pos,temp);
    BOOST_TEST(r[2] == dens, tt::tolerance(0.0001));

    // ---- test longitude
    pos=Cvec(0.,1.,0.);
    dens=a->Density(pos,temp);
    BOOST_TEST((3 * PI / 2) == temp, tt::tolerance(0.0001));

    // ---- test latitude
    pos=Cvec(0.0,0.,1.);
    dens=b->Density(pos,temp);
    BOOST_TEST(0. == temp, tt::tolerance(0.0001));

    pos=Cvec(2./sqrt(2),0.,2./sqrt(2));
    dens=b->Density(pos,temp);
    BOOST_TEST((PI / 4) == temp, tt::tolerance(0.0001));

    pos=Cvec(-1.00001/sqrt(2),0.,1.00001/sqrt(2));
    dens=b->Density(pos,temp);
    BOOST_TEST((-PI / 4) == temp, tt::tolerance(0.0001));

    pos=Cvec(0,2.8,0.);
    dens=b->Density(pos,temp);
    cout << "dens : " << dens << endl;
    cout << "temp : " << temp << endl;

}

BOOST_AUTO_TEST_SUITE_END()

