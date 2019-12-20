
#include "rtmiscfunc.h"
#include "Cvec.h"
#include "Cbasis.h"
#include "constant.h"


#define BOOST_TEST_MODULE rtmiscfunctest
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>


#define NBSAMP 3

struct rtmiscfunctest {

rtmiscfunctest() { 

}

  ~rtmiscfunctest() { 
     
}
  static const float xarray[NBSAMP];
  static const float yarray[NBSAMP];
  static const unsigned int nbsamp;

};
const float rtmiscfunctest::xarray[NBSAMP]={1,2,3};
const float rtmiscfunctest::yarray[NBSAMP]={-2,-1,0};
const unsigned int rtmiscfunctest::nbsamp=NBSAMP;

                  


BOOST_FIXTURE_TEST_SUITE(s, rtmiscfunctest)

  BOOST_AUTO_TEST_CASE(test_Cvec)
  {
      
   namespace tt = boost::test_tools;
   
   BOOST_TEST_MESSAGE("running rtmiscfunctest");
     
   
  
short flagout;
float idshift;
unsigned int id;
float xin;
// ---- test too low
xin=0.5;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
BOOST_TEST(-1 == flagout);

// ---- test too high
xin=3.5;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
BOOST_TEST(1 == flagout);

// ---- test limit low
xin=1.;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
BOOST_TEST(0 == id);
BOOST_TEST(0 == flagout);
BOOST_TEST(0. == idshift, tt::tolerance(0.0001));

// ---- test limit high
xin=3.;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
BOOST_TEST(2 == id);
BOOST_TEST(0 == flagout);
BOOST_TEST(0. == idshift, tt::tolerance(0.0001));

// ---- test within limits
xin=1.3;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
BOOST_TEST(0 == id);
BOOST_TEST(0 == flagout);
BOOST_TEST(0.3 == idshift, tt::tolerance(0.0001));
xin=2.7;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
BOOST_TEST(1 == id);
BOOST_TEST(0 == flagout);
BOOST_TEST(0.7 == idshift, tt::tolerance(0.0001));
xin=2.;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
BOOST_TEST(1 == id);
BOOST_TEST(0 == flagout);
BOOST_TEST(0. == idshift, tt::tolerance(0.0001));

// ---- test with an array of negative numbers
xin=-2.;
id=searchnearestindex02(xin,yarray,nbsamp,flagout,idshift);
BOOST_TEST(0 == id);
BOOST_TEST(0 == flagout);
BOOST_TEST(0. == idshift, tt::tolerance(0.0001));
xin=-1.2;
id=searchnearestindex02(xin,yarray,nbsamp,flagout,idshift);
BOOST_TEST(0 == id);
BOOST_TEST(0 == flagout);
BOOST_TEST(0.8 == idshift, tt::tolerance(0.0001));

// ---- test ChangetoDensityCoord
Cbasis nps(Cvec(0,0,0),0,0,PI/2),necoord2(Cvec(2,1,0),0,0,0),necoord3(Cvec(0,1,0),0,0,PI/2);
Cvec vs(1,0,0);
Cvec v=ChangetoDensityCoord(nps,vs);

BOOST_TEST( 0. == v.v[0], tt::tolerance(0.0001));
BOOST_TEST(-1. == v.v[1], tt::tolerance(0.0001));
BOOST_TEST( 0. == v.v[2], tt::tolerance(0.0001));


v=ChangetoDensityCoord(nps,vs,necoord2);

BOOST_TEST( 0. == v.v[0], tt::tolerance(0.0001));
BOOST_TEST(-1. == v.v[1], tt::tolerance(0.0001));
BOOST_TEST( 0. == v.v[2], tt::tolerance(0.0001));


v=ChangetoDensityCoord(nps,vs,necoord3);

BOOST_TEST(-2. == v.v[0], tt::tolerance(0.0001));
BOOST_TEST( 1. == v.v[1], tt::tolerance(0.0001));
BOOST_TEST( 0. == v.v[2], tt::tolerance(0.0001));
  
}

BOOST_AUTO_TEST_SUITE_END()

