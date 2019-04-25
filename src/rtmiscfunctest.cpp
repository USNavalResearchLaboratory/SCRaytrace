
#include "rtmiscfunctest.h"
#include "rtmiscfunc.h"
#include "Cvec.h"
#include "Cbasis.h"
#include "constant.h"

CPPUNIT_TEST_SUITE_REGISTRATION (rtmiscfunctest);


void rtmiscfunctest :: setUp (void)
{

}

void rtmiscfunctest :: tearDown (void)
{

}

void rtmiscfunctest :: testrtmiscfunc (void)
{
short flagout;
float idshift;
unsigned int id;
float xin;
// ---- test too low
xin=0.5;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL(-1,(int)flagout);

// ---- test too high
xin=3.5;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL(1,(int)flagout);

// ---- test limit low
xin=1.;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL(0,(int)id);
CPPUNIT_ASSERT_EQUAL(0,(int)flagout);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0,(double)idshift,1e-4);

// ---- test limit high
xin=3.;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL(2,(int)id);
CPPUNIT_ASSERT_EQUAL(0,(int)flagout);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,(double)idshift,1e-4);

// ---- test within limits
xin=1.3;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL((unsigned int) 0,id);
CPPUNIT_ASSERT_EQUAL(0,(int)flagout);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.3,(double)idshift,1e-4);
xin=2.7;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL(1,(int)id);
CPPUNIT_ASSERT_EQUAL(0,(int)flagout);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7,(double)idshift,1e-4);
xin=2.;
id=searchnearestindex02(xin,xarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL(1,(int)id);
CPPUNIT_ASSERT_EQUAL(0,(int)flagout);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,(double)idshift,1e-4);

// ---- test with an array of negative numbers
xin=-2.;
id=searchnearestindex02(xin,yarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL((unsigned int) 0,id);
CPPUNIT_ASSERT_EQUAL(0,(int)flagout);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,(double)idshift,1e-4);
xin=-1.2;
id=searchnearestindex02(xin,yarray,nbsamp,flagout,idshift);
CPPUNIT_ASSERT_EQUAL((unsigned int) 0,id);
CPPUNIT_ASSERT_EQUAL(0,(int)flagout);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.8,(double)idshift,1e-4);

// ---- test ChangetoDensityCoord
Cbasis nps(Cvec(0,0,0),0,0,PI/2),necoord2(Cvec(2,1,0),0,0,0),necoord3(Cvec(0,1,0),0,0,PI/2);
Cvec vs(1,0,0);
Cvec v=ChangetoDensityCoord(nps,vs);

CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.,(double)v.v[0],1e-4);
CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.,(double)v.v[1],1e-4);
CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.,(double)v.v[2],1e-4);


v=ChangetoDensityCoord(nps,vs,necoord2);

CPPUNIT_ASSERT_DOUBLES_EQUAL( 0,(double)v.v[0],1e-4);
CPPUNIT_ASSERT_DOUBLES_EQUAL(-1,(double)v.v[1],1e-4);
CPPUNIT_ASSERT_DOUBLES_EQUAL( 0,(double)v.v[2],1e-4);


v=ChangetoDensityCoord(nps,vs,necoord3);

CPPUNIT_ASSERT_DOUBLES_EQUAL(-2,(double)v.v[0],1e-4);
CPPUNIT_ASSERT_DOUBLES_EQUAL( 1,(double)v.v[1],1e-4);
CPPUNIT_ASSERT_DOUBLES_EQUAL( 0,(double)v.v[2],1e-4);





}
