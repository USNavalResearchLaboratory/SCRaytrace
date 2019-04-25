
#include "constant.h"
#include "Cvec.h"
#include "model67test.h"

CPPUNIT_TEST_SUITE_REGISTRATION (model67Test);


void model67Test :: setUp (void)
{
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

void model67Test :: tearDown (void)
{
	// finally delete objects
	delete a,b;
  delete[] cube,cubeb,cucube;
}

void model67Test :: testmodel67 (void)
{

// ---- test trilininterp function
float neout;
neout=trilininterp(0.,0.,0.,0,0,0,cucube,2,2);
CPPUNIT_ASSERT_DOUBLES_EQUAL(1,(double)neout,1e-3);

neout=trilininterp(1.,0.,0.,0,0,0,cucube,2,2);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0,(double)neout,1e-3);

neout=trilininterp(0.5,0.,0.,0,0,0,cucube,2,2);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,(double)neout,1e-3);

neout=trilininterp(0.,0.6,0.,0,0,0,cucube,2,2);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4,(double)neout,1e-3);

neout=trilininterp(0.,0.,0.3,0,0,0,cucube,2,2);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7,(double)neout,1e-3);


cucube[0]=0;
cucube[7]=1;

neout=trilininterp(0.,0.,0.,0,0,0,cucube,2,2);
CPPUNIT_ASSERT_DOUBLES_EQUAL(0,(double)neout,1e-3);

neout=trilininterp(1.,1.,1.,0,0,0,cucube,2,2);
CPPUNIT_ASSERT_DOUBLES_EQUAL(1,(double)neout,1e-3);



// ---- test initialization
// -- object a
CPPUNIT_ASSERT_EQUAL(sr,a->sr);
CPPUNIT_ASSERT_EQUAL(slon,a->slon);
CPPUNIT_ASSERT_EQUAL(slat,a->slat);

CPPUNIT_ASSERT_DOUBLES_EQUAL((double)r[0],(double)a->r[0],1e-3);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)lon[0],(double)a->lon[0],1e-3);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)lat[0],(double)a->lat[0],1e-3);

CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cube[offsetnele],(double)a->nele[0],1e-3);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cube[offsettemp],(double)a->temp[0],1e-3);

CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cube[offsetnele+2],(double)a->nele[2],1e-3);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cube[offsettemp+2],(double)a->temp[2],1e-3);

// -- object b
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cubeb[offsettemp],(double)b->temp[0],1e-3);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cubeb[offsettemp+2],(double)b->temp[2],1e-3);

CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cubeb[offsettemp+1+1*sr+0*sr*slon],(double)lat[0],1e-3);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cubeb[offsettemp+1+1*sr+1*sr*slon],(double)lat[1],1e-3);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)cubeb[offsettemp+1+1*sr+2*sr*slon],(double)lat[2],1e-3);




// ---- test radius
float temp=-1.;
float dens=-1.;
Cvec pos(0.,0.,1.);
dens=a->Density(pos,temp);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)r[0],(double)dens,1e-3);
pos=Cvec(0.,0.,2.);
dens=a->Density(pos,temp);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)r[1],(double)dens,1e-3);
pos=Cvec(0.,0.,3.);
dens=a->Density(pos,temp);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)r[2],(double)dens,1e-3);

// ---- test longitude
pos=Cvec(0.,1.,0.);
dens=a->Density(pos,temp);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)3*PI/2,(double)temp,1e-3);

// ---- test latitude
pos=Cvec(0.0,0.,1.);
dens=b->Density(pos,temp);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)0,(double)temp,1e-3);

pos=Cvec(2./sqrt(2),0.,2./sqrt(2));
dens=b->Density(pos,temp);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)PI/4,(double)temp,1e-3);

pos=Cvec(-1.00001/sqrt(2),0.,1.00001/sqrt(2));
dens=b->Density(pos,temp);
CPPUNIT_ASSERT_DOUBLES_EQUAL(-(double)PI/4,(double)temp,1e-3);


pos=Cvec(0,2.8,0.);
dens=b->Density(pos,temp);
cout << "dens : " << dens << endl;
cout << "temp : " << temp << endl;
//CPPUNIT_ASSERT_DOUBLES_EQUAL((double)PI/4,(double)temp,1e-3);


}

