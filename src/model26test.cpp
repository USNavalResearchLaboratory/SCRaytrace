// $Id: model26test.cpp,v 1.1 2009/03/17 14:45:12 thernis Exp $


#include "constant.h"
#include "Cvec.h"
#include "model26test.h"

CPPUNIT_TEST_SUITE_REGISTRATION (model26Test);


void model26Test :: setUp (void)
{
	// ---- create a model instance
	a=new CModel26;
 
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

void model26Test :: tearDown (void)
{
	// finally delete objects
	delete a;
  delete[] pparam;
}

void model26Test :: testmodel26 (void)
{


float dens=-1.;
Cvec pos(0.,0.,0.);
dens=a->Density(pos);
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)3.,(double)dens,1e-3);


pos=Cvec(0.,0.,1.);
dens=a->Density(pos);
cout << "dens : " << dens << endl;
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)2.,(double)dens,1e-3);

pos=Cvec(0.,0.,-1.5);
dens=a->Density(pos);
cout << "dens : " << dens << endl;
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)0.,(double)dens,1e-3);

pos=Cvec(1.49,1.49,1.49);
dens=a->Density(pos);
cout << "dens : " << dens << endl;
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)20.,(double)dens,1e-3);


pos=Cvec(0.,0.,.5);
dens=a->Density(pos);
cout << "dens : " << dens << endl;
CPPUNIT_ASSERT_DOUBLES_EQUAL((double)2.5,(double)dens,1e-3);




// ---- test trilininterp function
// float neout;
// neout=trilininterp(0.,0.,0.,0,0,0,cucube,2,2);
// CPPUNIT_ASSERT_DOUBLES_EQUAL(1,(double)neout,1e-3);

// ---- test initialization
// -- object a
//CPPUNIT_ASSERT_EQUAL(sr,a->sr);
//CPPUNIT_ASSERT_EQUAL(slon,a->slon);
//CPPUNIT_ASSERT_EQUAL(slat,a->slat);

}


// $Log: model26test.cpp,v $
// Revision 1.1  2009/03/17 14:45:12  thernis
// First commit
//
