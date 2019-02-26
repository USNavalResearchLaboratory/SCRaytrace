// $Id: physicstest.cpp,v 1.1 2010-09-01 19:56:55 thernis Exp $

#include "physicstest.h"
#include <iostream>
#include <string>
#include "Cvec.h"

#include "scene.h"
#include "CModelBase.h"

CPPUNIT_TEST_SUITE_REGISTRATION (PhysicsTest);

void PhysicsTest :: setUp (void)
{
    // ---- set up test environment (initializing objects)
    puv = physicsSelect(UV);

}

void PhysicsTest :: tearDown (void)
{
    // ---- cleaning destructor
    delete puv;
}

void PhysicsTest :: generalTests (void)
{
    string s1;
    bool flagok;
    float bto,bpo,neo;
    Scene *pscene;
    pscene=new Scene;
    
    base.setParentScene(pscene);
    
    unsigned int sx=256,sy=256;
    pscene->camera.setCCD(CCD(sx,sy));
    pscene->camera.setProjType(ARC);
    pscene->camera.setFovpix(0.01);
    pscene->camera.setCrpix(63.5,63.5);

    // -- LOS integration parameter definition
    
    int losnbp=20;
    float losrange[2]={190,230};
    pscene->los.setLOS(losnbp,losrange[0],losrange[1]);

    // -- density model
    
    int modelid=14;
    float *pmodparam;
    pmodparam=NULL;
    pscene->setDensityModel(14,pmodparam);
    pscene->setobs(Cbasis());
    
    s1=base.getPhysics();
    CPPUNIT_ASSERT (string("Physics Base").compare(s1)==0);
    
    flagok=base.computeRadiation(Cvec(0,0,0),1,1,bto,bpo,neo);
    CPPUNIT_ASSERT_EQUAL ((bool)1,flagok);

/*    pscene->setPhysics(FCOR);
    s1=pscene->getPhysics();
    CPPUNIT_ASSERT (string("F Corona").compare(s1)==0);*/
    
    puv->printParam();
    float phyparam=1;
    puv->setParam(&phyparam);
    puv->printParam();
    
    delete pscene;
    
    // -- test physicsVSF constructor
    PhysicsVSF vsf;

    // -- test physicsVSFVaryDist constructor
    PhysicsVSFVaryDist vsfvd;
    
    
//    pvsf = physicsSelect(VSF);
//      delete pvsf;

    
    
}

