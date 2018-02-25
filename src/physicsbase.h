// $Id: physicsbase.h,v 1.1 2010-09-01 19:54:30 thernis Exp $

#ifndef PHYSICSBASE_H
#define PHYSICSBASE_H

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include "Cvec.h"
#include "CModelBase.h"
             
//! Enumerate the different types of physics available
enum PhysicsType {THOMSON,UV,MIE,ISOTROPIC,VSF};

//! Forward declaration of class Scene
class Scene;
         
/**
Virtual class implementing the physics of the raytracing.
*/
class PhysicsBase{
  public:
    PhysicsBase() {physicsName="Physics Base";};

    virtual ~PhysicsBase() {};

    //! Compute the radiation for the given physics, geometry, model
    virtual bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout) {return 1;}

    virtual void getConstFactors(float &btf,float &bpf,float &nef)
    {
        btf=1.;
        bpf=1.;
        nef=1.;
    };
    string getPhysics(){return physicsName;}
    virtual void setParam(float *phyparam) {};
    virtual void printParam() {std::cout << "No PhyParam." << std::endl;};

    void setParentScene(Scene *pparentscene) {this->pparentscene=pparentscene;}
        
  protected:
    string physicsName;
    Scene *pparentscene;
};


//! physics selection function
PhysicsBase* physicsSelect(PhysicsType phytype);


#endif
