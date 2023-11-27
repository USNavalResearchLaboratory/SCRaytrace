/** \file physicsbase.h
 * \brief Base class for all the physics types.
 */

#ifndef PHYSICSBASE_H
#define PHYSICSBASE_H

#include <string>
using std::string;
#include <iostream>
using std::cout;
using std::endl;
#include "Cvec.h"
#include "CModelBase.h"
             
//! \brief Enumerate the different types of physics available
enum PhysicsType {THOMSON, UV, MIE, ISOTROPIC, VSF, VSFVARYDIST, VARIVSF};

//! Forward declaration of class Scene
class Scene;
         
/** \brief Virtual class implementing the physics of the raytracing.
*/
class PhysicsBase{
  public:
    PhysicsBase() {physicsName="Physics Base";};

    virtual ~PhysicsBase() {};

    
    //! Compute the radiation for the given physics, geometry, model    
    virtual bool computeRadiation(const Cvec &vs,   /**< [in] Point on the LOS, in Abs coordinates. */
                                  const float &r,   /**< [in] Distance LOS position to Sun center. */
                                  const float &rho, /**< [in] Impact parameter. */
                                  float &btout,     /**< [out] Total brightness. */
                                  float &bpout,     /**< [out] Polarized brightness. */
                                  float &neout)     /**< [out] Density. */
    {return 1;}

    
    //! Return the integration constant factor
    virtual void getConstFactors(float &btf,float &bpf,float &nef, float rho)
    {
        btf=1.;
        bpf=1.;
        nef=1.;
    };

    
    /** Perform Model density initialization tasks. Some model need the position of the observer.
     * \param vlosabs Line of sight unit vector, in the abs coordinate system.
     */
//     virtual void initDensityModel(const Cvec &vlosabs) {};
    
    
    string getPhysics(){return physicsName;}
    virtual void setParam(float *phyparam) {};
    virtual void printParam() {std::cout << "No PhyParam." << std::endl;};

    void setParentScene(Scene *pparentscene) {this->pparentscene=pparentscene;};
        
  protected:
    string physicsName;
    Scene *pparentscene;
};


//! physics selection function
PhysicsBase* physicsSelect(PhysicsType phytype);


#endif
