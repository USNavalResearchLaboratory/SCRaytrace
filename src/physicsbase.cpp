/** \file physicsbase.cpp
 * \brief Selection of the requested physics. Edit this file if a new physics is implemented.
 */

#include "physicsbase.h"
#include "physicsthomson.h"
#include "physicsuv.h"
#include "physicsmie.h"
#include "physicsisotropic.h"
#include "physicsvsf.h"
#include "physicsvsfvarydist.h"
#include "physicsvariablevsf.h"


PhysicsBase* physicsSelect(PhysicsType phytype)
{
  PhysicsBase *pphy;
  PhysicsThomson *ppthom;
  switch (phytype)
  {
    case THOMSON :
      ppthom = new PhysicsThomson;
      pphy = (PhysicsBase*) ppthom;
      break;
    case UV :
      PhysicsUV *ppuv;
      ppuv = new PhysicsUV;
      pphy = (PhysicsBase*) ppuv;
      break;
    case MIE :
      PhysicsMie *pmie;
      pmie = new PhysicsMie;
      pphy = (PhysicsBase*) pmie;
      break;
    case ISOTROPIC :
      PhysicsIsotropic *piso;
      piso = new PhysicsIsotropic;
      pphy = (PhysicsBase*) piso;
      break;
    case VSF :
      PhysicsVSF *pvsf;
      pvsf = new PhysicsVSF;
      pphy = (PhysicsBase*) pvsf;
      break;
    case VSFVARYDIST :
      PhysicsVSFVaryDist *pvsfvd;
      pvsfvd = new PhysicsVSFVaryDist;
      pphy = (PhysicsBase*) pvsfvd;
      break;
    case VARIVSF :
      PhysicsVariableVSF *pvvsf;
      pvvsf = new PhysicsVariableVSF;
      pphy = (PhysicsBase*) pvvsf;
      break;
   default :
      ppthom = new PhysicsThomson;
      pphy = (PhysicsBase*) ppthom;
  }
  return pphy;
}
