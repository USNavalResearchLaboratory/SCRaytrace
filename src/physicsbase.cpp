// $Id: physicsbase.cpp,v 1.1 2010-09-01 19:56:55 thernis Exp $

#include "physicsbase.h"
#include "physicsthomson.h"
#include "physicsuv.h"
#include "physicsmie.h"
#include "physicsisotropic.h"
#include "physicsvsf.h"
#include "physicsvsfvarydist.h"


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
   default :
      ppthom = new PhysicsThomson;
      pphy = (PhysicsBase*) ppthom;
  }
  return pphy;
}
