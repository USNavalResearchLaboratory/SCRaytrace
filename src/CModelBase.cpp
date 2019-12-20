/*! \file CModelBase.cpp 
 * \brief Implementation of the base class for all the density models.
 *  
 *  Add entry to this file when a new model is created.
 * 
 */



#include "CModelBase.h"

#include <iostream>
#include "models01to10.h"
#include "models11to20.h"
#include "models21to30.h"
#include "models31to40.h"
#include "models41to50.h"
#include "models51to60.h"
#include "models61to70.h"
#include "models71to80.h"
#include "CModel46.h"
#include "CModel48.h"


void CModelBase::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
	flagcase=1; // default parameter undefined
	return;
}


//! Model selection function
CModelBase* modelselect(int modelid)
{
	CModelBase *pmod;
 	CModel01 *pmodel1;
   // -- pick Ne density model func. pointer
  switch (modelid) {
  case 1 : 
	pmodel1 = new CModel01;
	pmod= (CModelBase*) pmodel1;
    return pmod;
    break;
   case 2 : 
	CModel02 *pmodel2;
	pmodel2 = new CModel02;
	pmod= (CModelBase*) pmodel2;
    return pmod;
    break;
   case 3 : 
	CModel03 *pmodel3;
	pmodel3 = new CModel03;
	pmod= (CModelBase*) pmodel3;
    return pmod;
    break;
   case 4 :
  	CModel04 *pmodel4;
	pmodel4 = new CModel04;
	pmod= (CModelBase*) pmodel4;
    return pmod;
    break;
  case 5 :
  	CModel05 *pmodel5;
	pmodel5 = new CModel05;
	pmod= (CModelBase*) pmodel5;
    return pmod;
    break;
  case 6 :
  	CModel06 *pmodel6;
	pmodel6 = new CModel06;
	pmod= (CModelBase*) pmodel6;
    return pmod;
    break;
  case 7 :
  	CModel07 *pmodel7;
	pmodel7 = new CModel07;
	pmod= (CModelBase*) pmodel7;
    return pmod;
    break;
  case 8 :
  	CModel08 *pmodel8;
	pmodel8 = new CModel08;
	pmod= (CModelBase*) pmodel8;
    return pmod;
    break;
  case 9 :
  	CModel09 *pmodel9;
	pmodel9 = new CModel09;
	pmod= (CModelBase*) pmodel9;
    return pmod;
    break;
  case 10 :
  	CModel10 *pmodel10;
	pmodel10 = new CModel10;
	pmod= (CModelBase*) pmodel10;
    return pmod;
    break;
  case 11 :
  	CModel11 *pmodel11;
	pmodel11 = new CModel11;
	pmod= (CModelBase*) pmodel11;
    return pmod;
    break;
  case 12 :
  	CModel12 *pmodel12;
	pmodel12 = new CModel12;
	pmod= (CModelBase*) pmodel12;
    return pmod;
    break;
  case 13 :
  	CModel13 *pmodel13;
	pmodel13 = new CModel13;
	pmod= (CModelBase*) pmodel13;
    return pmod;
    break;
  case 14 :
  	CModel14 *pmodel14;
	pmodel14 = new CModel14;
	pmod= (CModelBase*) pmodel14;
    return pmod;
    break;
  case 15 :
  	CModel15 *pmodel15;
	pmodel15 = new CModel15;
	pmod= (CModelBase*) pmodel15;
    return pmod;
    break;
  case 16 :
  	CModel16 *pmodel16;
	pmodel16 = new CModel16;
	pmod= (CModelBase*) pmodel16;
    return pmod;
    break;
  case 17 :
  	CModel17 *pmodel17;
	pmodel17 = new CModel17;
	pmod= (CModelBase*) pmodel17;
    return pmod;
    break;
  case 18 :
  	CModel18 *pmodel18;
	pmodel18 = new CModel18;
	pmod= (CModelBase*) pmodel18;
    return pmod;
    break;
  case 19 :
  	CModel19 *pmodel19;
	pmodel19 = new CModel19;
	pmod= (CModelBase*) pmodel19;
    return pmod;
    break;
  case 20 :
  	CModel20 *pmodel20;
	pmodel20 = new CModel20;
	pmod= (CModelBase*) pmodel20;
    return pmod;
    break;
  case 21 :
  	CModel21 *pmodel21;
	pmodel21 = new CModel21;
	pmod= (CModelBase*) pmodel21;
    return pmod;
    break;
  case 22 :
  	CModel22 *pmodel22;
	pmodel22 = new CModel22;
	pmod= (CModelBase*) pmodel22;
    return pmod;
    break;
  case 23 :
  	CModel23 *pmodel23;
	pmodel23 = new CModel23;
	pmod= (CModelBase*) pmodel23;
    return pmod;
    break;
  case 24 :
  	CModel24 *pmodel24;
	pmodel24 = new CModel24;
	pmod= (CModelBase*) pmodel24;
    return pmod;
    break;
  case 25 :
  	CModel25 *pmodel25;
	pmodel25 = new CModel25;
	pmod= (CModelBase*) pmodel25;
    return pmod;
    break;
  case 26 :
  	CModel26 *pmodel26;
	pmodel26 = new CModel26;
	pmod= (CModelBase*) pmodel26;
    return pmod;
    break;
  case 27 :
  	CModel27 *pmodel27;
	pmodel27 = new CModel27;
	pmod= (CModelBase*) pmodel27;
    return pmod;
    break;
  case 28 :
  	CModel28 *pmodel28;
	pmodel28 = new CModel28;
	pmod= (CModelBase*) pmodel28;
    return pmod;
    break;
  case 29 :
  	CModel29 *pmodel29;
	pmodel29 = new CModel29;
	pmod= (CModelBase*) pmodel29;
    return pmod;
    break;
  case 30 :
  	CModel30 *pmodel30;
	pmodel30 = new CModel30;
	pmod= (CModelBase*) pmodel30;
    return pmod;
    break;
  case 31 :
  	CModel31 *pmodel31;
	pmodel31 = new CModel31;
	pmod= (CModelBase*) pmodel31;
    return pmod;
    break;
  case 32 :
  	CModel32 *pmodel32;
	pmodel32 = new CModel32;
	pmod= (CModelBase*) pmodel32;
    return pmod;
    break;
  case 33 :
  	CModel33 *pmodel33;
	pmodel33 = new CModel33;
	pmod= (CModelBase*) pmodel33;
    return pmod;
    break;
  case 34 :
  	CModel34 *pmodel34;
	pmodel34 = new CModel34;
	pmod= (CModelBase*) pmodel34;
    return pmod;
    break;
  case 35 :
  	CModel35 *pmodel35;
	pmodel35 = new CModel35;
	pmod= (CModelBase*) pmodel35;
    return pmod;
    break;
  case 36 :
  	CModel36 *pmodel36;
	pmodel36 = new CModel36;
	pmod= (CModelBase*) pmodel36;
    return pmod;
    break;
  case 37 :
  	CModel37 *pmodel37;
	pmodel37 = new CModel37;
	pmod= (CModelBase*) pmodel37;
    return pmod;
    break;
  case 38 :
  	CModel38 *pmodel38;
	pmodel38 = new CModel38;
	pmod= (CModelBase*) pmodel38;
    return pmod;
    break;
  case 39 :
  	CModel39 *pmodel39;
	pmodel39 = new CModel39;
	pmod= (CModelBase*) pmodel39;
    return pmod;
    break;
  case 40 :
  	CModel40 *pmodel40;
	pmodel40 = new CModel40;
	pmod= (CModelBase*) pmodel40;
    return pmod;
    break;
  case 41 :
  	CModel41 *pmodel41;
	pmodel41 = new CModel41;
	pmod= (CModelBase*) pmodel41;
    return pmod;
    break;
  case 42 :
  	CModel42 *pmodel42;
	pmodel42 = new CModel42;
	pmod= (CModelBase*) pmodel42;
    return pmod;
    break;
  case 43 :
  	CModel43 *pmodel43;
	pmodel43 = new CModel43;
	pmod= (CModelBase*) pmodel43;
    return pmod;
    break;
  case 44 :
  	CModel44 *pmodel44;
	pmodel44 = new CModel44;
	pmod= (CModelBase*) pmodel44;
    return pmod;
    break;
  case 45 :
  	CModel45 *pmodel45;
	pmodel45 = new CModel45;
	pmod= (CModelBase*) pmodel45;
    return pmod;
    break;
  case 46 :
  	CModel46 *pmodel46;
	pmodel46 = new CModel46;
	pmod= (CModelBase*) pmodel46;
    return pmod;
    break;
  case 47 :
  	CModel47 *pmodel47;
	pmodel47 = new CModel47;
	pmod= (CModelBase*) pmodel47;
    return pmod;
    break;
  case 48 :
  	CModel48 *pmodel48;
	pmodel48 = new CModel48;
	pmod= (CModelBase*) pmodel48;
    return pmod;
    break;
 case 49 :
  	CModel49 *pmodel49;
	pmodel49 = new CModel49;
	pmod= (CModelBase*) pmodel49;
    return pmod;
    break;
 case 50 :
  	CModel50 *pmodel50;
	pmodel50 = new CModel50;
	pmod= (CModelBase*) pmodel50;
    return pmod;
    break;
 case 51 :
  	CModel51 *pmodel51;
	pmodel51 = new CModel51;
	pmod= (CModelBase*) pmodel51;
    return pmod;
    break;
 case 52 :
  	CModel52 *pmodel52;
	pmodel52 = new CModel52;
	pmod= (CModelBase*) pmodel52;
    return pmod;
    break;
 case 53 :
  	CModel53 *pmodel53;
	pmodel53 = new CModel53;
	pmod= (CModelBase*) pmodel53;
    return pmod;
    break;
 case 54 :
  	CModel54 *pmodel54;
	pmodel54 = new CModel54;
	pmod= (CModelBase*) pmodel54;
    return pmod;
    break;
 case 55 :
      CModel55 *pmodel55;
      pmodel55 = new CModel55;
      pmod= (CModelBase*) pmodel55;
      return pmod;
      break;
 case 56 :
      CModel56 *pmodel56;
      pmodel56 = new CModel56;
      pmod= (CModelBase*) pmodel56;
      return pmod;
      break;
    case 57 :
      CModel57 *pmodel57;
      pmodel57 = new CModel57;
      pmod= (CModelBase*) pmodel57;
      return pmod;
      break;
    case 58 :
      CModel58 *pmodel58;
      pmodel58 = new CModel58;
      pmod= (CModelBase*) pmodel58;
      return pmod;
      break;
    case 59 :
      CModel59 *pmodel59;
      pmodel59 = new CModel59;
      pmod= (CModelBase*) pmodel59;
      return pmod;
      break;
    case 60 :
      CModel60 *pmodel60;
      pmodel60 = new CModel60;
      pmod= (CModelBase*) pmodel60;
      return pmod;
      break;
    case 61 :
      CModel61 *pmodel61;
      pmodel61 = new CModel61;
      pmod= (CModelBase*) pmodel61;
      return pmod;
      break;
    case 62 :
      CModel62 *pmodel62;
      pmodel62 = new CModel62;
      pmod= (CModelBase*) pmodel62;
      return pmod;
      break;
    case 63 :
      CModel63 *pmodel63;
      pmodel63 = new CModel63;
      pmod= (CModelBase*) pmodel63;
      return pmod;
      break;
    case 64 :
      CModel64 *pmodel64;
      pmodel64 = new CModel64;
      pmod= (CModelBase*) pmodel64;
      return pmod;
      break;
    case 65 :
      CModel65 *pmodel65;
      pmodel65 = new CModel65;
      pmod= (CModelBase*) pmodel65;
      return pmod;
      break;
    case 66 :
      CModel66 *pmodel66;
      pmodel66 = new CModel66;
      pmod= (CModelBase*) pmodel66;
      return pmod;
      break;
    case 67 :
      CModel67 *pmodel67;
      pmodel67 = new CModel67;
      pmod= (CModelBase*) pmodel67;
      return pmod;
      break;
    case 68 :
      CModel68 *pmodel68;
      pmodel68 = new CModel68;
      pmod= (CModelBase*) pmodel68;
      return pmod;
      break;
    case 69 :
      CModel69 *pmodel69;
      pmodel69 = new CModel69;
      pmod= (CModelBase*) pmodel69;
      return pmod;
      break;
    case 70 :
      CModel70 *pmodel70;
      pmodel70 = new CModel70;
      pmod= (CModelBase*) pmodel70;
      return pmod;
      break;
    case 71 :
      CModel71 *pmodel71;
      pmodel71 = new CModel71;
      pmod= (CModelBase*) pmodel71;
      return pmod;
      break;
    case 72 :
      CModel72 *pmodel72;
      pmodel72 = new CModel72;
      pmod= (CModelBase*) pmodel72;
      return pmod;
      break;
    case 73 :
      CModel73 *pmodel73;
      pmodel73 = new CModel73;
      pmod= (CModelBase*) pmodel73;
      return pmod;
      break;
    case 74 :
      CModel74 *pmodel74;
      pmodel74 = new CModel74;
      pmod= (CModelBase*) pmodel74;
      return pmod;
      break;
    case 75 :
      CModel75 *pmodel75;
      pmodel75 = new CModel75;
      pmod= (CModelBase*) pmodel75;
      return pmod;
      break;
    case 76 :
      CModel76 *pmodel76;
      pmodel76 = new CModel76;
      pmod= (CModelBase*) pmodel76;
      return pmod;
      break;
    case 77 :
      CModel77 *pmodel77;
      pmodel77 = new CModel77;
      pmod= (CModelBase*) pmodel77;
      return pmod;
      break;
    case 78 :
      CModel78 *pmodel78;
      pmodel78 = new CModel78;
      pmod= (CModelBase*) pmodel78;
      return pmod;
      break;
    case 79 :
      CModel79 *pmodel79;
      pmodel79 = new CModel79;
      pmod= (CModelBase*) pmodel79;
      return pmod;
      break;
    case 80 :
      CModel80 *pmodel80;
      pmodel80 = new CModel80;
      pmod= (CModelBase*) pmodel80;
      return pmod;
      break;

    // -----------------------------------
    // |    REGISTER NEW DENSITIES HERE    |
    // -----------------------------------

  default : 
    std::cout << "Model ID out of range: model 1 used by default." << std::endl;
	pmodel1 = new CModel01;
	pmod= (CModelBase*) pmodel1;
    return pmod;
    break;
  }
}

moddefparam::moddefparam(std::string kw0,std::string def0,std::string desc0,std::string units0)
{
	keyword=kw0; 
	def=def0; 
	desc=desc0; 	
	units=units0;
}

