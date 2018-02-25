//
// File: Clos.cc
// $Id: Clos.cpp,v 1.3 2009/02/09 20:51:04 thernis Exp $

#include "Clos.h"


Clos::Clos()
{
	// TODO: put constructor code here
}


Clos::~Clos()
{
	// TODO: put destructor code here
}

Clos::Clos(int n,float st,float en) {
setLOS(n,st,en);

  }


Clos::Clos(const Clos &a)
{
   nbp=a.nbp;
   ds=a.ds;
   sstart=a.sstart;
   send=a.send;
}


void Clos::setLOS(int n,float st,float en)
{

    this->nbp=n;
    this->ds=(en-st)/(nbp);
    this->sstart=st+ds/2;
    this->send=en-ds/2;

}


  
/*
  * $Log: Clos.cpp,v $
  * Revision 1.3  2009/02/09 20:51:04  thernis
  * - Clean up the code
  * - Change CModel::Density prototype
  * - Update documentation
  * - Implement multi-threading using boost thread
  * - Add new models
  *
  * Revision 1.2  2007/05/14 17:19:40  thernis
  * Add CVS id and log in all files
  *
*/
