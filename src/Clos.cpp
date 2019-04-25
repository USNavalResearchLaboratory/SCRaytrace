
#include "Clos.h"


Clos::Clos()
{
	// TODO: put constructor code here
}


Clos::~Clos()
{
	// TODO: put destructor code here
}

//! Common constructor 
/*!
      \param n Number of steps along the LOS
      \param st Start distance from reference point (either Obs or plane of sky)
      \param en End distance from reference point (either Obs or plane of sky)
      \sa setLOS()
    */
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

    this->nbp = n;
    this->ds = (en-st) / nbp;
    this->sstart = st + ds/2;
    this->send = en - ds/2;

}

