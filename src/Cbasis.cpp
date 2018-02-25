//
// $Id: Cbasis.cpp,v 1.10 2009/04/13 20:52:31 thernis Exp $


#include "Cbasis.h"

Cbasis::Cbasis()
{
	o=Cvec(0,0,0);
  u=Cmat(1,0,0,0,1,0,0,0,1);
	ui=u;

}

Cbasis::Cbasis(const Cbasis &a)
{
    o=a.o;
    u=a.u;
    ui=a.ui;
}

Cbasis::Cbasis(const Cvec &cntr, const float &a, const float &b, const float &c)
{
    o=cntr;

    u=Cmat(c,3)*Cmat(b,2)*Cmat(a,1);
    ui=Cmat(-a,1)*Cmat(-b,2)*Cmat(-c,3);

}


Cbasis::Cbasis(const Cvec &cntr, const float &a, const float &b, const float &c,
               const float &Hlon,const float &Hlat,const float &Hrot)
{
    o=cntr;

    u=Cmat(c,3)*Cmat(b,2)*Cmat(a,1)*Cmat(Hlat,2)*Cmat(-Hlon,1)*Cmat(-Hrot,3);
    ui=Cmat(Hrot,3)*Cmat(Hlon,1)*Cmat(-Hlat,2)*Cmat(-a,1)*Cmat(-b,2)*Cmat(-c,3);

}



Cbasis::Cbasis(const Cvec &cntr,Cmat mat) {
    o=cntr;
    u=mat;
    ui=mat.inverse();

}


Cbasis::Cbasis (const float &lon,
        const float &lat,
        const float &height,
        const float &a,const float &b,const float &c) 
{
  u=Cmat(c,3)*Cmat(b,2)*Cmat(a,1)*Cmat(-lat,2)*Cmat(lon+PI,1);
  ui=Cmat(-lon-PI,1)*Cmat(lat,2)*Cmat(-a,1)*Cmat(-b,2)*Cmat(-c,3);

  o=Cmat(-lon-PI,1)*Cmat(lat,2)*Cvec(0.,0.,-height);

}



void Cbasis::setRotationPerAxis(const float &angA,const unsigned int &axisA,
                                const float &angB,const unsigned int &axisB,
                                const float &angC,const unsigned int &axisC)
{
  u=Cmat(angA,axisA)*Cmat(angB,axisB)*Cmat(angC,axisC);
  ui=Cmat(-angC,axisC)*Cmat(-angB,axisB)*Cmat(-angA,axisA);
}


Cbasis& Cbasis::operator = (const Cbasis &a) {
    o=a.o;
    u=a.u;
    ui=a.ui;

    return *this;
}




std::ostream& operator << (std::ostream& os,const Cbasis& c) {
  os << std::endl
      << "O : " << c.o << std::endl
      << "U : " << c.u << std::endl;

  return os;
}




Cvec ChangeBase(const Cvec &v,const Cbasis &baseinit,const Cbasis &basedest) { 
  return (basedest.u*(baseinit.o-basedest.o+(baseinit.ui*v)));

//  return ((basedest.u*(baseinit.o-basedest.o+(baseinit.ui*(v+baseinit.translation))))-basedest.translation);

}


/*
* $Log: Cbasis.cpp,v $
* Revision 1.10  2009/04/13 20:52:31  thernis
* - Cleanup not useful stuff.
* - Implement method to allow setting user defined rotation order.
*
* Revision 1.9  2009/03/06 21:21:19  thernis
* Add a translation within the coordinate system, not in the absolute coordinate system
*
* Revision 1.8  2009/02/09 20:51:02  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.7  2007/07/20 14:33:41  thernis
* Got it, there is a 180 deg shift, but not on Hlon
*
* Revision 1.6  2007/07/20 14:11:25  thernis
* Move back to no 180 deg shift: I'm just crazy
*
* Revision 1.5  2007/07/19 22:20:09  thernis
* Shift longitude for Carrington positionning so that the z axis of abs system points to Carrington longitude and latitude 0
*
* Revision 1.4  2007/07/19 19:27:04  thernis
* Change rotation order for the Carrington longitude and latitude
*
* Revision 1.3  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
