
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
}



