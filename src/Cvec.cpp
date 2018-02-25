//
// File: Cvec.cc
// $Id: Cvec.cpp,v 1.5 2009/02/09 20:51:08 thernis Exp $

#include "Cvec.h"
#include <cassert>

Cvec::Cvec() {
    v[0]=0.;
    v[1]=0.;
    v[2]=0.;
}

Cvec::Cvec(float *u) {
    v[0]=u[0];
    v[1]=u[1];
    v[2]=u[2];
}

Cvec::Cvec(const float &x,const float &y,const float &z) {
  v[0]=x;
  v[1]=y;
  v[2]=z;
}

Cvec::Cvec(const Cvec &a) {
    v[0]=a.v[0];
    v[1]=a.v[1];
    v[2]=a.v[2];
}

// ---- op overload
std::ostream& operator << (std::ostream& os,const Cvec& u) {
    os << u.v[0] << " , " << u.v[1] << " , " << u.v[2];
    return os;
}

float Cvec::norm() const {
    return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

float Cvec::mag() const {
    return(sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

float Cvec::magsqr() const {
  return(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
}

void Cvec::leftshift() {
    float tmp=v[0];
    v[0]=v[1];
    v[1]=v[2];
    v[2]=tmp;
}

void Cvec::rightshift() {
    float tmp=v[0];
    v[0]=v[2];
    v[2]=v[1];
    v[1]=tmp;
}
    
bool Cvec::IsNull() {
  return (v[0]==0. && v[1]==0. && v[2]==0.); 
}

bool Cvec::operator ==(const Cvec &b) const {
  return (v[0]==b.v[0] && v[1]==b.v[1] && v[2]==b.v[2]);
}
  
Cvec& Cvec::operator = (const Cvec &u) {
    v[0]=u.v[0];
    v[1]=u.v[1];
    v[2]=u.v[2];

    return *this;
}

Cvec Cvec::operator + (const Cvec &b) const {
    Cvec a;
    a.v[0]=v[0]+b.v[0];
    a.v[1]=v[1]+b.v[1];
    a.v[2]=v[2]+b.v[2];

    return a;
}

Cvec Cvec::operator - (const Cvec &b) const {
    Cvec a;
    a.v[0]=v[0]-b.v[0];
    a.v[1]=v[1]-b.v[1];
    a.v[2]=v[2]-b.v[2];
    return a;
}

Cvec Cvec::operator + (const float &a) {
  Cvec b;
  b.v[0]=v[0]+a;
  b.v[1]=v[1]+a;
  b.v[2]=v[2]+a;
  return b;
}

Cvec Cvec::operator - (const float &a) {
  Cvec b;
  b.v[0]=v[0]-a;
  b.v[1]=v[1]-a;
  b.v[2]=v[2]-a;
  return b;
}

Cvec Cvec::operator - () {
    Cvec b;
    b.v[0]=-v[0];
    b.v[1]=-v[1];
    b.v[2]=-v[2];
    return b;
}

Cvec Cvec::operator * (const float &a) const {
    Cvec b;
    b.v[0]=v[0]*a;
    b.v[1]=v[1]*a;
    b.v[2]=v[2]*a;
    return b;
}

Cvec operator * (float a,Cvec v) {
    Cvec b;
    b.v[0]=v.v[0]*a;
    b.v[1]=v.v[1]*a;
    b.v[2]=v.v[2]*a;
    return b;
}

Cvec Cvec::operator / (const float &a) const {
    Cvec b;
    b.v[0]=v[0]/a;
    b.v[1]=v[1]/a;
    b.v[2]=v[2]/a;
    return b;
}



Cvec& Cvec::operator += (Cvec b) {
    v[0]+=b.v[0];
    v[1]+=b.v[1];
    v[2]+=b.v[2];
    return *this;
}

float Cvec::operator [](const int &i) const {
    assert(i<3);
    return v[i];
}




/*
* $Log: Cvec.cpp,v $
* Revision 1.5  2009/02/09 20:51:08  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.4  2007/07/19 19:32:39  thernis
* Implement magsqr function
*
* Revision 1.3  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
