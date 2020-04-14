/*! \file Cvec.h 
 * \brief 3 elements vector class, with method to act on them.
 *
 */

#ifndef _CVEC_H_
#define _CVEC_H_

#include <ostream>
#include "Cmat.h"

//! 3 element vector.
class Cvec {
    friend class Cmat;
public:
    Cvec();
    //! Init with an array of 3 elements.
    Cvec(float*);
    //! Init with 3 elements.
    Cvec(const float &,const float &,const float &);
    //! Copy.
    Cvec(const Cvec &a);
    //! Print.
    friend std::ostream& operator << (std::ostream&,const Cvec&);
    float norm() const; //!< return magnitude
    float mag() const; //!< return magnitude
    float magsqr() const; //!< return square magnitude (no square root applied)
    void leftshift(); //!< [a,b,c] -> [b,c,a]
    void rightshift(); //!< [a,b,c] -> [c,a,b]
    bool IsNull(); //! test if null vector
    Cvec& operator = (const Cvec &u);
    Cvec operator + (const Cvec &b) const;
    Cvec operator + (const float &a);
    Cvec operator - (const float &a);
    Cvec operator - (const Cvec &b) const;
    Cvec operator - ();
    Cvec operator * (const float &a) const;
    friend Cvec operator * (float a,Cvec v);
    Cvec operator / (const float &a) const;
    Cvec& operator += (Cvec b);
    float operator [](const int &i) const;
    bool operator ==(const Cvec &b) const;
    
    //! vector product
    friend Cvec pvect(Cvec a,Cvec b);
    //! distance point - straight line
    friend float psldist(Cvec psl,Cvec vsl,Cvec p);
    //! scalar product
    friend float pscal(Cvec a,Cvec b);
    //! Orthogonal projection of a point on a straight line
    friend Cvec orthoproj(Cvec psl,Cvec vsl,Cvec p);
    //! Orthogonal projection of a point on a plane
    friend Cvec orthoprojpointplane(Cvec pplane,Cvec vnormal,Cvec p);
    //! Intersection SL with a plane
    friend Cvec interslplane(const Cvec &psl,const Cvec &vsl,const Cvec &pplane,const Cvec &vnplane);


public:
    // Cvec variables
    float v[3];
};


//! vector product
inline Cvec pvect(Cvec a,Cvec b) {
    Cvec p;
    p.v[0]=a.v[1]*b.v[2]-a.v[2]*b.v[1];
    p.v[1]=a.v[2]*b.v[0]-a.v[0]*b.v[2];
    p.v[2]=a.v[0]*b.v[1]-a.v[1]*b.v[0];

    return p;
}

//! distance point - straight line
inline float psldist(Cvec psl,Cvec vsl,Cvec p) {
    return(pvect(psl-p,vsl).norm() / vsl.norm());
}

//! scalar product
inline float pscal(Cvec a,Cvec b) {
    return(a.v[0]*b.v[0]+a.v[1]*b.v[1]+a.v[2]*b.v[2]);
}

//! Orthogonal projection of a point on a straight line
inline Cvec orthoproj(Cvec psl,Cvec vsl,Cvec p) {
    Cvec r;
    float m;
    vsl=vsl/vsl.norm();
    m=pscal(p-psl,vsl);
    r=psl+(vsl * m);
    return r;
}

//! Orthogonal projection of a point on a plane
inline Cvec orthoprojpointplane(Cvec pplane,Cvec vnormal,Cvec p) {
    Cvec r;
    float m;
    vnormal=vnormal/vnormal.norm();
    m=pscal(vnormal,p-pplane);
    r=p-(vnormal * m);
    return r;
}

//! Intersection SL with a plane
inline Cvec interslplane(const Cvec &psl,const Cvec &vsl,const Cvec &pplane,const Cvec &vnplane) {

    float k=pscal(pplane-psl,vnplane)/pscal(vsl,vnplane);

    return psl+k*vsl;

}


//! vector product
//Cvec pvect(Cvec a,Cvec b);
//! distance point - straight line
//float psldist(Cvec psl,Cvec vsl,Cvec p);
//! scalar product
//float pscal(Cvec a,Cvec b);
//! Orthogonal projection of a point on a straight line
//Cvec orthoproj(Cvec psl,Cvec vsl,Cvec p);
//! Orthogonal projection of a point on a plane
//Cvec orthoprojpointplane(Cvec pplane,Cvec vnormal,Cvec p);
//! Intersection SL with a plane
//Cvec interslplane(const Cvec &psl,const Cvec &vsl,const Cvec &pplane,const Cvec &vnplane);

#endif	//_CVEC_H_


