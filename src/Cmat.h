/*! \file Cmat.h
 * \brief 3 x 3 matrices with various useful methods.
 *
 */

#ifndef _CMAT_H_
#define _CMAT_H_

#include <iostream>
#include <cmath>
#include "Cvec.h"

class Cvec;

//! 3 x 3 matrix manipulation and operations
class Cmat {
    friend class Cvec;

public:
    Cmat();
    //! Copy.
    Cmat(const Cmat &a);
    //! Init with a 3x3 array.
    Cmat(const float [3][3]);
    //! 
    Cmat(const float);
    Cmat(float,float,float,float,float,float,float,float,float);
    Cmat(const float,const int);

    // member function declarations...
    Cvec column(int col);
    Cvec row(int row);
    Cmat& rotmat(float a,int ax);
    float determinant();
    Cmat tranpose();
    float subdet(int i,int j);
    Cmat const inverse();

    friend std::ostream& operator << (std::ostream&,const Cmat&);
    Cmat& operator = (const Cmat &a);
    Cmat operator * (Cmat b);
    friend Cmat operator * (const Cmat &a,const Cmat &b);
    friend Cvec operator * (const Cmat &m,const Cvec &v);
    Cvec operator [](int col);
    bool operator ==(const Cmat &a) const;

public:
    // Cmat variables
    float m[3][3];

};


#endif  //_CMAT_H_

