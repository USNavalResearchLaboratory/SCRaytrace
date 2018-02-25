//
// File: Cbasis.h
// $Id: Cbasis.h,v 1.6 2009/04/13 20:52:31 thernis Exp $

#ifndef _CBASIS_H_
#define _CBASIS_H_

#include <iostream>
#include "constant.h"
#include "Cvec.h"
#include "Cmat.h"

//! Coordinate system definition and passage matrices.
class Cbasis {
public:
    //! Default constructor.
    Cbasis();
    //! Copy.
    Cbasis(const Cbasis &a);
    //! Definition of center and the 3 rotations around XYZ.
    //! The rotation order is Z then Y then X.
    Cbasis(const Cvec &cntr, const float &a, const float &b, const float &c);



    //! Center, the 3 rotations around XYZ, heliographic Lon and Lat projection, position angle of the North pole  counterclockwise
    Cbasis(const Cvec &cntr, const float &a, const float &b, const float &c,
               const float &Hlon,const float &Hlat,const float &Hrot);

    //! Center and direct tranform matrix
    Cbasis(const Cvec &cntr,Cmat mat);
    //! Carrington longitude, latitude and height: useful for defining the obs position
    Cbasis (const float &lon,
            const float &lat,
            const float &height,
            const float &a,const float &b,const float &c);

    //! set the origin of the basis
    void setCenter(const Cvec &cntr) {this->o=cntr;};

    //! set the rotation matrix using rotation angles and corresponding axis
    //! \param ang# angle of rotation, in radian, corresponding to the axis #
    //! \param axis# should be 1, 2 or 3. 1: axis X, 2: axis Y, 3: axis Z
    void setRotationPerAxis(const float &angA,const unsigned int &axisA,
                            const float &angB,const unsigned int &axisB,
                            const float &angC,const unsigned int &axisC);

		//! Change coordinates of \a v from \a baseinit coordinate system to \a basedest coordinate system
    friend Cvec ChangeBase(const Cvec &v,const Cbasis &baseinit,const Cbasis &basedest);
    
    //! Print 
    friend std::ostream& operator << (std::ostream&,const Cbasis&);

    Cbasis& operator = (const Cbasis &a);
public:
    Cvec o; //!< origin of the coordinate system, always expressed in absolute coordinates
    Cmat u; //!< coord transform matrix from absolute to this coordinate system
    Cmat ui; //!< coord transform matrix from this coordinate system to absolute

};


std::ostream& operator << (std::ostream& os,const Cbasis& c);
Cvec ChangeBase(const Cvec &v,const Cbasis &baseinit,const Cbasis &basedest);

#endif	//_CBASIS_H_


/*
* $Log: Cbasis.h,v $
* Revision 1.6  2009/04/13 20:52:31  thernis
* - Cleanup not useful stuff.
* - Implement method to allow setting user defined rotation order.
*
* Revision 1.5  2009/03/06 21:21:19  thernis
* Add a translation within the coordinate system, not in the absolute coordinate system
*
* Revision 1.4  2009/02/09 20:51:02  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.3  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
