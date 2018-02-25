//
// File: CControlPoint.h
// $Id: CControlPoint.h,v 1.3 2009/02/09 20:50:58 thernis Exp $

#ifndef _CCONTROLPOINT_H_
#define _CCONTROLPOINT_H_

#include <iostream>
#include "Cvec.h"

//! Control point class for the Nurbs
class CControlPoint
{
	public:
		CControlPoint();
		 ~CControlPoint();
	
		// CControlPoint interface	
	
		CControlPoint(Cvec p0,float w0);
		CControlPoint(float a,float b,float c,float w0);

	Cvec p; //!< coordinate of the control point
	float w; //!< weight of the control point
	
	protected:
		// CControlPoint variables
	
		// TODO: add member variables...
	
};


#endif	//_CCONTROLPOINT_H_


/*
* $Log: CControlPoint.h,v $
* Revision 1.3  2009/02/09 20:50:58  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.2  2007/05/14 17:19:39  thernis
* Add CVS id and log in all files
*
*/

