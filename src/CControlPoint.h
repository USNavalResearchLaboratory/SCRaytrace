/*! \file CControlPoint.h
 * \brief Control point class for NURBS.
 *
 */

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

