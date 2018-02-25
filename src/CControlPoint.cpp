//
// File: CControlPoint.cc
// $Id: CControlPoint.cpp,v 1.3 2009/02/09 20:50:58 thernis Exp $

#include "CControlPoint.h"

CControlPoint::CControlPoint()
{
	p=Cvec(0,0,0);
	w=0;

}


CControlPoint::~CControlPoint()
{
	// TODO: put destructor code here
}

CControlPoint::CControlPoint(Cvec p0,float w0)
{
	p=p0;
	w=w0;
}

CControlPoint::CControlPoint(float a,float b,float c,float w0)
{
	p=Cvec(a,b,c);
	w=w0;
}


/*
* $Log: CControlPoint.cpp,v $
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
