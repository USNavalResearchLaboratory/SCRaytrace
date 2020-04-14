
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

