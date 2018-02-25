#include <cmath>
#include "camera.h"
#include "Cvec.h"

CCD::CCD()
{
  sxpix=0;sypix=0;
  sxmm=1.;symm=1.;
}

CCD::CCD(const unsigned int &sxpix0,const unsigned int &sypix0,const float &sxmm0,const float &symm0)
{
  sxpix=sxpix0;sypix=sypix0;
  sxmm=sxmm0;symm=symm0;
}

CCD::CCD(const unsigned int &sxpix0,const unsigned int &sypix0)
{
  sxpix=sxpix0;sypix=sypix0;
}

CCD::CCD(const CCD &ccd0)
{
  sxpix=ccd0.sxpix;sypix=ccd0.sypix;
  sxmm=ccd0.sxmm;symm=ccd0.symm;
}

const CCD &CCD::operator=(const CCD &ccd0)
{
  if (&ccd0 != this) {
    sxpix=ccd0.sxpix;sypix=ccd0.sypix;
    sxmm=ccd0.sxmm;symm=ccd0.symm;
  }
  return *this;
}

bool CCD::operator==(const CCD &a) const
{
  return (sxpix==a.sxpix && sypix==a.sypix && sxmm==a.sxmm && symm==a.symm);
}


void CCD::setSizePix(const unsigned int &sxpix0,const unsigned int &sypix0)
{
  sxpix=sxpix0;sypix=sypix0; 
}
    
unsigned int CCD::getSizePixX() const
{
  return sxpix;
}
unsigned int CCD::getSizePixY() const
{
  return sypix;
}
void CCD::setSizemm(const float &sxmm0,const float &symm0)
{
  sxmm=sxmm0;symm=symm0;
}
float CCD::getSizemmX() const
{
  return sxmm;
}
float CCD::getSizemmY() const
{
  return symm;
}

