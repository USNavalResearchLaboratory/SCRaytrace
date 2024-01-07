
#include <cmath>
#include "camera.h"
#include "Cvec.h"

Detector::Detector()
{
  sxpix=128;sypix=128;
  sxmm=1.;symm=1.;
}

Detector::Detector(const unsigned int &sxpix0,const unsigned int &sypix0,const float &sxmm0,const float &symm0)
{
  sxpix=sxpix0;sypix=sypix0;
  sxmm=sxmm0;symm=symm0;
}

Detector::Detector(const unsigned int &sxpix0,const unsigned int &sypix0)
{
  sxpix=sxpix0;sypix=sypix0;
  sxmm=1.;symm=1.;
}

Detector::Detector(const Detector &detector0)
{
  sxpix=detector0.sxpix;sypix=detector0.sypix;
  sxmm=detector0.sxmm;symm=detector0.symm;
}

const Detector &Detector::operator=(const Detector &detector0)
{
  if (&detector0 != this) {
    sxpix=detector0.sxpix;sypix=detector0.sypix;
    sxmm=detector0.sxmm;symm=detector0.symm;
  }
  return *this;
}

bool Detector::operator==(const Detector &a) const
{
  return (sxpix==a.sxpix && sypix==a.sypix && sxmm==a.sxmm && symm==a.symm);
}


void Detector::setSizePix(const unsigned int &sxpix0,const unsigned int &sypix0)
{
  sxpix=sxpix0;sypix=sypix0; 
}
    
unsigned int Detector::getSizePixX() const
{
  return sxpix;
}
unsigned int Detector::getSizePixY() const
{
  return sypix;
}
void Detector::setSizemm(const float &sxmm0,const float &symm0)
{
  sxmm=sxmm0;symm=symm0;
}
float Detector::getSizemmX() const
{
  return sxmm;
}
float Detector::getSizemmY() const
{
  return symm;
}

