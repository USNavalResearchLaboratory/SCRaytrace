
#include <cmath>
#include "sun.h"
#include "constant.h"




Sun::Sun()
{
  radius=696000e3;
  u=0.7;
  basis=new Cbasis(Cvec(0.,0.,0.), 0., 0., 0.);
  calcConstfactor();
}


Sun::~Sun()
{
  delete basis;
}



float Sun::getRadius()
{
  return radius;
}

float Sun::getLimbDarkening()
{
  return u;
}

void Sun::setLimbDarkening(float limbdarkeningcoeff)
{
  this->u=limbdarkeningcoeff;
  calcConstfactor();
}

Cbasis* Sun::getPosition()
{
  return basis;
}

