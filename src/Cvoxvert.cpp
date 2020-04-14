
#include "Cvoxvert.h"
#include <cassert>
Cvoxvert::Cvoxvert()
{
      // -- vert 0
  vert.push_back(Cvec(0,0,0));
    // -- vert 1
  vert.push_back(Cvec(1,0,0));
    // -- vert 2
  vert.push_back(Cvec(1,1,0));
    // -- vert 3
  vert.push_back(Cvec(0,1,0));
    // -- vert 4
  vert.push_back(Cvec(0,0,1));
    // -- vert 5
  vert.push_back(Cvec(1,0,1));
    // -- vert 6
  vert.push_back(Cvec(1,1,1));
    // -- vert 7
  vert.push_back(Cvec(0,1,1));

}


Cvoxvert::~Cvoxvert()
{
  vert.clear();
}

Cvec Cvoxvert::operator [](const int &i) const {
    assert(i>=0 && i<8);
    return vert[i];
}
