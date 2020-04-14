/*! \file Cvoxvert.h
 * \brief Voxel vertex.
 *
 */

#ifndef CVOXVERT_H
#define CVOXVERT_H
#include <vector>
#include "Cvec.h"

//! Voxel vertice position, in voxel coord system
class Cvoxvert{
public:
    Cvoxvert();

    ~Cvoxvert();
    Cvec operator [](const int &i) const;

    std::vector<Cvec> vert;
    
};

#endif

