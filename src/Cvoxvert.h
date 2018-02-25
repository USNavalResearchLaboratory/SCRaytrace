//
// $Id: Cvoxvert.h,v 1.3 2009/02/09 20:51:09 thernis Exp $

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


/*
* $Log: Cvoxvert.h,v $
* Revision 1.3  2009/02/09 20:51:09  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.2  2007/05/14 17:19:40  thernis
* Add CVS id and log in all files
*
*/
