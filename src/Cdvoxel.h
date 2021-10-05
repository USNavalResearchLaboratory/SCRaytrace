//
// File: Cdvoxel.h
//

#ifndef _CDVOXEL_H_
#define _CDVOXEL_H_

#include <iostream>
#include <vector>
#include "Cvec.h"
#include "Cmat.h"
#include "Cbasis.h"
#include "Cvoxvert.h" 

const int POSFACENORMAL[8][3]={{1,3,4},
      {0,2,5},
      {1,3,6},
      {0,2,7},
      {0,5,7},
      {1,4,6},
      {2,5,7},
      {3,4,6}}; //!< pos of the 3 adjacent face normal vector from the closest vert
const int FACEID[8][3]={{2,3,1},
        {4,3,1},
        {5,4,1},
        {5,2,1},
        {6,2,3},
        {6,4,3},
        {6,5,4},
        {6,5,2}}; //! the 3 adjacent face id number
const int BACKFACEID[8][3]={{4,5,6},
          {2,5,6},
          {3,2,6},
          {3,4,6},
          {1,4,5},
          {1,2,5},
          {1,3,2},
          {1,3,4}}; //! the 3 backside faces
const int COORDID[8][3][2]={{{1,2},{0,2},{0,1}},
            {{1,2},{0,2},{0,1}},
            {{0,2},{1,2},{0,1}},
            {{0,2},{1,2},{0,1}},
            {{0,1},{1,2},{0,2}},
            {{0,1},{1,2},{0,2}},
            {{0,1},{0,2},{1,2}},
            {{0,1},{0,2},{1,2}}}; //! 2 axis id of the 3 adjacent faces
const int OPPOSITEVERT[8]={6,7,4,5,2,3,0,1};

const Cvec ADJEDGE[8][3]={{Cvec(1,0,0),Cvec(0,1,0),Cvec(0,0,1)}, // 1,3,4
  {Cvec(1,0,0),Cvec(0,1,0),Cvec(0,0,1)},  // 0,2,5
  {Cvec(0,1,0),Cvec(1,0,0),Cvec(0,0,1)},  // 1,3,6
  {Cvec(0,1,0),Cvec(1,0,0),Cvec(0,0,1)},  // 0,2,7
  {Cvec(0,0,1),Cvec(1,0,0),Cvec(0,1,0)},  // 0,5,7
  {Cvec(0,0,1),Cvec(1,0,0),Cvec(0,1,0)},  // 1,4,6
  {Cvec(0,0,1),Cvec(0,1,0),Cvec(1,0,0)},  // 2,5,7
  {Cvec(0,0,1),Cvec(0,1,0),Cvec(1,0,0)}}; // 3,4,6


class Cvoxvert;

const Cvoxvert voxvert;

//! Voxel class for splatting integration
class Cdvoxel {
public:
    Cdvoxel();
    virtual ~Cdvoxel();

    //! Cdvoxel interfaces
    Cdvoxel(Cbasis b,float s,float d,Cbasis o);
    Cdvoxel(Cvec cntr,Cmat u,float s,float d,Cbasis o);
        
    //! Compute the LOSes for the vertices from the obs POV, origin on a vertice of the voxel
    void CalcVertLOS();
    //! Recalculate the LOSes for the vertices if the origin is changed.
    void UpdateVertLOS();
    //! Compute the density for a given LOS, do the work in the vox coord sys
    int GetThruLength(const Cvec &vlosobs,Cvec &vcntrbase,float &thrulength,const int &);
    //! Compute the density for a LOS passing through the voxel center
    void GetThruLengthCenter(float &thrulength,const int &);
    //! Find the closest vertice from the obs POV
    void CalcClosestVert();
    //! Find the closest vertice from the obs POV, different method
    void CalcClosestVert2();
    //! Update the position of the base cood sys origin and recompute closest vert
    void UpdateOriginandDensity(const Cvec&,const float&);
    //! Display the vertice position
    const void outputVert();
    //! Assigment operator
    const Cdvoxel &operator=(const Cdvoxel&);
    //! Print voxel parameters
    friend std::ostream& operator << (std::ostream&,const Cdvoxel&);

    Cbasis base; //!< position and orientation of the voxel with respect to absolute coord sys
    float side; //!< side lenght of the voxel, in Rsun
    float volume; //!< volume of the voxel, in cm^3
    Cvec voxcentervox; //! center of the voxel in voxel coordinates
    Cvec voxcenterobs0; //! initial center of the voxel in obs coord, with no translation
    Cvec voxcenterobs; //! center of the voxel in obs coord
    float dens;	//!< density in the voxel
    Cbasis obs; //!< pos and ori of the observer
    std::vector<Cvec> vertobs0; //!< pos of the initial 8 vertices, with no translation
    std::vector<Cvec> vertobs; //!< pos of the 8 vertices with respect to obs
    Cvec vtrans; //! pos of the vox origin in obs coord: useful to speed up vertice recomputation
    int posclosestvert; //!< points to the closest vertice from the obs pos
    //float vertmindist2obs; //!< dist from obs to closest vertice
    Cvec obsovox;
    Cvec obsovoxnorm;

  protected:
    // Cdvoxel variables

    // TODO: add member variables...

};


#endif	//_CDVOXEL_H_


