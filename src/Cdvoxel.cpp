//
// File: Cdvoxel.cc
//

#include "Cdvoxel.h"
#include <cmath>
#include <cassert>
#include "constant.h"


Cdvoxel::Cdvoxel() {}

Cdvoxel::Cdvoxel(Cbasis b,float s,float d,Cbasis o) {
    base=b;
    side=s;
    volume=side*side*side;
    dens=d;
    obs=o;
    voxcentervox=Cvec(side/2,side/2,side/2);
    obsovox=base.u*(obs.o-base.o);
    obsovoxnorm=obsovox/side;
    CalcVertLOS();
    CalcClosestVert2();
}

Cdvoxel::Cdvoxel(Cvec cntr,Cmat u,float s,float d,Cbasis o) {
    base=Cbasis(cntr,u);
    side=s;
    volume=side*side*side;
    dens=d;
    obs=o;
    voxcentervox=Cvec(side/2,side/2,side/2);
    obsovox=base.u*(obs.o-base.o);
    obsovoxnorm=obsovox/side;
    CalcVertLOS();
    CalcClosestVert2();
}

const Cdvoxel &Cdvoxel::operator=(const Cdvoxel &a) {
  if (&a != this) {
    base=a.base;
    side=a.side;
    volume=a.volume;
    dens=a.dens;
    obs=a.obs;
    obsovox=a.obsovox;
    obsovoxnorm=a.obsovoxnorm;
    voxcentervox=a.voxcentervox;
    voxcenterobs0=a.voxcenterobs0;
    voxcenterobs=a.voxcenterobs;
    posclosestvert=a.posclosestvert;
    vertobs0=a.vertobs0;
    vertobs=a.vertobs;
    vtrans=a.vtrans;
  }
  return *this;
}


void Cdvoxel::CalcVertLOS() {
    // ---- Pos of the vert in the Obs coord sys
    vertobs0.clear();
    vtrans=obs.u*base.o;
    Cvec voobs=obs.u*obs.o;
    Cmat Mbase2obs=obs.u*base.ui;
    
    // -- vert 0
    //vertobs0.push_back(vtrans-obs.u*obs.o);
    //vertobs0.push_back(obs.u*(base.o-obs.o));
    vertobs0.push_back(voobs);
    // -- vert 1
    //vertobs0.push_back(obs.u*((base.o+base.ui*Cvec(side,0,0))-obs.o));
    vertobs0.push_back(Mbase2obs*Cvec(side,0,0)-voobs);
    // -- vert 2
    //vertobs0.push_back(obs.u*((base.o+base.ui*Cvec(side,side,0))-obs.o));
    vertobs0.push_back(Mbase2obs*Cvec(side,side,0)-voobs);
   // -- vert 3
    //vertobs0.push_back(obs.u*((base.o+base.ui*Cvec(0,side,0))-obs.o));
    vertobs0.push_back(Mbase2obs*Cvec(0,side,0)-voobs);
    // -- vert 4
    //vertobs0.push_back(obs.u*((base.o+base.ui*Cvec(0,0,side))-obs.o));
    vertobs0.push_back(Mbase2obs*Cvec(0,0,side)-voobs);
   // -- vert 5
    //vertobs0.push_back(obs.u*((base.o+base.ui*Cvec(side,0,side))-obs.o));
    vertobs0.push_back(Mbase2obs*Cvec(side,0,side)-voobs);
   // -- vert 6
    //vertobs0.push_back(obs.u*((base.o+base.ui*Cvec(side,side,side))-obs.o));
    vertobs0.push_back(Mbase2obs*Cvec(side,side,side)-voobs);
    // -- vert 7
    //vertobs0.push_back(obs.u*((base.o+base.ui*Cvec(0,side,side))-obs.o));
    vertobs0.push_back(Mbase2obs*Cvec(0,side,side)-voobs);
    
    voxcenterobs0=Mbase2obs*voxcentervox-voobs;

    // -- translated verlos
    vertobs.clear();
    for (unsigned int i=0;i<vertobs0.size();i++) vertobs.push_back(vertobs0[i]+vtrans);
    voxcenterobs=voxcenterobs0+vtrans;


}
    
void Cdvoxel::UpdateVertLOS() {
    // ---- recalc pos of the vert in the Obs coord sys
  vtrans=obs.u*base.o;
  for (unsigned int i=0;i<vertobs0.size();i++) vertobs[i]=vertobs0[i]+vtrans;
  voxcenterobs=voxcenterobs0+vtrans;

  /*
  Cvec vtransnew=obs.u*base.o;
  Cvec vt=vtransnew-vtrans;
  for (int i=0;i<vertobs.size();i++) vertobs[i]+=vt;
  voxcenterobs+=vt;
  vtrans=vtransnew;
  */
}
    

//! Compute the density for a given LOS, do the work in the vox coord sys
int Cdvoxel::GetThruLength(const Cvec &vlosobs,Cvec &vcntrbase,float &thrulength,const int &displayflag ) {
  
  // -- compute vlos in vox coord
  Cvec vlosvox=base.u*obs.ui*vlosobs;
  
    // ---- Compute intersection of the LOS with the 3 adjacent faces of the closest vertice.
    //    and find in which one the LOS enters
  int faceenterpos=-1;
  Cvec PTinter[3]; //!< The intersection points with the 3 adj faces
  for (int i=0;i<3;i++) {

    PTinter[i]=interslplane(obsovoxnorm,vlosvox,
                            voxvert[posclosestvert],ADJEDGE[posclosestvert][i]);
    
    if (displayflag) {
      std::cout << "-----------------" << std::endl;
      std::cout << "PTinter["<<i<<"]   : " << PTinter[i] << std::endl;
      std::cout << "ADJEDGE[posclosestvert][i] : " << ADJEDGE[posclosestvert][i] << std::endl;
      std::cout << "FACEID : " << FACEID[posclosestvert][i] << std::endl;
      std::cout << "COORDID[posclosestvert][i][0,1] : " << COORDID[posclosestvert][i][0] <<", "<< COORDID[posclosestvert][i][1] << std::endl;
      
    }
    
    
    if (PTinter[i][COORDID[posclosestvert][i][0]] >= 0 && PTinter[i][COORDID[posclosestvert][i][0]] < 1 &&
        PTinter[i][COORDID[posclosestvert][i][1]] >= 0 && PTinter[i][COORDID[posclosestvert][i][1]] < 1) {
      faceenterpos=i;
      break; // exit the loop when the face is found
        }

  }
  
  if (displayflag) {
    std::cout << "vlosobs        : " << vlosobs << std::endl;
    std::cout << "vlosvox        : " << vlosvox << std::endl;
    std::cout << "obsovox        : " << obsovox << std::endl;
    std::cout << "posclosestvert : " << posclosestvert << std::endl;
    std::cout << "voxvert[posclosestvert]    : " << voxvert[posclosestvert] << std::endl;
    std::cout << "faceenterpos : " << faceenterpos << std::endl;
    if (faceenterpos >= 0)
      std::cout << "FACEID[faceenterpos] : " << FACEID[posclosestvert][faceenterpos] << std::endl;
  }
    
    // -- return if outside the voxel
  if (faceenterpos < 0)
    return 1;

    // ---- find where the LOS exits the voxel
  Cvec PTinterback[3]; //!< The intersection points with the 3 backside faces
  int oppvert=OPPOSITEVERT[posclosestvert];
  int faceexitpos=-1;
  for (int i=0;i<3;i++) {

    PTinterback[i]=interslplane(obsovoxnorm,vlosvox,
                                voxvert[oppvert],ADJEDGE[oppvert][i]);

    if (PTinterback[i][COORDID[oppvert][i][0]] >= 0 && PTinterback[i][COORDID[oppvert][i][0]] < 1 &&
        PTinterback[i][COORDID[oppvert][i][1]] >= 0 && PTinterback[i][COORDID[oppvert][i][1]] < 1) {
      faceexitpos=i;
      break; // exit the loop when the face is found
        }
  }

  if (faceexitpos < 0)
    return 1;

    // ---- compute central pos of the LOS inside the voxel
  vcntrbase=(PTinter[faceenterpos]+PTinterback[faceexitpos])*side/2;
 
    // ---- compute length of LOS inside the voxel
  thrulength=
      (PTinter[faceenterpos]-PTinterback[faceexitpos]).mag()*side;
  
  if (displayflag) {
    std::cout << "vcntrbase : " << vcntrbase << std::endl;
    std::cout << "thrulength : " << thrulength << std::endl;
  }
  
  return 0;
}
    

//! Compute the density for a LOS passing through the voxel center
void Cdvoxel::GetThruLengthCenter(float &thrulength,const int &displayflag) {
  // -- compute vlos in vox coord
  Cvec vlosvox=voxcentervox-obsovox;
      
  // ---- Compute intersection of the LOS with the 3 adjacent faces of the closest vertice.
  //    and find in which one the LOS enters
  int faceenterpos=-1;
  Cvec PTinter[3]; //!< The intersection points with the 3 adj faces
  for (int i=0;i<3;i++) {

    PTinter[i]=interslplane(obsovoxnorm,vlosvox,
                            voxvert[posclosestvert],ADJEDGE[posclosestvert][i]);
    
    if (displayflag) {
      std::cout << "-----------------" << std::endl;
      std::cout << " interslplane : " << interslplane(obsovoxnorm,vlosvox,
          voxvert[posclosestvert],ADJEDGE[posclosestvert][i]) << std::endl;
      std::cout << "obsovoxnorm  : " << obsovoxnorm << std::endl;
      std::cout << "PTinter["<<i<<"]   : " << PTinter[i] << std::endl;
      std::cout << "ADJEDGE[posclosestvert][i] : " << ADJEDGE[posclosestvert][i] << std::endl;
      std::cout << "FACEID : " << FACEID[posclosestvert][i] << std::endl;
      std::cout << "COORDID[posclosestvert][i][0,1] : " << COORDID[posclosestvert][i][0] <<", "<< COORDID[posclosestvert][i][1] << std::endl;
      
    }
    
    
    if (PTinter[i][COORDID[posclosestvert][i][0]] >= 0 && PTinter[i][COORDID[posclosestvert][i][0]] < 1 &&
        PTinter[i][COORDID[posclosestvert][i][1]] >= 0 && PTinter[i][COORDID[posclosestvert][i][1]] < 1) {
      faceenterpos=i;
      break; // exit the loop when the face is found
        }

  }
  
  thrulength=PTinter[faceenterpos].mag()*side*2;

  if (displayflag) {
    std::cout << "PTinter[faceenterpos]       : " << PTinter[faceenterpos] << std::endl;
    std::cout << "thrulength       : " << thrulength << std::endl;
    std::cout << "vlosvox        : " << vlosvox << std::endl;
    std::cout << "obsovox        : " << obsovox << std::endl;
    std::cout << "posclosestvert : " << posclosestvert << std::endl;
    std::cout << "voxvert[posclosestvert]    : " << voxvert[posclosestvert] << std::endl;
    std::cout << "faceenterpos : " << faceenterpos << std::endl;
    if (faceenterpos >= 0)
      std::cout << "FACEID[faceenterpos] : " << FACEID[posclosestvert][faceenterpos] << std::endl;
  }

  return;
}



void Cdvoxel::CalcClosestVert() {
    // ---- find the voxel vertice the closest from the obs
    float vertmindist2obs=vertobs[0].mag();
    posclosestvert=0;
    for (int i=1;i<8;i++) {
        float junk=vertobs[i].mag();
        if (junk < vertmindist2obs) {
            vertmindist2obs=junk;
            posclosestvert=i;
        }
    }
}


void Cdvoxel::CalcClosestVert2() {
    // ---- compute orientation of vector obs.o to voxel.o at the center of the voxelin the vox coord sys
  Cvec OoOv=base.u*(base.o+(base.ui*voxcentervox)-obs.o);
    
    // ---- find the largest coordinate absolute value
    int lorder[3]={0,1,2};
    int swap;
    if (fabs(OoOv[0]) < fabs(OoOv[1])) {
      lorder[0]=1;lorder[1]=0;
    }
    if (fabs(OoOv[lorder[0]]) < fabs(OoOv[2])) {
        swap=lorder[0];
        lorder[0]=2;
        lorder[2]=swap;
    }
    if (fabs(OoOv[lorder[1]]) < fabs(OoOv[lorder[2]])) {
        swap=lorder[1];
        lorder[1]=lorder[2];
        lorder[2]=swap;
    }

    //std::cout << "lorder : " << lorder[0] << ", " << lorder[1] << ", " << lorder[2] << std::endl;
    
    // ---- find the entrance side
    switch (lorder[0]) {

    case 0 :
        // -- face 2 or 4
        if (OoOv[0] > 0) {
            // -- face 2
            if (lorder[1] == 1) {
                if (OoOv[1] > 0)  {
                    // -- edge 0-4
                    if (OoOv[2] > 0)
                      posclosestvert=0;
                    else
                        posclosestvert=4;
                } else {
                    // -- edge 3-7
                    if (OoOv[2] > 0)
                        posclosestvert=3;
                    else
                        posclosestvert=7;
                }
            } else {
                // then lorder[1]==2
                if (OoOv[2] > 0)  {
                    // -- edge 0-3
                    if (OoOv[1] > 0)
                        posclosestvert=0;
                    else
                        posclosestvert=3;
                } else {
                    // -- edge 4-7
                    if (OoOv[1] > 0)
                        posclosestvert=4;
                    else
                        posclosestvert=7;
                }
            }
        } else {
            // -- then face 4
            if (lorder[1] == 1) {
                if (OoOv[1] > 0)  {
                    // -- edge 1-5
                    if (OoOv[2] > 0)
                        posclosestvert=1;
                    else
                        posclosestvert=5;
                } else {
                    // -- edge 2-6
                    if (OoOv[2] > 0)
                        posclosestvert=2;
                    else
                        posclosestvert=6;
                }
            } else {
                // then lorder[1]==2
                if (OoOv[2] > 0)  {
                    // -- edge 1-2
                    if (OoOv[1] > 0)
                        posclosestvert=1;
                    else
                        posclosestvert=2;
                } else {
                    // -- edge 5-6
                    if (OoOv[1] > 0)
                        posclosestvert=5;
                    else
                        posclosestvert=6;
                }
            }
        }
        break;

    case 1 :
        // -- face 3 or 5
        if (OoOv[1] > 0) {
            // -- face 3
            if (lorder[1] == 0) {
                if (OoOv[0] > 0)  {
                    // -- edge 0-4
                    if (OoOv[2] > 0)
                        posclosestvert=0;
                    else
                        posclosestvert=4;
                } else {
                    // -- edge 1-5
                    if (OoOv[2] > 0)
                        posclosestvert=1;
                    else
                        posclosestvert=5;
                }
            } else {
                // then lorder[1]==2
                if (OoOv[2] > 0)  {
                    // -- edge 0-1
                    if (OoOv[0] > 0)
                        posclosestvert=0;
                    else
                        posclosestvert=1;
                } else {
                    // -- edge 4-5
                    if (OoOv[0] > 0)
                        posclosestvert=4;
                    else
                        posclosestvert=5;
                }
            }
        } else {
            // -- then face 5
            if (lorder[1] == 0) {
                if (OoOv[0] > 0)  {
                    // -- edge 3-7
                    if (OoOv[2] > 0)
                        posclosestvert=3;
                    else
                        posclosestvert=7;
                } else {
                    // -- edge 2-6
                    if (OoOv[2] > 0)
                        posclosestvert=2;
                    else
                        posclosestvert=6;
                }
            } else {
                // then lorder[1]==2
                if (OoOv[2] > 0)  {
                    // -- edge 2-3
                    if (OoOv[0] > 0)
                        posclosestvert=3;
                    else
                        posclosestvert=2;
                } else {
                    // -- edge 6-7
                    if (OoOv[0] > 0)
                        posclosestvert=7;
                    else
                        posclosestvert=6;
                }
            }
        }

        break;
    case 2 :
        // -- face 1 or 6
        if (OoOv[2] > 0) {
            // -- face 3
            if (lorder[1] == 0) {
                if (OoOv[0] > 0)  {
                    // -- edge 0-3
                    if (OoOv[1] > 0)
                        posclosestvert=0;
                    else
                        posclosestvert=3;
                } else {
                    // -- edge 1-2
                    if (OoOv[1] > 0)
                        posclosestvert=1;
                    else
                        posclosestvert=2;
                }
            } else {
                // then lorder[1]==1
                if (OoOv[1] > 0)  {
                    // -- edge 0-1
                    if (OoOv[0] > 0)
                        posclosestvert=0;
                    else
                        posclosestvert=1;
                } else {
                    // -- edge 3-2
                    if (OoOv[0] > 0)
                        posclosestvert=3;
                    else
                        posclosestvert=2;
                }
            }
        } else {
            // -- then face 6
            if (lorder[1] == 0) {
                if (OoOv[0] > 0)  {
                    // -- edge 4-7
                    if (OoOv[1] > 0)
                        posclosestvert=4;
                    else
                        posclosestvert=7;
                } else {
                    // -- edge 5-6
                    if (OoOv[1] > 0)
                        posclosestvert=5;
                    else
                        posclosestvert=6;
                }
            } else {
                // then lorder[1]==1
                if (OoOv[1] > 0)  {
                    // -- edge 4-5
                    if (OoOv[0] > 0)
                        posclosestvert=4;
                    else
                        posclosestvert=5;
                } else {
                    // -- edge 6-7
                    if (OoOv[0] > 0)
                        posclosestvert=7;
                    else
                        posclosestvert=6;
                }
            }
        }
    }
}


void Cdvoxel::UpdateOriginandDensity(const Cvec &vo,const float &d) {
  base.o=vo;
  dens=d;
  UpdateVertLOS();
  obsovox=base.u*(obs.o-base.o);
  obsovoxnorm=obsovox/side;
  CalcClosestVert2();
}


std::ostream& operator << (std::ostream& os,const Cdvoxel& c) {
    os << std::endl
    << "Center : " << c.base.o << std::endl
    << "Matrix : " << c.base.u << std::endl
    << "Side : " << c.side << std::endl
        << "Volume : " << c.volume << std::endl
        << "Density : " << c.dens << std::endl
    << "Obs.o : " << c.obs.o << std::endl
    << "Obs.u : " << c.obs.u << std::endl
    << "Pos Closest Vert : " << c.posclosestvert << std::endl;
    for (int i=0;i<2;i++)
        std::cout << "FACEID["<<i<<"] : " << FACEID[c.posclosestvert][i] << std::endl;
    std::cout << "FACEID[2] : " << FACEID[c.posclosestvert][2];
    return os;
}


const void Cdvoxel::outputVert(){
  for (unsigned int i=0;i<vertobs.size();i++) 
    std::cout << i<< " : " << vertobs[i] << std::endl;
}


Cdvoxel::~Cdvoxel() {
    // TODO: put destructor code here
}


