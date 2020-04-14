
#include "models31to40.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "Cvec.h"
#include "constant.h"
#include "rtmiscfunc.h"

//-----------------------------------------------
// -- density  31
// Graduated Curved Cylindrical Shell Model
float CModel31::Density(const Cvec &v) {
  
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  
  //COMMON /params/rh,dh,d0,rc,rb,db,nemin,nemax
  float nel = 0;
  float nemax = pparam[0]; //1.5;		//non-constant density
  float nemin = pparam[1]; //1;
  //int makecloud = 0;

  // Input Parameters:
  float rh = pparam[2]; //14;      // height of leading edge of cylinder at the nose
  float dh  = pparam[3]; //4;     // diameter of cylinder at the top
  float d0 = pparam[4]; //0.3;   // 1/2 thickness of shell
  float rc = pparam[5]; //3;    // radius of curvature of the cylinder
  float rb = pparam[6]; //2;   // height of base of cylinder
  float db = pparam[7]; //.1; // diameter of cylinder at the base
  
  float del_d = (dh-db)/(rh-rb);
  float r0 = rh-dh*0.5-rc;             // height of center of curvature of flux rope
  float r1 = sqrt(r0*r0-(rc+dh*.5)*(rc+dh*.5));   // height of transition of conical legs to curved cylinder
  float theta = atan2(rc,r1);         // half angle of width of curved cylinder
//  float r2 = r1*cos(theta);         // projected height of transition

//  float temp = rc+d0*.5;                               
//  float thetamax = atan2(rc+d0,r1);// max angle including shell
//  float th = atan2(x,y);	 // angle of point relative to axis

  float dely1 = 0;	      // distance of the point to shell middle on far side
  float dely2 = 0;	     // distance of the point to shell middle on near side
  float delx1 = 0;	    // same for x
  float delx2 = 0;         //
  //float rdist = 0;        //
  
  ;if (y < rb || y > rh)  return(nel);	      // beyond the model
  if (r < r1) {
    // in cone region (legs)

    float rmid = y*tan(theta);	   // distance of middle of shell to central axis
    if (x < 0) rmid = -rmid;
    float rcirc = 0.5*(db+(r-rb)*del_d);      // radius of cylinder cross section
    float rdist = sqrt((x-rmid)*(x-rmid)+z*z);
    float delr = fabs(rdist-rcirc);
    if (delr <= d0) {
      nel=nemin + (r-rb)*(nemax-nemin)/(rh-rb);
    }
    // if (makecloud){
//       if (nel > 0) {
// 	ofstream myFile("output31.txt",ios::app);
// 	myFile<<x<<" "<<y<<" "<<z<<" "<<nel<<endl;  
// 	myFile.close();
//       }
//     }
    return (nel);

    dely1 = fabs(rmid+rcirc-x);	 // distance of the point to shell middle on far side
    dely2 = fabs(rmid-rcirc-x);	// distance of the point to shell middle on near side
    delx1 = fabs(z-rcirc);      // same for x
    delx2 = fabs(z+rcirc);     // same for x

  } else {		    // in hemispherical region

    float rdist  = sqrt (x*x+(y-r0)*(y-r0)); // distance from center of curvature
    float phi = atan2(x,y-r1);
    float ht = rc*cos(phi)+r1; //new (bad): (rc+(d0/2)*cos(phi)+r0)
    float rcirc = 0.5*(db+(ht-rb)*del_d);  // radius of cylinder cross section
    float xc = rc*sin(phi);
    rdist = sqrt((x-xc)*(x-xc)+(y-ht)*(y-ht)+z*z); // distance of the point to radius of circle
    float delr = fabs(rdist-rcirc);
    //dely2 = abs(rc-rcirc-rdist);  // distance of the point to shell middle on near side
    //delx1 = abs(z-rcirc);        // same for x
    //delx2 = abs(z+rcirc);       // same for x
    if (delr <= d0)  nel=nemin + (r-rb)*(nemax-nemin)/(rh-rb);

    // if (makecloud){
//       if (nel > 0) {
// 	ofstream myFile("output31.txt",ios::app);
// 	myFile<<x<<" "<<y<<" "<<z<<" "<<nel<<endl;  
// 	myFile.close();
//       }
//     }
    return (nel);
  }
//  int inshell = 0;
//  float delx = 0;
//  float dely = 0;
  
  if ((delx1 <= d0 || delx2 <= d0) && (dely1 <= d0 || dely2 <= d0)) 
    nel = nemin + (nemax-nemin)*(y-rb)/(rh-rb);

  std::cout<<"SOMETHING WENT WRONG, WE SHOULDN'T BE HERE"<<std::endl;
  // if (makecloud){
//     if (nel > 0) {
//       ofstream myFile("output31.txt",ios::app);
//       myFile<<x<<" "<<y<<" "<<z<<" "<<nel<<endl;  
//       myFile.close();
//     }
//   }
  return(nel);
}
void CModel31::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel31::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Graduated Curved Cylindrical Shell Model","",""));	

	vp.push_back(moddefparam("nemax","1.5","Ne max","cm^-3"));	
	vp.push_back(moddefparam("nemin","1.","Ne min","cm^-3"));	
	vp.push_back(moddefparam("rh","14.","Height of leading edge of cylinder at the nose","Rsun"));	
	vp.push_back(moddefparam("dh","4.","Diameter of cylinder at the top","Rsun"));	
	vp.push_back(moddefparam("d0","0.3","Half thickness of shell","Rsun"));	
	vp.push_back(moddefparam("rc","3.","Radius of curvature of the cylinder","Rsun"));	
	vp.push_back(moddefparam("rb","2.","Height of base of cylinder","Rsun"));	
	vp.push_back(moddefparam("db","0.1","Diameter of cylinder at the base","Rsun"));	
	
	return;
}



// -- density 32
// Streamer belt simulation with source surface field map:
float CModel32::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.55) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- param for the neutral sheet map
  int sang=(int) *pparam,slat=(int) *(pparam+1);

  // -- streamer thickness parameters
  //    w0=d1, u0=d2
  float w0=*(pparam+2),u0=*(pparam+3);

  // ---- get the distance from the point
  //      to the nearest neighbors
  // -- nearest neighbor seeking range
  int srang=15,srlat=15; // in pix

  // -- point to the neutral sheet map
  float *pnsheetmap;
  pnsheetmap=pparam+4;

  float thetannpos,phinnpos,dist,val;

  int nnok=wherennsmoothed(pnsheetmap,sang,slat,srang,srlat,(phi*RADEG),(theta*RADEG),&phinnpos,&thetannpos,&dist,&val);

  if (nnok == 0) return 0.;

  // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
  //float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  //float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};
  float coef[4]={-1,-2,-3,-4};
  float nel=0;
  // -- streamer half thickness
  //float w0=0,u0=0; 

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    nel+=c[i]*powr;
    //w0+=cw[i]*powr;
    //u0+=cu[i]*powr;
    //nel+=c[i]*pow(r,coef[i]);
 }

  float dcross=w0/(2*u0);

  // -- theta is given by the distance of the requested 
  //    point to the neutral line
  theta=dist*DTOR;
  
  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }
  
  return nel;

}
void CModel32::initParam(float *pparam)
{
this->pparam=pparam;
}

// -- density 33
// Tube shell model
float CModel33::Density(const Cvec &v) {
 float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
 float nel = 0;
 float d0 = pparam[0]; //0.7; //FULL thickness of shell
 float rb = pparam[1]; //2.55; //dist to bottom of structure
 if (r <= rb) return nel;
 
 float rbidon,phi,theta;
 cvcoord(x,y,z,&rbidon,&phi,&theta);

 // -- Oyz is plane of symmetry
 x=fabs(x);
 // -- deal only with first quadrant
 if (y < 0) return nel;

 float alpha=pparam[2]*DTOR; //30.*DTOR; // angle between axis and foot [rad]
 float rf=pparam[3]; //10.; // dist junction line-circle [Rsun]

 //area between the following two distances from skeleton will have non-zero Ne density
 //float innerDist = pparam[X]; //1.1 //distance from skeleton to inner radius of shell
 //float outerDist = pparam[X]; //1.5 //distance from skeleton to outer radius of shell

 // -- center of the circle depending on the alpha and rf param
 // rs : dist solar center - shell skeleton circle center
 float rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
 // rc : radius shell skeleton circle
 float rc=rf*tan(alpha);
 float rh = rs+rc; //max height of structure???

 //Ne densities at min and max heights (range is from rb to rh)
 float nemin = 1e4;
 //nemin=density21(v,pparam,temperature);
 float nemax = 1e5;
 nemax=nemin;
 float ratio = pparam[4]; //ratio of tube radius to height //.2;
 
 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan(alpha)-rs;
 
 if (testside < 0) {
   // -- foot side
   //return 0.; // Tres bourrin : no feet !

   // dist point - foot
   float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // constant density for the moment

   //if (r < d0) return(1e4); 
   //else return(0);

   if (r < v.norm()*ratio+d0 && r > v.norm()*ratio) {
     nel = nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb);
   }

   return(nel);

   //if (r < outerDist && r > innerDist) {return 1e4;} else return 0.;

 } else {
   // -- shell side
   // angular position
   float rrr,beta,ttt;
   cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
   // dist point - circle
   Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);
   // constant density for the moment
   
   //if (vpc.norm() < ratio) return(1e4);
   //else return(0);

   if (vpc.norm() < v.norm()*ratio+d0 && vpc.norm() > v.norm()*ratio) {
     nel = nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb);
   }

   return(nel);   

   //if (vpc.norm() < outerDist && vpc.norm() > innerDist) {return 1e4;} else return 0.;
 }


 return 0;
} 
void CModel33::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel33::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Tube shell model.","",""));	

	vp.push_back(moddefparam("d0","0.7","FULL thickness of shell.","Rsun"));	
	vp.push_back(moddefparam("rb","2.55","Dist to bottom of structure.","Rsun"));	
	vp.push_back(moddefparam("alpha","30.","Angle between axis and foot.","deg"));	
	vp.push_back(moddefparam("rf","10","Dist junction line-circle.","Rsun"));	
	vp.push_back(moddefparam("ratio","0.2","ratio of tube radius to height","Rsun"));	
	
	return;
}


// -- density 34
// Test of a simple density rope 
float CModel34::Density(const Cvec &v) {

  //float nehomo=density20(v,pparam,temperature);

  //float returnne=nehomo;
  float returnne=0;

  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.01) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- Oyz is plane of symmetry 
  x=fabs(x);
  // -- deal only with first quadran
  if (y < 0) return returnne;

  float alpha=pparam[0]; //30.*DTOR; // angle between axis and foot [rad]
  float rf=pparam[1]; //3.; // dist junction line-circle [Rsun]
  float tubethick=pparam[2]; //.5; // thickness radius of the rope [Rsun]
  float constneout=pparam[3]; //1e7; // output constant density
  float skinthick=pparam[4]; //0.1 // thickness of the skin


  // -- center of the circle depending on the alpha and rf param
  // rs : dist solar center - shell skeleton circle center
  float rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
  // rc : radius shell skeleton circle
  float rc=rf*tan(alpha);

  // -- on which side are we ? foot or shell ?
  float testside=y+x*tan(alpha)-rs;
  if (testside < 0) {
    // -- foot side
    // dist point - foot
    float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

    // constant density for the moment
    if (r < tubethick && r >= skinthick ) return constneout+returnne; else return returnne;

  } else {
    // -- shell side
    // angular position
    float rrr,beta,ttt;
    cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
    
    // dist point - circle
    Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);
    
    // constant density for the moment
    if (vpc.norm() < tubethick && vpc.norm() >= skinthick) return constneout+returnne; else return returnne;
  }
  return returnne;
}
void CModel34::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel34::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Simple density rope.","",""));	
	
	vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));	
	vp.push_back(moddefparam("rf","3","dist junction line-circle","Rsun"));	
	vp.push_back(moddefparam("tubethick","0.5","thickness radius of the rope","Rsun"));	
	vp.push_back(moddefparam("constneout","1e7","Constant density","cm^-3"));	
	vp.push_back(moddefparam("skinthick","0.1","thickness of the skin","Rsun"));	
	
	return;
}


// -- density 35
// Density and temperature cubes from S.Gibson and B.Low CME progam model.
// 
float CModel35::Density(const Cvec &v,float &temperature) {

  float r,phi,theta;
  cvcoord(v.v[0],v.v[1],v.v[2],&r,&phi,&theta);

  float NeighborDensCube[2][2][2];
  float NeighborTempCube[2][2][2];

  // ---- get the model parameters and point cubes first element
  float *rco,*phico,*thetaco,*denscube,*tempcube;

  rco=&pparam[0];
  phico=&pparam[3];
  thetaco=&pparam[6];
  denscube=&pparam[9];
  tempcube=&pparam[9+((int) rco[2])*((int) phico[2])*((int) thetaco[2])];


  // ---- find the nearest neighbors for tri linear interpolation
  int mr=(int)(-rco[0]/ rco[1] + r / rco[1] + 1);
  if (mr >= rco[2] || mr <= 0) return 0.;

  int mtheta12[2];
  int mtheta=(int)(- thetaco[0]/ thetaco[1] + theta / thetaco[1] + 1);
  if (mtheta >= thetaco[2]) {
    mtheta12[0]=(int) thetaco[2]-1;
    mtheta12[1]=mtheta12[0];
  } else if (mtheta <= 0) {
    mtheta12[0]=0;
    mtheta12[1]=0;
  } else {
    mtheta12[0]=mtheta-1;
    mtheta12[1]=mtheta;   
  }
  int mphi12[2];
  int mphi=(int)(- phico[0]/ phico[1] + phi / phico[1] + 1);
  if (mphi >= phico[2]) {
    mphi12[0]=(int) phico[2]-1;
    mphi12[1]=0;
  } else if (mphi <= 0) {
    mphi12[0]=0;
    mphi12[1]=1;
  } else {
    mphi12[0]=mphi-1;
    mphi12[1]=mphi;
  }


  // ---- fill in the nearest neighbor cube for interpolation
  mr--;
  mtheta--;
  mphi--;

  int xxx,yyy,zzz;
  for (int ii=0;ii<2;ii++) {
    xxx=(mr+ii);
    for (int jj=0;jj<2;jj++) {
      yyy=mphi12[jj]*((int)(rco[2]));
      for (int kk=0;kk<2;kk++) {
	zzz=mtheta12[kk]*((int)(rco[2]) * ((int)(phico[2])));
	NeighborDensCube[ii][jj][kk]=*(denscube+(xxx+yyy+zzz));
	NeighborTempCube[ii][jj][kk]=*(tempcube+(xxx+yyy+zzz));
      }
    }
  }

  // ---- compute intervals for interpolation
  float rint[2]={rco[0] + rco[1] * mr , rco[0] + rco[1] * (mr+1)};
  float phiint[2]={phico[0] + phico[1] * mphi , phico[0] + phico[1] * (mphi+1)};
  float thetaint[2]={thetaco[0] + thetaco[1]*mtheta,
		     thetaco[0] + thetaco[1]*(mtheta+1)};

  // ---- perform trilinear interpolation
  float neinterp=trilininterp(r,phi,theta,rint,phiint,thetaint,NeighborDensCube);
  temperature=trilininterp(r,phi,theta,rint,phiint,thetaint,NeighborTempCube);

  return neinterp;

}
void CModel35::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel35::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Density and temperature cubes from S.Gibson and B.Low CME progam model.","",""));	
	
	vp.push_back(moddefparam("filedens","''","Density file from gibsoncmewrapper.pro",""));	
	vp.push_back(moddefparam("filetemp","''","Temperature file from gibsoncmewrapper.pro",""));	
	vp.push_back(moddefparam("modparam","","",""));	
	vp.push_back(moddefparam("","load_gibsondens,$filedens,$filetemp,modparam","Load the density and temperature from S.Gibson's model.",""));	

	return;
}

// -- density 36
// Streamer belt simulation with source surface field map:
float CModel36::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  
  if (r <= 1.05) return 0.;

  // -- now deal with orthoradial shape
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- param for the neutral sheet map
  int sang=(int) pparam[0],slat=(int) pparam[1];

  // ---- get the distance from the point
  //      to the nearest neighbors
  // -- nearest neighbor seeking range
  int srang=30,srlat=30; // in pix

  // -- point to the neutral sheet map
  float *pnsheetmap;
  pnsheetmap=pparam+2;

  float thetannpos,phinnpos,val;//dist;
  float thetamg;
  float nel;

  int nnok=wherennsmoothed(pnsheetmap,sang,slat,srang,srlat,(phi*RADEG),(theta*RADEG),&phinnpos,&thetannpos,&thetamg,&val);

  if (nnok == 0) thetamg=DTOR*40; 
  else thetamg*=DTOR;
  // if (nnok == 0) thetamg=DTOR*((float) srang * 360 / (float) sang); 
//   else thetamg*=DTOR;
  
  float cp[3]={0.14e7,8.02e7,8.12e7};
  float dp[3]={2.8,8.45,16.87};
  float ccs[3]={1.41e7,16.42e7,61.90e7};
  float dcs[3]={3.39,5.14,12.67};
  float gam[3]={16.3,10.0,43.20};
  float delt[3]={0.5,7.31,7.52};
  
  float np=0,ncs=0,w=0;
  
  for(int i=0;i<=2;i++) {
    np+=cp[i]*pow(r,-dp[i]);
    ncs+=ccs[i]*pow(r,-dcs[i]);
    w+=gam[i]*pow(r,-delt[i]);
  }  
  
  float expo=pow(thetamg/(w*DTOR),2);
  if (expo > 1E1) nel=np; else nel=np+(ncs-np)*exp(-expo);
  
  return nel;
}
void CModel36::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel36::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=2;
	vp.push_back(moddefparam("","Streamer belt simulation with source surface field map, user resolution, Guhata. orthoradial.","",""));	
	
	vp.push_back(moddefparam("crot","1912L","Carrington rotation number",""));	
	vp.push_back(moddefparam("modparam","","",""));	
	vp.push_back(moddefparam("","rdtxtmagmap,nsheetmap,crot=$crot","Load PFSS map from WSO web site.",""));	
	vp.push_back(moddefparam("","modparam=reform(360,181,nsheetmap,n_elements(nsheetmap))","Format the array of parameters.",""));	
	
	return;
}



// -- density 37
// Streamer belt simulation with source surface field map:
float CModel37::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.55) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- param for the neutral sheet map
  int sang=(int) pparam[1],slat=(int) pparam[2];

  // -- point to the neutral sheet map
  float *pnsheetmap;
  pnsheetmap=pparam+3;

  float dist,val;

  // ---- get the distance from the point
  //      here it's given by the magnetic field
  int nnok=getPosOnSSMap(pnsheetmap,sang,slat,(phi*RADEG),(theta*RADEG),dist,val);


  if (nnok == 0) return 0.;

  // -- updated on June 15th using pminimizer02.pro
  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};
  float coef[4]={-1,-2,-3,-4};
  float nel=0;
  // -- streamer half thickness
  float w0=0,u0=0; 

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    nel+=c[i]*powr;
    w0+=cw[i]*powr;
    u0+=cu[i]*powr;
 }

  float dcross=w0/(2*u0);

  // ---- multiply the mag field by a coeff
  theta=dist*DTOR*pparam[0]; // ca c'est tres bourrin !
  
  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }

  return nel;
}
void CModel37::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel37::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=2;
	vp.push_back(moddefparam("","Streamer belt simulation with source surface field map, distance from neutral line is given by the magnetic field value.","",""));	
	
	vp.push_back(moddefparam("sthick","1.","streamer belt thickness coefficient",""));	
	vp.push_back(moddefparam("crot","1912L","Carrington rotation number",""));	
	vp.push_back(moddefparam("modparam","","",""));	
	vp.push_back(moddefparam("","rdtxtmagmap,nsheetmap,crot=$crot","Load PFSS map from WSO web site.",""));	
	vp.push_back(moddefparam("","modparam=reform($sthick,360,181,nsheetmap,n_elements(nsheetmap))","Format the array of parameters.",""));	
	
	return;
}


// -- density 38
// Slab model
// see 11 May 2004
// see 14 June 2004
// see 31 Jan 2005
//
float CModel38::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  float sechang=pparam[4];

  if (phi < (PI/2-sechang) or phi > (PI/2+sechang)) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
                                                // see 7 May 2004
  // *2 correction fudge factor (don't know why that works !)
 //  float cw[4]={17.79/4.71,-85.93/4.71,138.32/4.71,391.45/4.71}; // valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};

  float coef[4]={-1,-2,-3,-4};
  float nel=0;

  // -- streamer half thickness
  float w0=0,u0=0; 
  
  float *pp;
  pp=pparam;

  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    //    nel+=c[i]*powr;
    nel+=(*(pp++))*powr;
    w0+=cw[i]*powr;
    u0+=cu[i]*powr;
 }

  float dcross=w0/(2*u0);

  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }

  // ---- FO modulation profile
  // the *pmodul should contain 120 samples 
  // representing the angles from 60 to 120 deg.,
  // with a range dynamic from 0 to 1
  // 
  float maxprofpos=pparam[5];

  float aaa=(maxprofpos*(sechang-PI/2))/(2*sechang);
  float bbb=maxprofpos/(2*sechang);
  
  int pos=(int) (aaa+bbb*phi);

  if (pos < 0) pos=0;
  if (pos > maxprofpos) pos=(int) maxprofpos;

  nel*=*(pparam+pos+6); // -- the profile starts at the pos 6 of pparam

  return nel;
}
void CModel38::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel38::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
// - pparam[0..3] : coeff of the Ne polynomial
// - pparam[4] : slab half angle
// - pparam[5] : subscript max of the modulation profile (size - 1)
// - pparam[6] : first element of the modulation profile
	flagcase=0;
	vp.push_back(moddefparam("","Streamer slab model with parameters useful for fitting.","",""));	
	vp.push_back(moddefparam("c","[1.34e4,1.15e6,-6.022e6,5.577e7]","Coeff of the Ne polynomial.",""));	
	vp.push_back(moddefparam("alpha","0.52","slab half angle","rad"));	
	vp.push_back(moddefparam("smodprof","119","subscript max of the modulation profile (size - 1)",""));	
	vp.push_back(moddefparam("pm","","",""));	
	vp.push_back(moddefparam("","pm=[$c,$alpha,$smodprof,(cos(findgen($smodprof +1)/2)+1.1)/2.1]","",""));	
	return;
}




// -- density 39
// Test of a simple helix
float CModel39::Density(const Cvec &v) {

  // pos projection point on x axis
  Cvec op=orthoproj(Cvec(0,0,0),Cvec(1,0,0),v);

  // distance projection - sun center
  float dpc=op.norm(); 

  // point of the helix
  float aa=2,bb=1.;
  Cvec ph=Cvec(dpc*bb,aa*cos(dpc*0.8),aa*sin(dpc*0.8));

  // dist helix center - point
  Cvec vr=ph-v;
 
  if (vr.norm() < 0.5) return 1e7; else return 0;




  //float nehomo=density20(v,pparam,temperature);

  //float returnne=nehomo;
  float returnne=0;

  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 1.01) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- Oyz is plane of symmetry 
  x=fabs(x);
  // -- deal only with first quadran
  if (y < 0) return returnne;

  float alpha=pparam[0]; //30.*DTOR; // angle between axis and foot [rad]
  float rf=pparam[1]; //3.; // dist junction line-circle [Rsun]
  float tubethick=pparam[2]; //.5; // thickness radius of the rope [Rsun]
  float constneout=pparam[3]; //1e7; // output constant density


  // -- center of the circle depending on the alpha and rf param
  // rs : dist solar center - shell skeleton circle center
  float rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
  // rc : radius shell skeleton circle
  float rc=rf*tan(alpha);

  // -- on which side are we ? foot or shell ?
  float testside=y+x*tan(alpha)-rs;
  if (testside < 0) {
    // -- foot side

    // point on the foot
    Cvec op=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

    // helix position
    float aa=0.2,omega=0.1; // param of the helix
    Cvec hp=Cvec(aa*cos(op.norm()*omega),aa*sin(op.norm()*omega),0)+op;

    // dist point - helix
    Cvec vr=Cvec(x,y,z)-hp;

    // constant density for the moment
    if (vr.norm() < tubethick) return constneout+returnne; else return returnne;

  } else {
    // -- shell side
    // angular position
    float rrr,beta,ttt;
    cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
    
    // dist point - circle
    Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);
    
    // constant density for the moment
    if (vpc.norm() < tubethick) return constneout+returnne; else return returnne;
  }
  return returnne;
}
void CModel39::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel39::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Test of a simple helix.","",""));	

	vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));	
	vp.push_back(moddefparam("rf","3.","dist junction line-circle","Rsun"));	
	vp.push_back(moddefparam("tubethick","0.5","thickness radius of the rope","Rsun"));	
	vp.push_back(moddefparam("constneout","1e7","output constant density","cm^-3"));	
	
	return;
}





// -- density 40
// Streamer belt simulation with source surface field map
//
float CModel40::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.55) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);

  // -- param for the neutral sheet map
  int sang=(int) pparam[0],slat=(int) pparam[1];

  // -- point to the neutral sheet map
  float *pnsheetmap;
  pnsheetmap=pparam+3;

  float dist,val;

  // ---- get the distance from the point
  //      here it's given by the magnetic field
  int nnok=getPosOnSSMap(pnsheetmap,sang,slat,(phi*RADEG),(theta*RADEG),dist,val);

  if (nnok == 0) return 0.;
  
  // ---- get the height of the passed Ne layer
//  float layerheigth=pparam[2];

  // ---- compute Saito model (equatorial)
//  float c1=1.36e6;
//  float d1=2.14;
//  float c2=1.68e8;
//  float d2=6.13;

  // ---- compute ratio for scaling
  //if (val < 7e5 ) return 0.;
  //float Neratio=val/(c1*pow(layerheigth,-d1)+c2*pow(layerheigth,-d2));

  // ---- compute Ne and scale
  //float nel=Neratio*(c1*pow(r,-d1)+c2*pow(r,-d2));
  //float nel=(c1*pow(r,-d1)+c2*pow(r,-d2));


  float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
  float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};
  
  float coef[4]={-1,-2,-3,-4};
  float nel=0;
  
  // -- streamer half thickness
  float w0=0,u0=0; 
  
  for(int i=0;i<=3;i++) {
    float powr=pow(r,coef[i]);
    nel+=c[i]*powr;
    w0+=cw[i]*powr;
    u0+=cu[i]*powr;
 }

  float dcross=w0/(2*u0);

  if (fabs(theta) < dcross) {

    float ee=(theta*theta)/w0;
    if (ee > 1E1) return 0.;
    
    nel*=exp(-ee);  
  } else {

    float ee=fabs(theta/u0)-dcross/u0/2;
    if (ee > 1E1) return 0.;

    nel*=exp(-ee);

  }


  // ---- modulate by the layer value
  return val*nel;
}
void CModel40::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel40::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	//! - pparam[0] : longitude size in pix
	//! - pparam[1] : latitude size in pix
	//! - pparam[2] : height of the layer in Rsun
	//! - pparam[3] : Ne layer map (lon,lat)
	vp.push_back(moddefparam("","Streamer belt simulation with source surface field map, Saito radial scaled on the Ne layer map passed in parameter.","",""));	

	vp.push_back(moddefparam("pfssmapfile","'WSOpsffR250CR2012.fts'","!!! ACTUALLY THE MAP SHOULD BE NE NOT A MAG FIELD !!! THIS IS NOT A VALID EXAMPLE !!!",""));	
	vp.push_back(moddefparam("heightpfss","2.5","",""));	
	vp.push_back(moddefparam("mp","","",""));	
	vp.push_back(moddefparam("","if not file_exist($pfssmapfile) then begin","",""));	
	vp.push_back(moddefparam("","rtinitenv","",""));	
	vp.push_back(moddefparam("","progpath=getenv('RT_PATH')+get_delim()","",""));	
	vp.push_back(moddefparam("","endif else progpath=''","",""));	
	vp.push_back(moddefparam("","map=readfits(progpath + $pfssmapfile)","Load the pfss map",""));	
	vp.push_back(moddefparam("","mp=[(size(map,/dim))[0],(size(map,/dim))[1],$heightpfss,reform(map,n_elements(map))]","",""));	

	return;
}

