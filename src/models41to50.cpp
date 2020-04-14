
#include "models41to50.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "Cvec.h"
#include "constant.h"
#include "rtmiscfunc.h"

// -- density 41
// Tube shell model based on density33 but with an extra parameters that takes into account the Ne
// tube shell model
float CModel41::Density(const Cvec &v) {
  float x=fabs(v.v[0]),y=v.v[1],z=v.v[2],r=v.norm();
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
 float nemin;// = 1e4;
 //nemin=density21(v,pparam,temperature);
 float nemax;// = 1e5;

 nemin=pparam[5];
 //nemax=pparam[6];
 nemax=nemin;
 float ratio = pparam[4]; //ratio of tube radius to height //.2;
 
 // -- coeff for the decreasing of the front 
 float wrapcoeff=pparam[6]; //0.2;  // between 0 and 1

 // -- edge slope
 float edgeslope=pparam[7]; //-4; edge smoothing

 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan(alpha)-rs;

 if (testside < 0) {
   // ---- foot side
   //return 0.; // Tres bourrin : no feet !

   // -- dist point - foot
   float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // -- projection point on foot
   Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // -- angle (point to orthoproj on foot) - plane of the structure
   float gammafactor=acos((Cvec(x,y,0.)-vpf).norm() / r);

   // -- check on which side of the foot we are
   if ((y * sin(alpha)) < (x * cos(alpha))) gammafactor = PI - gammafactor;
   gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);

   // -- take into account the skin thickness
   if (r <= v.norm()*ratio && r > v.norm()*ratio-d0) {

     // -- edge factor
     float edgefactor=1.+((v.norm()*ratio)-r)*edgeslope;
     if (edgefactor < 0) edgefactor=0.;

     // -- compute Ne
     if (gammafactor > 0) nel = edgefactor*gammafactor*(nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb)); else nel=0;
   }

   return(nel);

   //if (r < outerDist && r > innerDist) {return 1e4;} else return 0.;

 } else {
   // -- shell side
   // angular position
   float rrr,beta,ttt;
   cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
   // -- dist point - circle
   Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);

   // -- angle (point to orthoproj on foot) - plane of the structure
   float gammafactor=acos((Cvec(x,y,0.) - Cvec(rc*cos(beta),(rs+rc*sin(beta)),0.)).norm() / vpc.norm());

   // -- check if inside or outside the circle
   if ((x*x+(y-rs)*(y-rs)) > (rc*rc)) gammafactor = PI - gammafactor;
   gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);

   
   // -- take into account the skin thickness   
   if (vpc.norm() <= v.norm()*ratio && vpc.norm() > v.norm()*ratio-d0) {
     // -- edge factor
     float edgefactor=1.+((v.norm()*ratio)-vpc.norm())*edgeslope;
     if (edgefactor < 0) edgefactor=0.;

     // -- compute Ne  
     if (gammafactor > 0) nel = edgefactor*gammafactor*(nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb)); else nel=0;
   }

   return(nel);   

   //if (vpc.norm() < outerDist && vpc.norm() > innerDist) {return 1e4;} else return 0.;
 }


 return 0;
} 
void CModel41::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel41::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=2;
	vp.push_back(moddefparam("","Tube shell model based on  model 33.","",""));	
	vp.push_back(moddefparam("d0","0.7","FULL thickness of shell","Rsun"));	
	vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));	
	vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));	
	vp.push_back(moddefparam("rf","10.","dist junction line-circle","Rsun"));	
	vp.push_back(moddefparam("ratio","0.2","ratio of tube radius to height",""));	
	vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));	
	vp.push_back(moddefparam("wrapcoeff","0.2","Front to back decreasing coeff.",""));	
	vp.push_back(moddefparam("edgeslope","-4","edge smoothing",""));	
	
	return;
}



// -- density 42
// Tube shell model based on density33 but with an extra parameter
//      
float CModel42::Density(const Cvec &v) {
  //  float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  // Tres tres bourrin: I've shifted the order of the axis 
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0;
  float d0 = pparam[0]; //0.7; //FULL thickness of shell
  float rb = pparam[1]; //2.55; //dist to bottom of structure
  if (rr <= rb) return nel;

 // -- Oxz is plane of symmetry
 //x=abs(x);
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
 float nemin;// = 1e4;
 //nemin=density21(v,pparam,temperature);
 float nemax;// = 1e5;

 nemin=pparam[5];
 //nemax=pparam[6];
 nemax=nemin;
 float ratio = pparam[4]; //ratio of tube radius to height //.2;
 
 // -- coeff for the decreasing of the front 
 float wrapcoeff=pparam[6]; //0.2;  // between 0 and 1

 // -- edge slope
 float edgeslope=pparam[7]; //-4; edge smoothing

 // -- transition region of thickness at the foot-circle junction
 float ratio1=ratio;
 float ratio2=pparam[8]; // [0.35]
 float stiffness=pparam[9];  // [1.] stiffness of the sigmoidal transition 

 float transit=1./(1.+exp(-stiffness*(rr-rf)));
 ratio=ratio1*(1-transit)+ratio2*transit;

 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan(alpha)-rs;

 if (testside < 0) {
   // ---- foot side
   //return 0.; // Tres bourrin : no feet !

   // -- dist point - foot
   float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // -- take into account the skin thickness
   if (r <= v.norm()*ratio && r > v.norm()*ratio-d0) {

     // -- projection point on foot
     Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     // -- angle (point to orthoproj on foot) - plane of the structure
     float gammafactor=acos((Cvec(x,y,0.)-vpf).norm() / r);
     // -- check on which side of the foot we are
     if ((y * sin(alpha)) < (x * cos(alpha))) gammafactor = PI - gammafactor;
     gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);
     
     // -- edge factor
     float edgefactor=1.+((v.norm()*ratio)-r)*edgeslope;
     if (edgefactor < 0) edgefactor=0.;

     // -- compute Ne
     if (gammafactor > 0) nel = edgefactor*gammafactor*(nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb)); else nel=0;
   }

   return(nel);

   //if (r < outerDist && r > innerDist) {return 1e4;} else return 0.;

 } else {
   // -- shell side
   // angular position
   float rrr,beta,ttt;
   cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
   // -- dist point - circle
   Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);
   
   // -- take into account the skin thickness   
   if (vpc.norm() <= v.norm()*ratio && vpc.norm() > v.norm()*ratio-d0) {

     // -- angle (point to orthoproj on foot) - plane of the structure
     float gammafactor=
       acos((Cvec(x,y,0.) - Cvec(rc*cos(beta),
			       (rs+rc*sin(beta)),
			       0.)).norm() / vpc.norm());
     // -- check if inside or outside the circle
     if ((x*x+(y-rs)*(y-rs)) > (rc*rc)) gammafactor = PI - gammafactor;
     gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);

     // -- edge factor
     float edgefactor=1.+((v.norm()*ratio)-vpc.norm())*edgeslope;
     if (edgefactor < 0) edgefactor=0.;

     // -- compute Ne  
     if (gammafactor > 0) nel = edgefactor*gammafactor*(nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb)); else nel=0;
   }

   return(nel);   

   //if (vpc.norm() < outerDist && vpc.norm() > innerDist) {return 1e4;} else return 0.;
 }


 return 0;
} 
void CModel42::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel42::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=2;
	vp.push_back(moddefparam("","Tube shell model.","",""));	
	vp.push_back(moddefparam("d0","0.7","FULL thickness of shell","Rsun"));	
	vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));	
	vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));	
	vp.push_back(moddefparam("rf","10.","dist junction line-circle","Rsun"));	
	vp.push_back(moddefparam("ratio","0.2","ratio of tube radius to height",""));	
	vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));	
	vp.push_back(moddefparam("wrapcoeff","0.2","Front to back decreasing coeff.",""));	
	vp.push_back(moddefparam("edgeslope","-4","edge smoothing",""));	
	vp.push_back(moddefparam("ratio2","0.2","",""));	
	vp.push_back(moddefparam("stiffness","1.","stiffness of the sigmoidal transition",""));	

	return;
}


// Tube with an helix

float CModel43::Density(const Cvec &v) {
  //  float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  // Tres tres bourrin: I've shifted the order of the axis 
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0;
  float d0 = pparam[0]; //0.7; //FULL thickness of shell
  float rb = pparam[1]; //2.55; //dist to bottom of structure
  if (rr <= rb) return nel;

 // -- Oxz is plane of symmetry
 //x=abs(x);
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
// float rh=rs+rc; // max height of structure???

 float halfdist=rf+(rc*(alpha+PI/2)) ; // half distance on skeleton at half hight 

 //Ne densities at min and max heights (range is from rb to rh)
 float nemin;// = 1e4;
 //nemin=density21(v,pparam,temperature);
 float nemax;// = 1e5;

 nemin=pparam[5];
 //nemax=pparam[6];
 nemax=nemin;
 float ratio = pparam[4]; //ratio of tube radius to height //.2;
 
 // -- coeff for the decreasing of the front 
// float wrapcoeff=pparam[6]; //0.2;  // between 0 and 1

 // -- edge slope
// float edgeslope=pparam[7]; //-4; edge smoothing

 // -- transition region of thickness at the foot-circle junction
 float ratio1=ratio;
 float ratio2=pparam[8]; // [0.35]
 float stiffness=pparam[9];  // [1.] stiffness of the sigmoidal transition 

 float omegah=1./3; // (radian per Rsun)
 float helphase=PI/2.; // phase in radian 
 float helthick=40.*DTOR; // thickness of the helix (radian)
 float r0=3.; // radius of the tube

 float transit=1./(1.+exp(-stiffness*(rr-rf)));
 ratio=ratio1*(1-transit)+ratio2*transit;

 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan(alpha)-rs;

 if (testside < 0) {
   // ---- foot side
   //return 0.; // Tres bourrin : no feet !

   // -- dist point - foot
   float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // -- take into account the skin thickness
   if (r <= r0 && r > r0-d0) {

     // -- projection point on foot
     Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

     // -- angle of the helix depending on the point 
     //    projection on the foot height
     float skeletondist=((v.v[1] < 0.) ? 2*halfdist-vpf.norm() : vpf.norm());
     float thetahelix=(skeletondist*omegah+helphase)-floor((skeletondist*omegah+helphase)/TWOPI)*TWOPI; // modulo two pi

     // -- angle (point to orthoproj on foot) - plane of the structure
     float gammafactor=acos((Cvec(x,y,0.)-vpf).norm() / r);
     // -- check on which side of the foot we are
     if ((y * sin(alpha)) < (x * cos(alpha))) gammafactor = PI - gammafactor;
     if (z < 0) gammafactor=TWOPI - gammafactor;
     //gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);
     
     // -- edge factor
     //float edgefactor=1.;
     //float edgefactor=1.+((v.norm()*ratio)-r)*edgeslope;
     //if (edgefactor < 0) edgefactor=0.;

     // -- compute Ne
     float angdiff=gammafactor-thetahelix;
     if (angdiff < 0) angdiff+=TWOPI;
     float angdiff2=angdiff-PI;
     if (angdiff2 < 0) angdiff2+=TWOPI;
     if ((angdiff <= helthick) || (angdiff2 <= helthick)) nel = nemin; else nel=0;
   }

   return(nel);

   //if (r < outerDist && r > innerDist) {return 1e4;} else return 0.;

 } else {
   //return 0.;
   // -- shell side
   // angular position
   float rrr,beta,ttt;
   cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
   // -- dist point - circle
   Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);

   // -- take into account the skin thickness
   if (vpc.norm() <= r0  && vpc.norm() > r0-d0) {

     // Distance (foot-circle junction) to projected point on circle
     float totdist=rf+rc*(beta+alpha+((beta > PI ) ? -TWOPI : 0 ));
     
     // -- angle of the helix depending on the point 
     //    projection on the foot height
     float skeletondist=((v.v[1] < 0.) ? 2*halfdist-totdist : totdist);
     float thetahelix=(skeletondist*omegah+helphase)-floor((skeletondist*omegah+helphase)/TWOPI)*TWOPI;

     // -- angle (point to orthoproj on circle) - plane of the structure
     // float gammafactor=acos((Cvec(x,y,0.)-vpf).norm() / r);
    float gammafactor=
       acos((Cvec(x,y,0.) - Cvec(rc*cos(beta),
			       (rs+rc*sin(beta)),
			       0.)).norm() / vpc.norm());
     // -- check if inside or outside the circle
     if ((x*x+(y-rs)*(y-rs)) > (rc*rc)) gammafactor = PI - gammafactor;
      if (z < 0) gammafactor=TWOPI - gammafactor;
      //gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);

     // -- edge factor
//     float edgefactor=1.+((v.norm()*ratio)-vpc.norm())*edgeslope;
//     if (edgefactor < 0) edgefactor=0.;

    // -- compute Ne
      float angdiff=gammafactor-thetahelix;
      if (angdiff < 0) angdiff+=TWOPI;
      float angdiff2=angdiff-PI;
      if (angdiff2 < 0) angdiff2+=TWOPI;
      if ((angdiff <= helthick) || (angdiff2 <= helthick)) nel = nemin; else nel=0;

   }

   return(nel);   

   //if (vpc.norm() < outerDist && vpc.norm() > innerDist) {return 1e4;} else return 0.;
 }


 return 0;
} 
void CModel43::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel43::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=2;
	vp.push_back(moddefparam("","Tube shell with an helix.","",""));	
	vp.push_back(moddefparam("d0","0.7","FULL thickness of shell","Rsun"));	
	vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));	
	vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));	
	vp.push_back(moddefparam("rf","10.","dist junction line-circle","Rsun"));	
	vp.push_back(moddefparam("ratio","0.2","ratio of tube radius to height",""));	
	vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));	
	vp.push_back(moddefparam("wrapcoeff","0.2","Front to back decreasing coeff.",""));	
	vp.push_back(moddefparam("edgeslope","-4","edge smoothing",""));	
	vp.push_back(moddefparam("ratio2","0.2","",""));	
	vp.push_back(moddefparam("stiffness","1.","stiffness of the sigmoidal transition",""));	

	return;
}




// -- density 44
// Helicoidal ribbon
float CModel44::Density(const Cvec &v) {
  //  float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  //  !! Tres tres bourrin: I've shifted the order of the axis 
  // -- Oxz is plane of symmetry
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0;
  float rb = pparam[0]; //2.55; //dist to bottom of structure
  if (rr <= rb) return nel;

 // -- deal only with first quadrant
 if (y < 0) return nel;

 float alpha=pparam[1]*DTOR; //30.*DTOR; // angle between axis and foot [Deg]
 float rf=pparam[2]; //10.; // dist junction line-circle [Rsun]

 // -- center of the circle depending on the alpha and rf param
 // rs : dist solar center - shell skeleton circle center
 float rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
 // rc : radius shell skeleton circle
 float rc=rf*tan(alpha);
// float rh=rs+rc; // max height of structure???

 float halfdist=rf+(rc*(alpha+PI/2)) ; // half distance on skeleton at half hight 

 //Ne densities at min and max heights (range is from rb to rh)
 float nemin=pparam[3];// = 1e5;

 float omegah=pparam[4]; //1./3; // (radian per Rsun)
 float helphase=pparam[5]; // PI/2.; // phase in radian 

 float xhalfthick=pparam[6]; //2.5; // half thickness of the ribon, X axis
 float yhalfthick=pparam[7]; //0.3; // half thickness of the ribon, Y axis

 // max extend of the ribbon
 float r0=sqrt(xhalfthick*xhalfthick+yhalfthick*yhalfthick); 

 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan(alpha)-rs;

 if (testside < 0) {
   // ---- foot side
   //return 0.; // Tres bourrin : no feet !

   // -- dist point - foot
   float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // -- take into account the skin thickness
   if (r <= r0) {

     // -- projection point on foot
     Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

     // -- angle of the helix depending on the point 
     //    projection on the foot height
     float skeletondist=((v.v[1] < 0.) ? 2*halfdist-vpf.norm() : vpf.norm());
     float thetahelix=(skeletondist*omegah+helphase); //-floor((skeletondist*omegah+helphase)/TWOPI)*TWOPI; // modulo two pi

     // projection of the point on the axis of the slab
     float xrib=pscal(vpf-Cvec(x,y,z),Cvec(cos(alpha)*cos(thetahelix),-sin(alpha)*cos(thetahelix),-sin(thetahelix)));
     float yrib=pscal(vpf-Cvec(x,y,z),Cvec(cos(alpha)*sin(thetahelix),-sin(alpha)*sin(thetahelix),cos(thetahelix)));

     // -- test if within the slab ribbon range and return the Ne
     if ((fabs(xrib) <= xhalfthick) && (fabs(yrib) <= yhalfthick)) nel = nemin; else nel=0;

   }
   return(nel);

 } else {
   // -- shell side
   // angular position
   float rrr,beta,ttt;
   cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
   // -- dist point - circle
   Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);

   // -- take into account the skin thickness
   if (vpc.norm() <= r0 ) {

     // Distance (foot-circle junction) to projected point on circle
     float totdist=rf+rc*(beta+alpha+((beta > PI ) ? -TWOPI : 0 ));
     
     // -- angle of the helix depending on the point 
     //    projection on the foot height
     float skeletondist=((v.v[1] < 0.) ? 2*halfdist-totdist : totdist);
     float thetahelix=(skeletondist*omegah+helphase);//-floor((skeletondist*omegah+helphase)/TWOPI)*TWOPI;

     // projection of the point on the axis of the slab
     float xrib=pscal(Cvec(rc*cos(beta),rs+rc*sin(beta),0)-Cvec(x,y,z),Cvec(cos(beta)*cos(thetahelix),sin(beta)*cos(thetahelix),-sin(thetahelix)));
     float yrib=pscal(Cvec(rc*cos(beta),rs+rc*sin(beta),0)-Cvec(x,y,z),Cvec(cos(beta)*sin(thetahelix),sin(beta)*sin(thetahelix),cos(thetahelix)));

     // -- test if within the slab ribbon range and return the Ne
     if ((fabs(xrib) <= xhalfthick) && (fabs(yrib) <= yhalfthick)) nel = nemin; else nel=0;

   }
   return(nel);   
 }
 return 0;
} 
void CModel44::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel44::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Twisted or helicoidal ribbon.","",""));	
	vp.push_back(moddefparam("rb","1.5","starting point of the structure","Rsun"));	
	vp.push_back(moddefparam("alpha","0.52","angle between axis and feet","rad"));	
	vp.push_back(moddefparam("rf","10","dist junction line/circle","Rsun"));	
	vp.push_back(moddefparam("nemin","1e5","Ne min","cm^-3"));	
	vp.push_back(moddefparam("omegah","0.33","pulsation of the helix","rad/Rsun"));	
	vp.push_back(moddefparam("helphase","1.5708","phase of the helix","rad"));	
	vp.push_back(moddefparam("xhalfthick","2.5","half thickness of the ribbon on X axis","Rsun"));	
	vp.push_back(moddefparam("yhalfthick","0.3","half thickness of the ribbon on Y axis","Rsun"));	

	return;
}



// -- density 45
// Tube shell model based on density42
//      
float CModel45::Density(const Cvec &v) {
  //  float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  // Tres tres bourrin: I've shifted the order of the axis 
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0; // output density

  if (rr <= rb) return nel;

 // -- Oxz is plane of symmetry
 // -- deal only with first quadrant
 if (y < 0) return nel;

 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan(alpha)-rs;
 float diffinout;
 
 // -- compute cross section radius
 float rcross=(v.norm()-h0)*ratio;

 if (testside < 0) {
   // ---- foot side
   // -- dist point - foot
   float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // -- take into account the skin thickness
   // no density further than 3 sigmas
   if (r < rcross+skinsigmafr*3 && r > rcross-skinsigmain*3) {
     // -- projection point on foot
     Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     // -- angle (point to orthoproj on foot) - plane of the structure
     float gammaang=(Cvec(x,y,0.)-vpf).norm() / r;
	if (gammaang < -1.) gammaang=-1;
	if (gammaang > 1.) gammaang=1;
     float gammafactor=acos(gammaang);
     // -- check on which side of the foot we are
     if ((y * sin(alpha)) < (x * cos(alpha))) gammafactor = PI - gammafactor;
     gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);
     
     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit
     //    also check the side to allow asymmetry
     diffinout=rcross-r;
     nel=nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

     // -- compute Ne
     //if (gammafactor > 0) nel = gammafactor*(nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb)); else nel=0;
     if (gammafactor > 0) nel *= gammafactor; else nel=0;
   }

   return(nel);

   //if (r < outerDist && r > innerDist) {return 1e4;} else return 0.;

 } else {
   // -- shell side
   // angular position
   float rrr,beta,ttt;
   cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
   // -- dist point - circle
   Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);
   
   // -- take into account the skin thickness   
   // no density further than 3 sigmas
   if (vpc.norm() < rcross+skinsigmafr*3 && vpc.norm() > rcross-skinsigmain*3) {

     // -- angle (point to orthoproj on foot) - plane of the structure
	float gammaang=(Cvec(x,y,0.) - Cvec(rc*cos(beta),
			       (rs+rc*sin(beta)),
			       0.)).norm() / vpc.norm();
	if (gammaang < -1.) gammaang=-1;
	if (gammaang > 1.) gammaang=1;
     float gammafactor=acos(gammaang);

     // -- check if inside or outside the circle
     if ((x*x+(y-rs)*(y-rs)) > (rc*rc)) gammafactor = PI - gammafactor;
     gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);

     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit 
     //    also check the side to allow assymetry
     diffinout=rcross-vpc.norm();
     nel=nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

     // -- compute Ne  
     //if (gammafactor > 0) nel = gammafactor*(nemin + (v.norm()-rb)*(nemax-nemin)/(rh-rb)); else nel=0;
     if (gammafactor > 0) nel *= gammafactor; else nel=0;
   }
   return(nel);   
 }
 return 0;
} 
//! Inititialization of the parameters for the Model 45
void CModel45::initParam(float* pparam) {
	rb = pparam[0]; //2.55; //dist to bottom of structure
	alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
	rf=pparam[2]; //10.; // dist junction line-circle [Rsun]
 	ratio = pparam[3]; //ratio of tube radius to height //.2
	nemin=pparam[4];
    // -- coeff for the decreasing of the front 
	wrapcoeff=pparam[5]; //0.2;  // between 0 and 1
	h0=pparam[6]; // offset height aspect ratio: usually = 1
	// -- center of the circle depending on the alpha and rf param
	// rs : dist solar center - shell skeleton circle center
	rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
	// rc : radius shell skeleton circle
	rc=rf*tan(alpha);
	// sigma of the skin gaussian profile (in and front to all assymetry)
	skinsigmain=pparam[8]; // -- inner sigma
	skinsigmafr=pparam[9]; // -- front sigma
}

void CModel45::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Tube shell model, Gaussian skin profile.","",""));	
	vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));	
	vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));	
	vp.push_back(moddefparam("rf","10.","dist junction line-circle","Rsun"));	
	vp.push_back(moddefparam("ratio","0.2","ratio of tube radius to height",""));	
	vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));	
	vp.push_back(moddefparam("wrapcoeff","0.2","Front to back decreasing coeff",""));	
	vp.push_back(moddefparam("h0","1.","Offset height for the aspect ratio: 1. is at the solar limb.","Rsun"));	
	vp.push_back(moddefparam("stiffness","1.","stiffness of the sigmoidal transition: UNUSED",""));	
	vp.push_back(moddefparam("skinsigmain","0.01","inner sigma",""));	
	vp.push_back(moddefparam("skinsigmafr","0.1","front sigma",""));	

	return;
}



// -- density 47
// Tube shell model based on density45.
float CModel47::Density(const Cvec &v) {
  //  float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  // Tres tres bourrin: I've shifted the order of the axis 
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0; // output density
  if (rr <= rb) return nel;
  
 // -- Oxz is plane of symmetry
 // -- deal only with first quadrant
 if (y < 0) return nel;

 float transit=1./(1.+exp(-stiffness*(rr-rf)));
 ratio=ratio1*(1-transit)+ratio2*transit;

 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan(alpha)-rs;
 float diffinout;
 float heightfactor=0;
 float eee=0;
 float neout=0;
 if (testside < 0) {
   // ---- foot side
   // -- dist point - foot
   float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));

   // -- take into account the skin thickness
   // no density further than 3 sigmas
   if (r < v.norm()*ratio+skinsigmafr*3 && r > v.norm()*ratio-skinsigmain*3) {
     // -- projection point on foot
     Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
	 	   
	 // -- linear decreasing of Ne
	eee=vpf.norm()/halfdist-1;
	eee*=eee;
	eee/=ldecdens;
	if (eee < 5e1) heightfactor=exp(-eee); // avoid floating underflow

	// -- angle (point to orthoproj on foot) - plane of the structure
     float gammaang=(Cvec(x,y,0.)-vpf).norm() / r;
	if (gammaang < -1.) gammaang=-1;
	if (gammaang > 1.) gammaang=1;
     float gammafactor=acos(gammaang);
     // -- check on which side of the foot we are
     if ((y * sin(alpha)) < (x * cos(alpha))) gammafactor = PI - gammafactor;
     gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);
     
     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit
     //    also check the side to allow assymetry
     diffinout=v.norm()*ratio-r;
     neout=neinit*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

     // -- compute Ne
     if (gammafactor > 0) nel = gammafactor*neout; else nel=0;
   }

   return(nel*heightfactor);

 } else {
   // -- shell side
   // angular position
   float rrr,beta,ttt;
   cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
   // -- dist point - circle
   Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);
   
   // -- take into account the skin thickness   
   // no density further than 3 sigmas
   if (vpc.norm() < v.norm()*ratio+skinsigmafr*3 && vpc.norm() > v.norm()*ratio-skinsigmain*3) {

	// Distance (foot-circle junction) to projected point on circle
	float totdist=rf+rc*(beta+alpha+((beta > PI ) ? -TWOPI : 0 ));  
	// -- linear decreasing of Ne
	eee=totdist/halfdist-1;
	eee*=eee;
	eee/=ldecdens;
	if (eee < 5e1) heightfactor=exp(-eee); // avoid floating underflow

     // -- angle (point to orthoproj on foot) - plane of the structure
	float gammaang=(Cvec(x,y,0.) - Cvec(rc*cos(beta),
			       (rs+rc*sin(beta)),
			       0.)).norm() / vpc.norm();
	if (gammaang < -1.) gammaang=-1;
	if (gammaang > 1.) gammaang=1;
     float gammafactor=acos(gammaang);

     // -- check if inside or outside the circle
     if ((x*x+(y-rs)*(y-rs)) > (rc*rc)) gammafactor = PI - gammafactor;
     gammafactor=((gammafactor/PI)-wrapcoeff)/(1.-wrapcoeff);

     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit 
     //    also check the side to allow assymetry
     diffinout=v.norm()*ratio-vpc.norm();
     neout=neinit*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

     // -- compute Ne  
     if (gammafactor > 0) nel = gammafactor*neout; else nel=0;
   }
   return(nel*heightfactor);   
 }
 return 0;
} 
// Inititialization of the parameters for the Model 47
void CModel47::initParam(float* pparam)
{
	rb = pparam[0]; //2.55; //dist to bottom of structure
	alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
	rf=pparam[2]; //10.; // dist junction line-circle [Rsun]
 // -- center of the circle depending on the alpha and rf param
 // rs : dist solar center - shell skeleton circle center
	rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
 // rc : radius shell skeleton circle
	rc=rf*tan(alpha);
 
	halfdist=rf+(rc*(alpha+PI/2)) ; // half distance on skeleton at half hight 

	ratio1 = pparam[3]; //ratio of tube radius to height //.2;
	neinit=pparam[4]; // 30000

 // -- coeff for the decreasing of the front 
	wrapcoeff=pparam[5]; //0.2;  // between 0 and 1
 // -- transition region of thickness at the foot-circle junction
	//ratio1=ratio;
	ratio2=pparam[6]; // [0.35]
	stiffness=pparam[7];  // [1.] stiffness of the sigmoidal transition 

 // sigma of the skin gaussian profile (in and front to all assymetry)
	skinsigmain=pparam[8]; // -- inner sigma // 0.2
	skinsigmafr=pparam[9]; // -- front sigma // 0.3
	ldecdens=pparam[10]; // -- sigma^2 decreasing of density with height // 0.04
	
}
void CModel47::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Tube shell model, asymmetric Gaussian skin profile.","",""));	
	vp.push_back(moddefparam("rb","1.5","dist to bottom of structure","Rsun"));	
	vp.push_back(moddefparam("alpha","0.523599","angle between axis and foot","rad"));	
	vp.push_back(moddefparam("rf","5.","dist junction line-circle","Rsun"));	
	vp.push_back(moddefparam("ratio1","0.2","ratio of tube radius to height",""));	
	vp.push_back(moddefparam("neinit","30000","max electron density","e-/cm3"));	
	vp.push_back(moddefparam("wrapcoeff","0.2","transition region of thickness at the foot-circle junction",""));	
	vp.push_back(moddefparam("ratio2","0.35","",""));	
	vp.push_back(moddefparam("stiffness","1.","stiffness of the sigmoidal transition",""));	
	vp.push_back(moddefparam("skinsigmain","0.2","trailing edge sigma",""));	
	vp.push_back(moddefparam("skinsigmafr","0.3","leading edge sigma",""));	
	vp.push_back(moddefparam("ldecdens","0.04","sigma^2 decreasing of density with height",""));	
	
	return;
}







// -- density 49
// Streamer slab model
// see notes from 10 nov 2005
float CModel49::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  if (phi < 1.0472 or phi > 2.0944) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
                                                // see 7 May 2004
  // *2 correction fudge factor (don't know why that works !)
   // -- updated on June 15th using pminimizer02.pro
  float c[5]={0.,1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
//  float cw[4]={17.79/4.71,-85.93/4.71,138.32/4.71,391.45/4.71}; // valid from 2.5 to 20 Rsun
  //float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
    
//	float cw[4]={110.252*DTOR*DTOR,-520.597*DTOR*DTOR,974.667*DTOR*DTOR,0};
float cw[5]={  0., 101.596 *DTOR*DTOR ,   -368.355  *DTOR*DTOR  ,  237.046  *DTOR*DTOR  , 1051.72*DTOR*DTOR };
  // see 10 May 2004
  //float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};

	//float cu[2]={1.65217*DTOR,0.936756*DTOR};
	float cu[5]={4.57523 *DTOR  , -29.5290   *DTOR  ,  106.808   *DTOR  ,  6.36859  *DTOR ,   -407.585*DTOR};
  float coef[5]={0,-1,-2,-3,-4};
  float nel=0;

  // -- streamer half thickness
  float w0=0;
  float u0=0;//cu[0]+cu[1]*r; 
  
  for(int i=0;i<=4;i++) {
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

 
  return nel;
}
void CModel49::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=4;
	vp.push_back(moddefparam("","Streamer slab model.","",""));	
	return;
}






// -- density 50
// Slab model with modulation of the slab FO. 
// see 11 May 2004
// see 14 June 2004
// see 31 Jan 2005
// see 14 Nov 2005
//
float CModel50::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2],r=v.norm();
  if (r <= 2.57) return 0.;
  float rbidon,phi,theta;
  cvcoord(x,y,z,&rbidon,&phi,&theta);
  
  float sechang=pparam[5];

  if (phi < (PI/2-sechang) or phi > (PI/2+sechang)) return 0.;
  //if (y < 0) return 0.;
  
  //  float c[4]={-1168765.3   ,    20505397.,  -1.0056558e+08 ,  1.8379559e+08};
  //float c[4]={7.27e3*2,5.899e5*2,-2.781e6*2,2.932e7*2}; // valid from 2.5 to 20 Rsun
                                                // see 7 May 2004
  // *2 correction fudge factor (don't know why that works !)
 //  float cw[4]={17.79/4.71,-85.93/4.71,138.32/4.71,391.45/4.71}; // valid from 2.5 to 20 Rsun
 float cw[5]={  0., 119.093 *DTOR*DTOR ,   -591.955  *DTOR*DTOR  ,  1121.48  *DTOR*DTOR  , -47.338 * DTOR*DTOR }; // values from Nov 2005
// float cw[4]={17.79*DTOR*DTOR,-85.93*DTOR*DTOR,138.32*DTOR*DTOR,391.45*DTOR*DTOR}; // valid from 2.5 to 20 Rsun
                                            // see 10 May 2004
  //float cu[4]={60.9441*DTOR   ,  -542.687*DTOR  ,    1889.58*DTOR   ,  -2289.55*DTOR};
float cu[5]={4.40467 *DTOR  , -31.8336   *DTOR ,  155.048 *DTOR , -206.623  *DTOR ,   -130.65*DTOR}; // values from Nov 2005

  float coef[5]={0,-1,-2,-3,-4};
  float nel=0;

  // -- streamer half thickness
  float w0=0,u0=0; 
  
  float *pp;
  pp=pparam;

  for(int i=0;i<=4;i++) {
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
  float maxprofpos=pparam[6];

  float aaa=(maxprofpos*(sechang-PI/2))/(2*sechang);
  float bbb=maxprofpos/(2*sechang);
  
  int pos=(int) (aaa+bbb*phi);

  if (pos < 0) pos=0;
  if (pos > maxprofpos) pos=(int) maxprofpos;

  nel*=*(pparam+pos+7); // -- the profile starts at the pos 7 of pparam

  return nel;
}
void CModel50::initParam(float *pparam)
{
this->pparam=pparam;
}
void CModel50::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=3;
	vp.push_back(moddefparam("","Streamer slab model, parameter passing for minimization.","",""));	
	return;
}


