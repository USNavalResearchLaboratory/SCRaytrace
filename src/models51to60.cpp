
/** \file models51to60.cpp
 * \brief Model 51 to 60
 */


#include "models51to60.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "Cvec.h"
#include "constant.h"
#include "rtmiscfunc.h"

// -- density 51
/*
** cme_demi_sphere - Compute octree for CME half sphere model.
**
** Implemented for comparison between Marseille and Arnaud's renderer
**
**     Computes CME electron density at (x,y,z) with respect to
**    jet axis (Ox)
**
**
**      Created on 14/09/2005 by M. Burtin and F. Saez
*************************************************************************/
float CModel51::Density(const Cvec &v) {
  float x=v.v[0],y=v.v[1],z=v.v[2];
//	_PARAMS *
//cme_demi_sphere (const double x,    /* get jetval at x,y,z */
//          const double y,
//          const double z,
//          const MODELE_PARAMS * jet,
//          PT_PARAMS * result, const LIM_PARAMS * lim)

  double rprime, rprime2, d, jetval;
  //const float pi = 3.14159265358979;

	float alpha,r,dr;
	alpha=pparam[0];
	r=pparam[1];
	dr=pparam[2];
	
  /* square of distance of (x,y,z) from Sun center
   */
  rprime2 = x*x + y*y + z*z;
  rprime2 = sqrt (rprime2);

  /* compute radius of the sphere */

  d = r * tan(alpha*(PI/180)/2.);

  /* if (x,y,z) is inside Sun
   *         then no cme value
   */

  jetval = 0.0;

  if (rprime2 >= 1.5)
  {
      /* definition of half space */

    if (x >= r)
    {    
      /* distance  from CME center
       */
      rprime = sqrt ((r-x)*(r-x) + y*y + z*z);
       
      if (fabs(rprime-d) <= dr/2.)
      {
          jetval = 1;
      }
        }            
  }                

  //result->computed = TRUE;
  //result->elec_density = jetval;
  //result->sun_distance = rprime2;
  //result->neutral_eq = 0;
  return jetval;
}
void CModel51::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	

	flagcase=0;
	vp.push_back(moddefparam("","Half shell from F.Saez","",""));	
	vp.push_back(moddefparam("alpha","30.","Cone angle","deg"));	
	vp.push_back(moddefparam("rh","5.","Height of the hemisphere","Rsun"));	
	vp.push_back(moddefparam("dr","0.5","Thickness of the shell","Rsun"));	

	return;
}



// -- density 52
// Tube shell model based on density47
// With decreasing of density depending on height.
float CModel52::Density(const Cvec &v) {
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
// Inititialization of the parameters for the Model 52
void CModel52::initParam(float* pparam)
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
void CModel52::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Tube shell model based on model 47.","",""));	
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



// -- density 53
// Tube shell model based on density45
//      
float CModel53::Density(const Cvec &v) {
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
  
  if (testside < 0) {
    // ---- foot side
    // -- dist point - foot
    float r=psldist(Cvec(0,0,0),cva,Cvec(x,y,z));

    // -- projection point on foot
    Cvec vpf=orthoproj(Cvec(0,0,0),cva,Cvec(x,y,z));

    // -- compute cross section radius
    float rcross=(vpf.norm()-h0)*ratio;
    // -- take into account the skin thickness
    //    no density further than 3 sigmas
    if (r < (rcross + skinsigmafr * 3) && r > (rcross - skinsigmain * 3)) {
           
      // -- gaussian profile: note it's not exactly the gaussian 
      //    probability function but it speeds up the code a tiny bit
      //    also check the side to allow asymmetry
      diffinout=rcross-r;
      nel=nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

    }
    
    return(nel);
        
  } else {
    // -- shell side
    // angular position
    float rrr,beta,ttt;
    cvcoord(x,y-rs,0,&rrr,&beta,&ttt);
//     cvcoord(v.v[0],v.v[1]-rs,0,&rrr,&beta,&ttt);
   // -- dist point - circle
    Cvec vpc(x-rc*cos(beta),y-(rs+rc*sin(beta)),z);
//  	  Cvec vpc(v.v[0]-rc*cos(beta),v.v[1]-(rs+rc*sin(beta)),v.v[2]);
    // -- projected height on circle skeleton
    Cvec vsk(rc*cos(beta),rs+rc*sin(beta),0);
    
    // -- compute cross section radius
    float rcross=(vsk.norm()-h0)*ratio;
    
    // -- take into account the skin thickness   
    //    no density further than 3 sigmas
	  float vpcnorm=vpc.norm();
    if (vpcnorm < (rcross+skinsigmafr*3) && vpcnorm > (rcross-skinsigmain*3)) {
            
      // -- gaussian profile: note it's not exactly the gaussian 
      //    probability function but it speeds up the code a tiny bit 
      //    also check the side to allow asymmetry
      diffinout=rcross-vpcnorm;
      nel=nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));
      
    }
    return(nel);   
  }
  return 0;
} 

void CModel53::initParam(float* pparam) {
  rb = pparam[0]; //2.55; //dist to bottom of structure
  alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
  rf=pparam[2]; //10.; // dist junction line-circle [Rsun]
  ratio = pparam[3]; //ratio of tube radius to height //.2
  nemin=pparam[4];
  // -- coeff for the decreasing of the front 
  h0=pparam[6]; // offset height aspect ratio: usually = 1
  // sigma of the skin gaussian profile (in and front to all assymetry)
  skinsigmain=pparam[8]; // -- inner sigma
  skinsigmafr=pparam[9]; // -- front sigma
	
  // -- center of the circle depending on the alpha and rf param
  // rs : dist solar center - shell skeleton circle center
  rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
  // rc : radius shell skeleton circle
  rc=rf*tan(alpha);

//	std::cout << "rs : " << rs << std::endl;
//	std::cout << "rc : " << rc << std::endl;
 	
	
  //cvo=Cvec(0,0,0);
  cva=Cvec(sin(alpha),cos(alpha),0);	
  //skinthresup=skinsigmafr*3;
  //skinthresdw=skinsigmain*3;
}

void CModel53::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
  flagcase=0;
  vp.push_back(moddefparam("","Tube shell model, Gaussian skin profile.","",""));	
  vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));	
  vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));	
  vp.push_back(moddefparam("rf","10.","dist junction line-circle","Rsun"));	
  vp.push_back(moddefparam("ratio","0.2","ratio of tube radius to height",""));	
  vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));	
  vp.push_back(moddefparam("wrapcoeff","0.2","Front to back decreasing coeff: OBSOLETE: LEFT FOR BACKWARD COMPATIBILITY WITH MODEL 45.",""));	
  vp.push_back(moddefparam("h0","1.","Offset height for the aspect ratio: 1. is at the solar limb.","Rsun"));	
  vp.push_back(moddefparam("stiffness","1.","stiffness of the sigmoidal transition: UNUSED",""));	
  vp.push_back(moddefparam("skinsigmain","0.01","inner sigma",""));	
  vp.push_back(moddefparam("skinsigmafr","0.1","front sigma",""));	
  
  return;
}


// -- density 54
//    CGS model, as published in ApJs Volume 194, Issue 2, article id. 33, 6 pp. (2011)
float CModel54::Density(const Cvec &v) {
  // -- I shift the order of the axis to put the plane of the loop in the O,x,y plane
  //    I take the absolute value of x to use O,y,z plane symmetry
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0; // -- output density
  float diffinout,R;

  // -- return nothing if below requested height
  if (y < 0 || rr <= rb) return nel;

  // -- construct x,y,z vector
  Cvec V_OP(x,y,z);
  
 // -- on which side are we ? foot or shell ?
 float testside=y+x*tan_alpha-b;

 if (testside < 0) {
   // ---- foot side	 
	 // -- projection on foot axis
     Cvec V_OQ=orthoproj(O,Vleg,V_OP);
	 // -- height cross section center on foot
     
	 // -- compute cross section radius
	 R=V_OQ.norm()*tan_delta;

	 // -- dist point - foot axis
	 float QP=(V_OP-V_OQ).norm();
     
	 // -- put some density in the axis if requested
	 if (QP <= thick) return neaxis;
	 
   // -- take into account the skin thickness
   //    no density further than 3 sigmas
   if (QP < (R+skinsigmafrcut) && QP > (R-skinsigmaincut)) {
     // -- projection point on foot
     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit
     //    also check the side to allow asymmetry
     diffinout=R-QP;
     nel=nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

   }

   return(nel);

 } else {
    // ---- shell side
    // -- angular position
    float beta=atan2(y-b,x);

    // -- compute center shell
	float X0=(rho+b*k2*sin(beta))/mk2;
	Cvec V_OC(X0*cos(beta),b+X0*sin(beta),0);
	
    // -- compute dist Cprime-P
    float CP=(V_OP-V_OC).norm();
    
    // -- put some density in the axis if requested
	if (CP <= thick) return neaxis;

	// -- compute radius of shell
	R=sqrt(rc1+X0*X0);

	// -- take into account the skin thickness   
    //    no density further than 3 sigmas
   if ((CP < R+skinsigmafrcut) && (CP > R-skinsigmaincut)) {

     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit 
     //    also check the side to allow assymetry
     diffinout=R-CP;
     nel=nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

   }
   return(nel);
 }
 return 0;
} 

void CModel54::initParam(float* pparam) {
	rb = pparam[0]; //2.55; //dist to bottom of structure // 
	alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
    h=pparam[2];  // leg height [Rsun]
    kappa = pparam[3]; //ratio of tube radius to height //.2 // this is K
    nemin=pparam[4]; // electron density 
	thick=pparam[5]; // thickness of the skeleton: set to 0. if no display requested 
    
    neaxis = pparam[6]; // electron density in the skeleton (or axis). 
    
    // -- compute useful parameters
    k2=kappa*kappa;
    mk2=1-k2;
    tan_delta=kappa/sqrt(mk2);
    b=h/cos(alpha); 
    tan_alpha=tan(alpha);
	rho=h*tan_alpha;
	rc1=(b*b*k2-rho*rho)/mk2;
    
	// -- Thickness of the skin: pseudo-gaussian profile
	skinsigmain=pparam[8]; // -- inner sigma
	skinsigmafr=pparam[9]; // -- front sigma
    
    skinsigmaincut=skinsigmain*3;
    skinsigmafrcut=skinsigmafr*3;
    
    // -- Origin
    O=Cvec(0,0,0);
    
    // -- vector director of the legs
    Vleg=Cvec(sin(alpha),cos(alpha),0);
}

void CModel54::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{	
	flagcase=0;
	vp.push_back(moddefparam("","Tube shell model, Gaussian skin profile.","",""));	
	vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));	
	vp.push_back(moddefparam("alpha","0.52","angle between Oy and leg axis","rad"));	
	vp.push_back(moddefparam("h","6.","leg height","Rsun"));	
	vp.push_back(moddefparam("kappa","0.4","aspect ratio",""));	
	vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));	
	vp.push_back(moddefparam("thick","0.1","Thickness of the skeleton axis.","Rsun"));
	vp.push_back(moddefparam("neaxis","0.","Electron density in the skeleton or axis","cm^-3"));
	vp.push_back(moddefparam("stiffness","0.","NOT USED",""));	
	vp.push_back(moddefparam("skinsigmain","0.1","inner sigma",""));	
	vp.push_back(moddefparam("skinsigmafr","0.1","front sigma",""));	

	return;
}


// -- density 55
//      
float CModel55::Density(const Cvec &v) {
  // float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  // I've shifted the order of the axis to have O,x,y in the plane of the loop
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0; // output density

  // -- return nothing if below requested height
  if (rr <= rb) return nel;

  // -- Oxz is plane of symmetry
  // -- deal only with first quadrant
  if (y < 0) return nel;

  // -- on which side are we ? foot or shell ?
  float testside=y+x*tan(alpha)-rs;

  if (testside < 0) {
    // ---- foot side          
    // -- dist point - foot axis
    float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     
    // -- put some density in the axis
    float foo=pow(r,exponent)/thick;
    if (foo < 80) return nemin*exp(-foo); else return 0;

    //if (r <= thick) return nemin; else return 0;

  } else {
    // -- shell side
    // angular position
    float rrr,beta;
    cvcoord(x,y-rs,&rrr,&beta);
     
    // -- compute center shell
    float x0=(rc+rs*k2*sin(beta))/mk2;
    //float x0=(rc+(rs-h0)*k2*sin(beta))/mk2;
    Cvec vCp(x0*cos(beta),rs+x0*sin(beta),0);
    // -- compute dist Cprime-P
    Cvec vCpP=Cvec(x,y,z)-vCp;

    float foo=pow(vCpP.norm(),exponent)/thick;
    if (foo < 80) return nemin*exp(-foo); else return 0;
    
    //if (vCpP.norm() <= thick) return nemin; else return 0;
    
  }
  return 0;
} 


void CModel55::initParam(float* pparam) {
  rb = pparam[0]; //2.55; //dist to bottom of structure
  alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
  rf=pparam[2]; //10.; // dist junction line-circle [Rsun]
  ratio = pparam[3]; //ratio of tube radius to height //.2
  k2=ratio*ratio;
  mk2=1-k2;
  
  nemin=pparam[4]; // electron density 
  thick=pparam[5]; // thickness of the skeleton: set to 0. if no display requested 
  exponent=pparam[6]; // 2; exponent of the tube profile : exp(-abs(x ^ exponent))  
  
  // -- center of the circle depending on the alpha and rf param
  // rs : dist solar center - shell skeleton circle center
  rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
  rs2=rs*rs;

  // rc : radius shell skeleton circle
  rc=rf*tan(alpha);
    
  rc1=(rs2*k2-rc*rc)/mk2;
}

void CModel55::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{   
  flagcase=0;
  vp.push_back(moddefparam("","Tube shell model, Gaussian skin profile.","","")); 
  vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));    
  vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));  
  vp.push_back(moddefparam("rf","6.","dist junction line-circle","Rsun"));    
  vp.push_back(moddefparam("ratio","0.4","ratio of tube radius to height",""));   
  vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));  
  vp.push_back(moddefparam("thick","0.5","Thickness of the skeleton axis.","Rsun"));  
  vp.push_back(moddefparam("exponent","2.","exponent of the skeleton cross section","")); 
  return;
}




// -- density 56
// Density 56 bended streamer
//      
float CModel56::Density(const Cvec &v) {
  // -- convert in spherical coord
  float r,phi,theta;
  cvcoord(v[1],v[2],v[0],&r,&phi,&theta);

  //if (phi < PI/4 || phi >= 3*PI/4) return 0.;
  if (phi > (PI/2+hang) || phi < (PI/2-hang)) return 0.;

  // -- compute the theta0
  float theta0s=theta-thetam/(1+exp((h-r)*s));
  
  // -- distance is heigth time angle
  float dist=fabs(theta0s*r);
      
  // -- put density only close to the streamer
  if (dist < hthick) {
    
    float c[4]={1.34e4,1.15e6,-6.022e6,5.577e7};// valid from 2.5 to 20 Rsun
    float coef[4]={-1,-2,-3,-4};
    float nel=0;
    for(int i=0;i<=3;i++) {
      float powr=pow(r,coef[i]);
      nel+=c[i]*powr;
    }
    
    return nel*nefactor; 
  } else  return 0.;
  
  return 0;
} 

// Inititialization of the parameters for the Model 56
void CModel56::initParam(float* pparam) {
  
  h=pparam[0]; //- transition height
  s=pparam[1]; // - transition stiffness
  thetam=pparam[2]; //- deviation height
  theta0=pparam[3]; // - position at the origin
  hthick=pparam[4]; // - half thickness of the streamer
  hang=pparam[5]; // - half angular width
  nefactor=pparam[6]; // - Ne scaling factor

}

void CModel56::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{   
  flagcase=0;
  vp.push_back(moddefparam("","test of bended streamer","","")); 
  vp.push_back(moddefparam("h","6.","height of the sigmoidal transition","Rsun"));    
  vp.push_back(moddefparam("s","1.","stiffness of the transition",""));  
  vp.push_back(moddefparam("thetam","0.174533","deviation","rad"));    
  vp.push_back(moddefparam("theta0","0.","angular position at the origin",""));   
  vp.push_back(moddefparam("hthick","0.5","thickness of the streamer","Rsun"));  
  vp.push_back(moddefparam("hang","1.5708","half angular width of the slab","rad"));  
  vp.push_back(moddefparam("nefactor","1.","Ne scaling factor",""));  
  return;
}




// -- density 57
// Density 57: constant and uniform density
float CModel57::Density(const Cvec &v)
{
  return unifdens;
} 
float CModel57::Density(const Cvec &v,float &temperature)
{
	temperature=uniftemp;
  return unifdens;
} 

// Inititialization of the parameters for the Model 57
void CModel57::initParam(float* pparam) {
  unifdens=pparam[0]; // - constant density, specified by the user
	uniftemp=pparam[1];
}

void CModel57::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","constant and uniform density","","")); 
  vp.push_back(moddefparam("unifdens","1.","Constant density","electron/cm^3"));
  vp.push_back(moddefparam("uniftemp","1.","Constant temperature","K"));
  return;
}


// -- density 58
// Blob
float CModel58::Density(const Cvec &v)
{
  Cvec p=c-v;
  if (p.mag() <= hwidth) return blobdens; else return 0.;
}

// Inititialization of the parameters for the Model 58
void CModel58::initParam(float* pparam) {
  
  blobdens=pparam[0]; // - constant density, specified by the user
  hwidth=pparam[1];
  //c=Cvec(pparam[2],pparam[3],pparam[4]);
  
  lon=pparam[2];
  lat=pparam[3];
  alt=pparam[4];
  
  c=carrington2cart(PI-lon,lat,alt);
  
}

void CModel58::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{   
  flagcase=0;
  vp.push_back(moddefparam("","Blob","","")); 
  vp.push_back(moddefparam("blobdens","1.","Density of the blob","electron/cm^3"));
  vp.push_back(moddefparam("hwidth","0.5","Half width","Rsun"));
  /*
  vp.push_back(moddefparam("x","0.","Carrington x pos","rad"));
  vp.push_back(moddefparam("y","10.","Carrington y pos","rad"));
  vp.push_back(moddefparam("z","1.","Carrington z pos","Rsun"));
  */
   
  vp.push_back(moddefparam("lon","0.","Carrington Longitude","rad"));
  vp.push_back(moddefparam("lat","0.","Carrington Latitude","rad"));
  vp.push_back(moddefparam("alt","10.","height","Rsun"));
  
  return;
}


// -- density 59
// Blob and tail
float CModel59::Density(const Cvec &v) {
  
  // ---- check if in the nucleus
  Cvec p=c-v;
  if (p.mag() <= hwidth) return blobdens;
  
  // ---- compute car longitude of the requested point
  float lonlatrad[3];
  cart2carrington(v,lonlatrad);
  lonlatrad[0]=PI-lonlatrad[0];
  
  // ---- compute corresponding time using longitude fit
  float tp=ComputeTime(lonlatrad[0]);

  // ---- compute position of the dust center using time difference
  Cvec pdcntr=c+speed*(timenucpos-tp)*v/v.mag();

  p=pdcntr-v;
  
  if (p.mag() <= (hwidth)) return (blobdens); else return 0.;
}


float CModel59::ComputeTime(const float &lon) 
{

  float tp=pc[0],lonexp=1.;
  for (unsigned int i=1;i<4;i++) {
    tp+=pc[i]*lonexp*lon;
    lonexp*=lon;
  }
  return tp;
}
      
      
// Inititialization of the parameters for the Model 59
void CModel59::initParam(float* pparam) {
  
  blobdens=pparam[0]; // - constant density, specified by the user
  hwidth=pparam[1];
  //c=Cvec(pparam[2],pparam[3],pparam[4]);
  
  lon=pparam[2];
  lat=pparam[3];
  alt=pparam[4];
  speed=pparam[5];
  
  
  c=carrington2cart(PI-lon,lat,alt);
  
  // polynomial coeff of comet Encke time in day vs longitude fit
  pc[0]=17.220938;pc[1]=-6.2051115;pc[2]=-1.3204229;pc[3]=0.42679229;
  
  
  timenucpos=ComputeTime(lon);

  std::cout << "timenucpos : " << timenucpos << std::endl;
  std::cout << "lon : " << lon << std::endl;

}

void CModel59::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{   
  flagcase=0;
  vp.push_back(moddefparam("","Blob","","")); 
  vp.push_back(moddefparam("blobdens","1.","Density of the blob","electron/cm^3"));
  vp.push_back(moddefparam("hwidth","0.5","Half width","Rsun"));
  /*
  vp.push_back(moddefparam("x","0.","Carrington x pos","rad"));
  vp.push_back(moddefparam("y","10.","Carrington y pos","rad"));
  vp.push_back(moddefparam("z","1.","Carrington z pos","Rsun"));
  */
   
  vp.push_back(moddefparam("lon","0.","Carrington Longitude","rad"));
  vp.push_back(moddefparam("lat","0.","Carrington Latitude","rad"));
  vp.push_back(moddefparam("alt","10.","height","Rsun"));
  vp.push_back(moddefparam("speed","10.","Ejection speed","Rsun/day"));
  //vp.push_back(moddefparam("timenucpos","1","Time corresponding to the nucleus pos.","Day"));
  


  return;
}




// -- density 60
// Tube shell model based on density54 but the density does not go down to 0
// within the shell
float CModel60::Density(const Cvec &v) {
  // Tres tres bourrin: I've shifted the order of the axis 
  float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0; // output density

    // -- return nothing if below requested height
  if (rr <= rb) return nel;

 // -- Oxz is plane of symmetry
 // -- deal only with first quadrant
  if (y < 0) return nel;

 // -- on which side are we ? foot or shell ?
  float testside=y+x*tan(alpha)-rs;
  float diffinout,rcross;

  if (testside < 0) {
   // ---- foot side     
     // -- projection on foot axis
    Cvec Vc=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     // -- height cross section center on foot
    float h=Vc.norm();

     // -- compute cross section radius
    float rcross=sqrt(h*h/(1./k2-1));
     
     //float rcross=sqrt((h-h0)*(h-h0)/(1./k2-1));
     
     // -- dist point - foot axis
    float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     
     // -- put some density in the axis if requested
    if (r <= thick) return 1e4;
     
   // -- take into account the skin thickness
   // no density further than 3 sigmas
    if (r > rcross+skinsigmafr*3) return 0;
     //if (r < rcross+skinsigmafr*3 && r > rcross-skinsigmain*3) {
     // -- projection point on foot
      Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     
     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit
     //    also check the side to allow asymmetry
      diffinout=rcross-r;
      if (diffinout >= 0) {
        if (r > rcross-skinsigmain*3)
          nel=nemin*(exp(-pow(diffinout/skinsigmain,2))+neinside)/(1.+neinside); else
              nel=nemin*neinside/(1.+neinside);
      } else {
        nel=nemin*exp(-pow(diffinout/skinsigmafr,2));
      }
    //}

    return(nel);

  } else {
   // -- shell side
   // angular position
    float rrr,beta;
    cvcoord(x,y-rs,&rrr,&beta);
   
    // -- compute center shell
    float x0=(rc+rs*k2*sin(beta))/mk2;
    //float x0=(rc+(rs-h0)*k2*sin(beta))/mk2;
    Cvec vCp(x0*cos(beta),rs+x0*sin(beta),0);
    // -- compute dist Cprime-P
    Cvec vCpP=Cvec(x,y,z)-vCp;
    if (vCpP.norm() <= thick) return 1e4;

    // -- compute radius of shell
    rcross=sqrt(rc1+x0*x0);

    // -- take into account the skin thickness   
   // no density further than 3 sigmas
    if (vCpP.norm() > rcross+skinsigmafr*3) return 0.;
    //if (vCpP.norm() < rcross+skinsigmafr*3 && vCpP.norm() > rcross-skinsigmain*3) {

     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit 
     //    also check the side to allow assymetry
      diffinout=rcross-vCpP.norm();
      if (diffinout >= 0) {
        if (vCpP.norm() > rcross-skinsigmain*3) 
          nel=nemin*(exp(-pow(diffinout/skinsigmain,2))+neinside)/(1.+neinside); else
              nel=nemin*neinside/(1.+neinside);
      } else {
        nel=nemin*exp(-pow(diffinout/skinsigmafr,2));
      }
    //}
    return(nel);
  }
  return 0;
} 
//! Inititialization of the parameters for the Model 45
void CModel60::initParam(float* pparam) {
  rb = pparam[0]; //2.55; //dist to bottom of structure
  alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
  rf=pparam[2]; //10.; // dist junction line-circle [Rsun]
  ratio = pparam[3]; //ratio of tube radius to height //.2
  k2=ratio*ratio;
  mk2=1-k2;
  nemin=pparam[4]; // electron density 
  thick=pparam[5]; // thickness of the skeleton: set to 0. if no display requested 
  h0=pparam[6]; // offset height aspect ratio: usually = 0
    // -- center of the circle depending on the alpha and rf param
    // rs : dist solar center - shell skeleton circle center
  rs=rf*(cos(alpha)+tan(alpha)*sin(alpha));
  rs2=rs*rs;
    //rs2=(rs-h0)*(rs-h0);

    // rc : radius shell skeleton circle
  rc=rf*tan(alpha);
    //rc2=rc*rc;
    
  rc1=(rs2*k2-rc*rc)/mk2;
    // sigma of the skin gaussian profile (in and front to all assymetry)
  skinsigmain=pparam[8]; // -- inner sigma
  skinsigmafr=pparam[9]; // -- front sigma
  neinside=pparam[10]; // -- inside the shell Ne

}

void CModel60::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{   
  flagcase=0;
  vp.push_back(moddefparam("","Tube shell model, Gaussian skin profile.","","")); 
  vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));    
  vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));  
  vp.push_back(moddefparam("rf","6.","dist junction line-circle","Rsun"));    
  vp.push_back(moddefparam("ratio","0.4","ratio of tube radius to height",""));   
  vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));  
  vp.push_back(moddefparam("thick","0.1","Thickness of the skeleton axis.","Rsun"));  
  vp.push_back(moddefparam("h0","0.","NOT USED","Rsun")); 
  vp.push_back(moddefparam("stiffness","0.","NOT USED",""));  
  vp.push_back(moddefparam("skinsigmain","0.1","inner sigma",""));    
  vp.push_back(moddefparam("skinsigmafr","0.1","front sigma",""));    
  vp.push_back(moddefparam("neinside","0.7","Residual Ne within the shell",""));    

  return;
}

