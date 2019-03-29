//
// $Id: models61to70.cpp,v 1.4 2010-09-01 15:33:00 thernis Exp $

#include "models61to70.h"
#include <iostream>
#include <cmath>
#include <vector>
#include "Cvec.h"
#include "constant.h"
#include "rtmiscfunc.h"

// -- density 61
// Density 61: Full spherical shell, cartesian coordinates
float CModel61::Density(const Cvec &v) {
	float r=(v-cntr).mag();
	if (r > radius) return 0.;
	if (r < radminusthick) return 0.;
	return dens;
} 

// Inititialization of the parameters
void CModel61::initParam(float* pparam) {
  dens=pparam[0]; // -- constant Ne for the shell
  radius=pparam[1]; // -- radius of the sphere
  thickness=pparam[2]; // -- Shell thickness
  cntr=Cvec(pparam[3],pparam[4],pparam[5]); // -- center
  radminusthick=radius-thickness;
}

void CModel61::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","Full sperical shell","","")); 
  vp.push_back(moddefparam("dens","1.","Density","electron/cm^3"));
  vp.push_back(moddefparam("radius","5.","Radius of the sphere","Rsun"));
  vp.push_back(moddefparam("thickness","1.","Shell thickness","Rsun"));
  vp.push_back(moddefparam("xcenter","0.","Center position X axis","Rsun"));
  vp.push_back(moddefparam("ycenter","0.","Center position Y axis","Rsun"));
  vp.push_back(moddefparam("zcenter","0.","Center position Z axis","Rsun"));
  return;
}





// Density 62: Spherical blob, spherical coordinates
float CModel62::Density(const Cvec &v) {
	float r=(c-v).mag();
	if (r > radius) return 0.;
	if (r < radminusthick) return 0.;
	return dens;
}

// Inititialization of the parameters
void CModel62::initParam(float* pparam) {
  dens=pparam[0];
  radius=pparam[1];
  thickness=pparam[2];
  lon=pparam[3];
  lat=pparam[4];
  alt=pparam[5];
  c=carrington2cart(PI-lon,lat,alt);
  radminusthick=radius-thickness;
}


void CModel62::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","Blob","","")); 
  vp.push_back(moddefparam("dens","1e5","Density of the blob","electron/cm^3"));
  vp.push_back(moddefparam("radius","5.","Radius","Rsun"));
  vp.push_back(moddefparam("thickness","1.","Half width","Rsun"));
  vp.push_back(moddefparam("lon","0.","Carrington Longitude","rad"));
  vp.push_back(moddefparam("lat","0.","Carrington Latitude","rad"));
  vp.push_back(moddefparam("alt","5.","height","Rsun"));
  return;
}



// Density 63: Full spherical shell, cartesian coordinates
float CModel63::Density(const Cvec &v)
{
	float r=v.mag();
	if (r > radius) return 0.;
	if (r < radminusthick) return 0.;
	return dens;
} 

// Inititialization of the parameters
void CModel63::initParam(float* pparam)
{
  dens=pparam[0]; // -- constant Ne for the shell
  radius=pparam[1]; // -- radius of the sphere
  thickness=pparam[2]; // -- Shell thickness
  radminusthick=radius-thickness;
}

void CModel63::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","Full sperical shell","","")); 
  vp.push_back(moddefparam("dens","1.","Density","electron/cm^3"));
  vp.push_back(moddefparam("radius","5.","Radius of the sphere","Rsun"));
  vp.push_back(moddefparam("thickness","1.","Shell thickness","Rsun"));
  return;
}

// Hollow cylinder
float CModel64::Density(const Cvec &v)
{
	//float rr=v.mag();

	//if (v[2] < 0.) return 0.;

	// -- compute theta
	float rrr,theta,xi,foo;
	//cvcoord(v.v[2],v.v[1],&rrr,&theta);

	theta=atan2(-v[1],v[2]);

	//if (theta > PI) theta-=TWOPI;

	//if (fabs(theta) > alpha) return 0.;

	// -- compute h
	float h=h0+b*theta*theta;

	
	// -- compute CM
	Cvec CM=v-Cvec(0.,-h*sin(theta),h*cos(theta));
	cvcoord(CM.v[2],CM.v[0],CM.v[1],&rrr,&xi,&foo);

//  	float cx2=(CM[1]*CM[1]+CM[2]*CM[2])/CM.magsqr();    //cos(xi);
//  	float cx=sqrt(cx2);

	float cx=cos(xi);
	float cx2=cx*cx;

 	float r=(h*k/(1.-k*k))*(k*cx+sqrt(1.-k*k*(1.-cx2)));

	if (flag == 1) {
		std::cout << "theta : " << theta << std::endl;
		std::cout << "h : " << h << std::endl;
		std::cout << "v : " << v << std::endl;
		std::cout << "OC : " << Cvec(0.,-h*sin(theta),h*cos(theta)) << std::endl;
		std::cout << "CM : " << CM << std::endl;
		std::cout << "xi : " << xi*RADEG << std::endl;
		std::cout << "r : " << r << std::endl;

		flag=0;
	}
//	if (v.mag() >= h && v.mag() <=h+thickness) return dens; else return 0.;


	//if (v.mag() >= h && v.mag() <= (h+thickness)) return dens; else return 0.;


	// -- compute xi
	//cvcoord(CM.v[2],CM.v[0],CM.v[1],&rrr,&xi,&foo);

	//xi=asin(CM[0]/CM.mag());

	// -- compute r



//	if (v.mag() >= h && v.mag() <= (h+thickness)) return dens; else return 0.;

	// -- compare
//	Cvec CM=v;
//	r=h;
	float CMmag=CM.mag();
	if ((CMmag >= r) && (CMmag <= r+thickness)) return dens; else return 0.;
//	if ((CMmag <= thickness)) return dens; else return 0.;
}

// Inititialization of the parameters
void CModel64::initParam(float* pparam) {
  hlead=pparam[0]; // -- height of the leading edge
  k=pparam[1]; // -- aspect ratio
  b=pparam[2]; // -- Flattening coeff of the tube axis
  alpha=pparam[3]; // -- half angle angular extent
  dens=pparam[4]; // -- constant Ne for the shell
  thickness=pparam[5]; // -- Shell thickness

//   mk=sqrt(1.-k*k);
//   h0=hlead*mk/(k+mk);

  h0=hlead*(1.-k);

	std::cout << "hlead : " << hlead << std::endl;
	std::cout << "k : " << k << std::endl;
	std::cout << "b : " << b << std::endl;
	std::cout << "alpha : " << alpha << std::endl;
	std::cout << "dens : " << dens << std::endl;
	std::cout << "thickness : " << thickness << std::endl;
//	std::cout << "mk : " << mk << std::endl;
	std::cout << "h0 : " << h0 << std::endl;
// 	std::cout << "Merde ! " << std::endl;

	flag=1;


}

void CModel64::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","Hollow cylinder","","")); 
  vp.push_back(moddefparam("hlead","5.","Leading edge height","Rsun"));
  vp.push_back(moddefparam("k","0.4","Aspect ratio",""));
  vp.push_back(moddefparam("b","1.","Flattening coeff of the axis","Rsun/rad^2"));
  vp.push_back(moddefparam("alpha","0.35","Half angular extent","rad"));
  vp.push_back(moddefparam("dens","1e5","Density","electron/cm^3"));
  vp.push_back(moddefparam("thickness","0.4","Shell thickness","Rsun"));
  return;
}

/*
void CModel64::wireFrame(float* pparam,float* vertparam,float* plist)
{

}
*/




// Model 65
// Cylindrical rod along the z axis.
float CModel65::Density(const Cvec &v)
{
float dist2z=sqrt(v[0]*v[0]+v[1]*v[1]);
if (dist2z < RodRadius && v[2] > 0) {
	float ne=ModelSaitoPolar->Density(v);
	for (unsigned int i=0;i<nbstrands;i++) {
		float dx=*(pstrandspos+i*2) - v[0];
		float dy=*(pstrandspos+i*2+1) - v[1];
		float dstrand=sqrt(dx*dx+dy*dy);
		if (dstrand < strandradius) return(20*ne);
	}

	return(ne);
	//return(ModelSaitoEquat->Density(v,pparam,temperature));
} 
return(0.);

//return(ModelSaitoPolar->Density(v,pparam,temperature));

}

// Inititialization of the parameters
void CModel65::initParam(float* pparam) {
// ---- init Saito polar and equatorial models
ModelSaitoPolar= new CModel20;
ModelSaitoPolar->initParam(pparam);
ModelSaitoEquat= new CModel21;
ModelSaitoEquat->initParam(pparam);

RodRadius=pparam[0]; // -- rod radius in Rsun
nbstrands=(unsigned int) (pparam[1]); // -- number of strands
strandradius=pparam[2];
pstrandspos=pparam+3; // points to the x,y array of strand position: must be 2 * nbstrands size

this->pparam=pparam;

}

// Destructor
CModel65::~CModel65()
{
    delete ModelSaitoPolar;
    delete ModelSaitoEquat;
}


void CModel65::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","Oz cylindrical rod","","")); 
  vp.push_back(moddefparam("RodRadius","0.5","Rod radius","Rsun"));

  return;
}




// Model 66
// Wavy neutral sheet with archimedian spiral: Jokipii model, ApJ 1981
float CModel66::Density(const Cvec &v)
{
float r=v.mag(); 
float phi0romegav=-phi0+r*OmegasunoverVwind;
float phi=atan2(-v[1],v[2]);
//float theta=asin(pow(v[0]/r,2)*(v[0] < 0 ? -1 : 1));
float theta=asin(v[0]/r);
float thetasheet=asin(sinalpha*sin(phi+phi0romegav));

if ((r * sin(fabs(theta-thetasheet))) <= SheetThickness) return ModelSaitoEquat->Density(v);

return 0.;

}
// Destructor
CModel66::~CModel66()
{
delete ModelSaitoEquat;
}
// Inititialization of the parameters
void CModel66::initParam(float* pparam) {
alpha=pparam[0]; //15.*DTOR;
sinalpha=sin(alpha);
phi0=pparam[1]; //0.;
Omegasun=2.9e-6 ; //[rad.s^-1] angular rotation velocity of the Sun
Vwind=pparam[2]; //4e5 ; //[m.s^-1] wind velocity
OmegasunoverVwind=696000e3*Omegasun/Vwind; // also convert the Rsun in m

SheetThickness=pparam[3];

ModelSaitoEquat= new CModel21;

this->pparam=pparam;

}

// Default parameters for IDL
void CModel66::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","Jokipii wavy sheet","","")); 
  vp.push_back(moddefparam("alpha","0.2618","Amplitude","rad"));
  vp.push_back(moddefparam("phi0","0.","Phase shift","rad"));
  vp.push_back(moddefparam("Vwind","4e5","Solar wind velocity","m.s^(-1)"));
  vp.push_back(moddefparam("SheetThickness","0.5","Sheet thickness","rad"));

  return;
}




// Model 67
// CCMC irregular spherical density cube
float CModel67::Density(const Cvec &v,float &temperature)
{

// ---- get r lon lat from v
float ri,loni,lati;

cart2rlonlat(v,ri,loni,lati);

// cout << "ri : " << ri << endl;
// cout << "loni : " << loni << endl;
// cout << "lati : " << lati << endl;

// -- check if within the r range
temperature=0.;
if (ri < r[0] || ri > r[sr-1]) return 0.;

// -- find r index
short flagout;
float idshiftr;
unsigned int idr=searchnearestindex02(ri,r,sr,flagout,idshiftr);
// unsigned int idr=searchnearestindex(ri,r,sr);

if (flagout != 0) return 0.;
// cout << "idr : " << idr << ", shift : "<< idshiftr<< endl;


// -- find lon index
float idshiftlon;
unsigned int idlon=searchnearestindex02(loni,lon,slon,flagout,idshiftlon);
// unsigned int idlon=searchnearestindex(loni,lon,slon);
if (flagout != 0) return 0.;
// cout << "idlon : " << idlon <<", shift : " << idshiftlon<< endl;

// -- find lat index
float idshiftlat;
unsigned int idlat=searchnearestindex02(lati,lat,slat,flagout,idshiftlat);
// unsigned int idlat=searchnearestindex(lati,lat,slat);
if (flagout != 0) return 0.;
// cout << "idlat : " << idlat << ",shift : "<<idshiftlat<<endl;

unsigned int loc=idr + idlon *sr + idlat * sr * slon;
// cout << "loc : " << loc << endl;


float neout=trilininterp(idshiftr,idshiftlon,idshiftlat,idr,idlon,idlat,nele,sr,slon);
temperature=trilininterp(idshiftr,idshiftlon,idshiftlat,idr,idlon,idlat,temp,sr,slon);

// float neout=*(nele + loc);
// temperature=*(temp + loc);

return neout;
}

float CModel67::Density(const Cvec &v)
{
float tempe=0.;
return this->Density(v,tempe);
}

// Inititialization of the parameters
void CModel67::initParam(float* pparam) {
sr=(unsigned int) pparam[0];
slon=(unsigned int) pparam[1];
slat=(unsigned int) pparam[2];
r=pparam+3;
lon=pparam+3+sr;
lat=pparam+3+sr+slon;
nele=pparam+3+sr+slon+slat;
temp=pparam+3+sr+slon+slat+(sr*slon*slat);

printvar(sr);
printvar(slon);
printvar(slat);
printvar(r[0]);
printvar(lon[0]);
printvar(lat[0]);
}




// Model 68
// Archimedian spiral
float CModel68::Density(const Cvec &v)
{
float r=sqrt(pow(v[1],2)+pow(v[2],2));
float theta=((r-a)/Vwind);
float nbturns=floor(theta / TWOPI);
theta-=TWOPI*nbturns;

float thetaspiral=atan2(v[1],v[2]);
if (thetaspiral < 0) thetaspiral+=TWOPI;

if (fabs(theta-thetaspiral) <= SheetThickness || fabs(theta-thetaspiral+TWOPI) <= SheetThickness || fabs(theta-thetaspiral-TWOPI) <= SheetThickness) return ModelSaitoEquat->Density(v);

return 0.;

}
// Destructor
CModel68::~CModel68()
{
delete ModelSaitoEquat;
}
// Inititialization of the parameters
void CModel68::initParam(float* pparam) {
a=pparam[0]; //15.*DTOR;

//Omegasun=2.9e-6 ; //[rad.s^-1] angular rotation velocity of the Sun
Vwind=pparam[1]/(2.9e-6 * 696000.); //0.;

//Vwind=pparam[2]; //4e5 ; //[m.s^-1] wind velocity
//OmegasunoverVwind=696000e3*Omegasun/Vwind; // also convert the Rsun in m

SheetThickness=pparam[2];

ModelSaitoEquat= new CModel21;

this->pparam=pparam;

}

// Default parameters for IDL
void CModel68::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{
  flagcase=0;
  vp.push_back(moddefparam("","Archimedean Spiral: r=a+b theta","","")); 
  vp.push_back(moddefparam("a","0","parameter a","Rsun"));
  vp.push_back(moddefparam("Vwind","500.","Solar wind speed","Km/s"));
  vp.push_back(moddefparam("SheetThickness","0.5","Sheet thickness","rad"));

  return;
}




// -- density 69
float CModel69::Density(const Cvec &v) {
  //  float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  // Tres tres bourrin: I've shifted the order of the axis 
  float x=v.v[1],y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0; // output density

    // -- return nothing if below requested height
  if (rr <= rb) return nel;

 // -- Oxz is plane of symmetry
 // -- deal only with first quadrant
 if (y < 0) return nel;

 // -- on which side are we ? foot or shell ?
 float testside=y+fabs(x)*tan(alpha)-rs;
 float diffinout,rcross;

 if (testside < 0) {
   // ---- foot side     
     // -- projection on foot axis
     Cvec Vc=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(fabs(x),y,z));
     // -- height cross section center on foot
     float h=Vc.norm();

     // -- compute cross section radius
     float rcross=sqrt(h*h/(1./k2-1));
     
     //float rcross=sqrt((h-h0)*(h-h0)/(1./k2-1));
     
     // -- dist point - foot axis
     float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(fabs(x),y,z));
     
     // -- put some density in the axis if requested
     if (r <= thick) return 1e4;
     
   // -- take into account the skin thickness
   // no density further than 3 sigmas
   if (r < (rcross+skinsigmafr*3) && r > (rcross-skinsigmain*3)) {
     // -- projection point on foot
     //Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     
     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit
     //    also check the side to allow asymmetry
     diffinout=rcross-r;
     
     // -- modulation of the Ne
     float polang=atan2(Vc[0],Vc[1]);
     float ModulationFactor=cos(polang*180./T_modulation);
     
     nel=(1.-A_modulation*ModulationFactor*ModulationFactor)*nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

   }

   return(nel);

 } else {
   // -- shell side
   // angular position
   float rrr,beta;
   //cvcoord(x,y-rs,&rrr,&beta);
   beta=atan2(y-rs,x);
   
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
   if (vCpP.norm() < rcross+skinsigmafr*3 && vCpP.norm() > rcross-skinsigmain*3) {

     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit 
     //    also check the side to allow assymetry
     diffinout=rcross-vCpP.norm();
     
     // -- modulation of the Ne
     float polang=atan2(vCp[0],vCp[1]);
     float ModulationFactor=cos(polang*180./T_modulation);
     
     nel=(1.-A_modulation*ModulationFactor*ModulationFactor)*nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

   }
   return(nel);   
 }
 return 0;
} 

void CModel69::initParam(float* pparam) {
    rb = pparam[0]; //2.55; //dist to bottom of structure // 
    alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
    rf=pparam[2]; //10.; // dist junction line-circle [Rsun] // this is h_f
    ratio = pparam[3]; //ratio of tube radius to height //.2 // this is K
    k2=ratio*ratio;
    mk2=1-k2;
    nemin=pparam[4]; // electron density 
    thick=pparam[5]; // thickness of the skeleton: set to 0. if no display requested 
    T_modulation=pparam[6]; // period of the Ne modulation.
    A_modulation=pparam[7]; // amplitude of the Ne modulation.

    // -- center of the circle depending on the alpha and rf param
    // rs : dist solar center - shell skeleton circle center
    rs=rf*(cos(alpha)+tan(alpha)*sin(alpha)); // this is h
    rs2=rs*rs;

    // rc : radius shell skeleton circle
    rc=rf*tan(alpha); // this is rho
    
    rc1=(rs2*k2-rc*rc)/mk2; // 
    // sigma of the skin gaussian profile (in and front to all assymetry)
    skinsigmain=pparam[8]; // -- inner sigma
    skinsigmafr=pparam[9]; // -- front sigma

}

void CModel69::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{   
    flagcase=0;
    vp.push_back(moddefparam("","GCS model with sin variation of Ne.","","")); 
    vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));    
    vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));  
    vp.push_back(moddefparam("rf","6.","dist junction line-circle","Rsun"));    
    vp.push_back(moddefparam("ratio","0.4","ratio of tube radius to height",""));   
    vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));  
    vp.push_back(moddefparam("thick","0.1","Thickness of the skeleton axis.","Rsun"));  
    vp.push_back(moddefparam("T_modulation","0.","Period of the Ne modulation","cycle")); 
    vp.push_back(moddefparam("A_modulation","1.","Amplitude of the Ne modulation",""));  
    vp.push_back(moddefparam("skinsigmain","0.1","inner sigma",""));    
    vp.push_back(moddefparam("skinsigmafr","0.1","front sigma",""));    

    return;
}








// -- density 70
float CModel70::Density(const Cvec &v) {
  //  float x=v.v[0],y=fabs(v.v[1]),z=v.v[2],rr=v.norm();
  // Tres tres bourrin: I've shifted the order of the axis 
  //float x=fabs(v.v[1]),y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float x=v.v[1],y=v.v[2],z=v.v[0],rr=sqrt(x*x+y*y+z*z);
  float nel = 0; // output density

    // -- return nothing if below requested height
  if (rr <= rb) return nel;

 // -- Oxz is plane of symmetry
 // -- deal only with first quadrant
 if (y < 0) return nel;

 // -- on which side are we ? foot or shell ?
 float testside=y+fabs(x)*tan(alpha)-rs;
 float diffinout,rcross;

 if (testside < 0) {
   // ---- foot side     
     // -- projection on foot axis
     Cvec Vc=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(fabs(x),y,z));
     // -- height cross section center on foot
     float h=Vc.norm();

     // -- compute cross section radius
     float rcross=sqrt(h*h/(1./k2-1));
     
     //float rcross=sqrt((h-h0)*(h-h0)/(1./k2-1));
     
     // -- dist point - foot axis
     float r=psldist(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(fabs(x),y,z));
     
     // -- put some density in the axis if requested
     if (r <= thick) return 1e4;
     
   // -- take into account the skin thickness
   // no density further than 3 sigmas
   if (r < (rcross+skinsigmafr*3) && r > (rcross-skinsigmain*3)) {
     // -- projection point on foot
     //Cvec vpf=orthoproj(Cvec(0,0,0),Cvec(sin(alpha),cos(alpha),0),Cvec(x,y,z));
     
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
   float rrr,beta;
   //cvcoord(x,y-rs,&rrr,&beta);
   
   beta=atan2(y-rs,x);
   
   
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
   if (vCpP.norm() < rcross+skinsigmafr*3 && vCpP.norm() > rcross-skinsigmain*3) {

     // -- gaussian profile: note it's not exactly the gaussian 
     //    probability function but it speeds up the code a tiny bit 
     //    also check the side to allow assymetry
     diffinout=rcross-vCpP.norm();
     
     // -- modulation of the Ne
     // float polang=atan2(vCp[1],vCp[0]);
     if (x < 0. && beta < 0.) beta+=TWOPI;
     float ModulationFactor=cos(beta*180./T_modulation);
     
     nel=(1.-A_modulation*ModulationFactor*ModulationFactor)*nemin*exp(-pow(diffinout/((diffinout >= 0) ? skinsigmain : skinsigmafr),2));

   }
   return(nel);   
 }
 return 0;
} 

void CModel70::initParam(float* pparam) {
    rb = pparam[0]; //2.55; //dist to bottom of structure // 
    alpha=pparam[1];//*DTOR; //30.*DTOR; // angle between axis and foot [rad]
    rf=pparam[2]; //10.; // dist junction line-circle [Rsun] // this is h_f
    ratio = pparam[3]; //ratio of tube radius to height //.2 // this is K
    k2=ratio*ratio;
    mk2=1-k2;
    nemin=pparam[4]; // electron density 
    thick=pparam[5]; // thickness of the skeleton: set to 0. if no display requested 
    T_modulation=pparam[6]; // period of the Ne modulation.
    A_modulation=pparam[7]; // amplitude of the Ne modulation.

    // -- center of the circle depending on the alpha and rf param
    // rs : dist solar center - shell skeleton circle center
    rs=rf*(cos(alpha)+tan(alpha)*sin(alpha)); // this is h
    rs2=rs*rs;

    // rc : radius shell skeleton circle
    rc=rf*tan(alpha); // this is rho
    
    rc1=(rs2*k2-rc*rc)/mk2; // 
    // sigma of the skin gaussian profile (in and front to all assymetry)
    skinsigmain=pparam[8]; // -- inner sigma
    skinsigmafr=pparam[9]; // -- front sigma

}

void CModel70::dumpDefaultParamForIDL(std::vector<moddefparam>& vp,int& flagcase)
{   
    flagcase=0;
    vp.push_back(moddefparam("","GCS model with sin variation of Ne.","","")); 
    vp.push_back(moddefparam("rb","2.55","dist to bottom of structure","Rsun"));    
    vp.push_back(moddefparam("alpha","0.52","angle between axis and foot","rad"));  
    vp.push_back(moddefparam("rf","6.","dist junction line-circle","Rsun"));    
    vp.push_back(moddefparam("ratio","0.4","ratio of tube radius to height",""));   
    vp.push_back(moddefparam("nemin","1e6","Ne","cm^-3"));  
    vp.push_back(moddefparam("thick","0.1","Thickness of the skeleton axis.","Rsun"));  
    vp.push_back(moddefparam("T_modulation","0.","Period of the Ne modulation","cycle")); 
    vp.push_back(moddefparam("A_modulation","1.","Amplitude of the Ne modulation",""));  
    vp.push_back(moddefparam("skinsigmain","0.1","inner sigma",""));    
    vp.push_back(moddefparam("skinsigmafr","0.1","front sigma",""));    

    return;
}




