/***************************************************************************
 *            raytrace.h
 *
 *  $Id: raytrace.h,v 1.16 2010-09-01 19:51:18 thernis Exp $
 ****************************************************************************/

/** \file raytrace.h
 * \brief Contains all the wrappers for IDL
 */
#ifndef RAYTRACE_H
#define RAYTRACE_H


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "constant.h"
#include "cuvemission.h"
#include "Cvec.h"
#include "Clos.h"
#include "CModelBase.h"
#include "Cdvoxel.h"
#include "rtmiscfunc.h"
#include "ModelPosition.h"

//! Raytracing program parameters, for parameter passing to raytracemain
struct rtparam
{
	int *pmodelid; //!< Model number
	int *pis; //!< X size of the image
	int *pjs; //!< Y size of the image
	float *pfovpix; //!< Pixel resolution [rad]
	float *pobspos; //!< Observer position [Rsun]
	float *pobsang; //!< Observer camera orientation [rad]
	float *pnepos; //!< Electron density center positon [Rsun]
	float *pneang; //!< Electron density orientation [rad]
	int *plosnbp; //!< LOS number of integration steps
	float *plosrange; //!< LOS integration range [Rsun]
	float *pmparam; //!< Model parameter array pointer
	int *pquiet; //!< Flag: Disable display of program progression
	float *pbtot; //!< Total brightness image pointer
	float *pbpol; //!< Polarized brightness image pointer
	float *pnetot; //!< Integrated electron density image pointer
	float *pcrpix; //!< Position of the Sun center on the image
	float *prhoim; //!< Image giving the impact parameter
	float *pmmlon; //!< Min and max longitude
	float *pmmlat; //!< Min and max latitude
	float *prrr; //!< Projected radius on the Earth plane of the sky
	int *pneonly; //!< Flag to set for Ne only calculation
	int *proi; //!< Region of interest map
	int *ppofinteg; //!< Force integration centered on a plane of integration
	float *ppoiang; //!< Plane of integration angle orientation
	float *Hlonlat; //!< Heliographic lon and lat of disk center in rad
	float *poccrad; //!< Radius of the occulter.
	float *padapthres; //!< Adaptative Simpson's integration difference threshold. No Simpson integration if = 0. 
	int *pmaxsubdiv; //!< Adaptative Simpson's integration maximun number of interval subdivision recursion. 
	float *plimbdark; //!< User limb darkening: default 0.58 
    float *protmat; //!< output final Ne rotation matrix
    float *obslonlat; //!< observer lon lat and height in carrington coordinate
    int *obslonlatflag; //!< set to 1 if user defined obslonlat: then obspos is ignored
    int *projtypecode; //!< type of projection: 1:ARC, 2:TAN, 3:SIN, 4:AZP
    float *pv2_1; //!< mu parameter for the AZP projection
    int *pfrontinteg; //!< Integration centered on the observer
    unsigned int uvinteg; //!< Use UV emission instead of Thomson scattering
	float disttofracmax; //!< Set to the fraction of max B to evaluate the integration distance to that threshold: the distance is returned in bpol
    float *pnerotcntr; //!< Electron density center positon within the nps coord system defined by nepos and neang [Rsun]
    float *pnerotang; //!< Electron density orientation within the nps coord system defined by nepos and neang [rad]
    float *pnetranslation; //!< Shift of the electron density, in the density coordinate system
	int *pnerotaxis; //!< Electron density axis ids corresponding to the nerotang rotation angles
	float *losDepthIn; //!< Start position (closer to observer) along the LOS when the density is not null
	float *losDepthOut; //!< End position (further away from observer) along the LOS when the density is not null
	int evalDepth; //!< Set to true in order to calculate losDepthIn and losDepthOut
};

//! Main raytracing routine in cartesian coordinates
int raytracemain(rtparam fp);
//! Wrapper of raytracemain called by IDL
extern "C" int raytracewl(int argc, void **argv);
//! Raytracing using WCS coordinate system
extern "C" int rtraytracewcs(int argc, void **argv);
//! Raytracing program that integrates only a cirular profile
extern "C" int rtwlcirc(int argc, void **argv);
//! raytracing program that integrates only along a straight line segment profile
extern "C" int rtwlseg(int argc, void **argv);
//! raytracing for only for a density cube
extern "C" int rtdenscube(int argc, void **argv);
//! return the carrington position on the solar limb for a given pixel: useful for EUVI
extern "C" int rtGetCarPosOnDisk(int argc, void **argv);


//! Compute ii,jj pix pos for the cube vertices
std::vector<std::vector<float> > voxelvert2iijj(std::vector<Cvec> &vert,const float &xstep,const float &ystep,const float &icntr,const float &jcntr);
//! Compute ii,jj pix pos for the cube vertices
std::vector<float> voxelvert2iijj(Cvec &vert,const float &xstep,const float &ystep,const float &icntr,const float &jcntr);

//! Convert LOS orientation in obs to pix coordinates
void losobs2xxyy(const Cvec &Vlosobs,float &xx,float &yy);
void losobs2xxyydebug(const Cvec &Vlosobs,float &xx,float &yy);



//! Convert pix coordinates to LOS orientation in obs
inline void xxyy2losobs(const float &xx,const float &yy,float &rrr,float &alpha) {
  rrr=sqrt(xx*xx+yy*yy);
  alpha=0;
  if (yy != 0) {
    alpha=atan(xx/yy);
    if (yy < 0)
      alpha=alpha+PI;
  }
}



//! Convert pix coordinates to LOS orientation in obs
inline Cvec xxyy2losobs(float const &xx,float const &yy) {
  float rrr,alpha;
  xxyy2losobs(xx,yy,rrr,alpha);
  return Cvec(sin(rrr)*sin(alpha),sin(rrr)*cos(alpha),cos(rrr));
}



//! Convert pix coordinates to LOS orientation in obs
inline void xxyy2losobs(const float &i,const float &j,float &rrr,float &alpha,const float *crpix,const float *fovpix) 
{
  float x=i-crpix[0];
  float y=j-crpix[1];  
  rrr=*fovpix*sqrt(x*x+y*y);
  alpha=atan2(y,x);
}

//! Convert pix coordinates to LOS orientation in obs
inline void xxyy2losobs(const float &i,const float &j,float &rrr,float &alpha,const float *crpix,const float *fovpix,const float *pc) 
{
  // ---- see Calabretta & Greisen, A&A, 2002
  float ic=i-crpix[0];
  float jc=j-crpix[1];
  
  float x=pc[0]*ic+pc[2]*jc; 
  float y=pc[1]*ic+pc[3]*jc;
   
  rrr=*fovpix*sqrt(x*x+y*y);
  alpha=atan2(y,x);
}




//! Convert LOS (in obs coord) into pix position
inline void losobs2xxyy(const Cvec &Vlosobs,const float *crpix,const float *fovpix,const float *pci,const int &projtypecode,const float &pv2_1,float &i, float &j)
{
  float phi=atan2( Vlosobs.v[0], Vlosobs.v[1]);
//  float theta0=acos( Vlosobs.v[2] / Vlosobs.norm());
  float theta0=atan2(sqrt(Vlosobs.v[0]*Vlosobs.v[0]+Vlosobs.v[1]*Vlosobs.v[1]),Vlosobs.v[2]);
  float theta=applyprojection(theta0,projtypecode,pv2_1);
  // see 10 Jul 2007
  float r=theta / (*fovpix);
  float x=r*cos(phi);
  float y=r*sin(phi);
  
  float ic=pci[0]*x+pci[2]*y; 
  float jc=pci[1]*x+pci[3]*y;

  i=ic+crpix[0];
  j=jc+crpix[1];
  
}

//! Convert pix coordinates to LOS orientation in obs
inline Cvec xxyy2losobs(const float &i,const float &j,float &rrr,float &alpha,const float *crpix,const float *fovpix,const float *pc,const int &projtypecode,const float &pv2_1)
{
  // ---- see Calabretta & Greisen, A&A, 2002
  float ic=i-crpix[0];
  float jc=j-crpix[1];
  
  float x=pc[0]*ic+pc[2]*jc; 
  float y=pc[1]*ic+pc[3]*jc;
   
  rrr=*fovpix*sqrt(x*x+y*y);
  alpha=atan2(y,x);
  
  rrr=applyinverseprojection(rrr,projtypecode,pv2_1);
  
  return Cvec(sin(rrr)*sin(alpha),sin(rrr)*cos(alpha),cos(rrr));
}

//! returns the determinant and the two roots (if det >= 0) of a 2nd deg polynomial. a0+a1x+a2x^2
inline float get2ndDegPolyRoots(const float &a0,const float &a1,const float &a2,float &root1,float &root2) {
  float det=a1*a1-4.*a0*a2;
  if (det < 0.) return det;
  
  float sqrtdet=sqrt(det);
  root1=-(a1+sqrtdet)/(2.*a2);
  root2=(-a1+sqrtdet)/(2.*a2);
  
  return det;
}

//! Compute impact parameter and LOS direction
inline void getimpactandlos(const float &xx,const float &yy,const Cbasis &obs,const Cbasis &abs,Cvec &vlosabs,float &rho) {
  float rrr;
  float alpha=0;
  xxyy2losobs(xx,yy,rrr,alpha);

  vlosabs=obs.ui * Cvec(sin(rrr)*sin(alpha),sin(rrr)*cos(alpha),cos(rrr)); // 3

    // -- compute rho : dist LOS - Sun cntr: impact parameter
  rho=psldist(obs.o,vlosabs,abs.o);
}




//! Compute impact parameter and LOS direction
inline void getimpactandlos(const float &xx,const float &yy,const Cbasis &obs,const Cbasis &abs,const int projtypecode,const float pv2_1,Cvec &vlosabs,float &rho) {
  float rrr;
  float alpha=0;
  
  // ---- compute LOS according to pixel position
  xxyy2losobs(xx,yy,rrr,alpha);
  // -- apply inverse WCS projection
  rrr=applyinverseprojection(rrr,projtypecode,pv2_1);

  vlosabs=obs.ui * Cvec(sin(rrr)*sin(alpha),sin(rrr)*cos(alpha),cos(rrr)); // 3

  // -- compute rho : dist LOS - Sun cntr: impact parameter
  rho=psldist(obs.o,vlosabs,abs.o);
}

//! Compute impact parameter and LOS direction
inline void getimpactandlos(const float &i,const float &j,const Cbasis &obs,const Cbasis &abs,const int projtypecode,const float pv2_1,Cvec &vlosabs,float &rho,const float *crpix,const float *fovpix,const float *pc) {
  float rrr;
  float alpha=0;
  
  // ---- get the los according to the pixel position
  xxyy2losobs(i,j,rrr,alpha,crpix,fovpix,pc);
  // -- apply inverse WCS projection
  rrr=applyinverseprojection(rrr,projtypecode,pv2_1);
  
  vlosabs=obs.ui * Cvec(sin(rrr)*sin(alpha),sin(rrr)*cos(alpha),cos(rrr));
  // -- compute rho : dist LOS - Sun cntr: impact parameter
  rho=psldist(obs.o,vlosabs,abs.o);
}

//! Compute impact parameter and LOS direction
inline void getimpactandlos(const float &i,const float &j,const Cbasis &obs,const Cbasis &abs,const int projtypecode,const float pv2_1,Cvec &vlosabs,float &rho,float &rrr,float &alpha,Cvec &vlosobs,const float *crpix,const float *fovpix,const float *pc) 
{
  
  // ---- get the los according to the pixel position
  xxyy2losobs(i,j,rrr,alpha,crpix,fovpix,pc);
	printvar(*fovpix);
	printvar(pc[0]);

	printvar(rrr);
	printvar(alpha);

  // -- apply inverse WCS projection
  rrr=applyinverseprojection(rrr,projtypecode,pv2_1);
  
  vlosobs=Cvec(sin(rrr)*sin(alpha),sin(rrr)*cos(alpha),cos(rrr));
  vlosabs=obs.ui * vlosobs;
  // -- compute rho : dist LOS - Sun cntr: impact parameter
  rho=psldist(obs.o,vlosabs,abs.o);

	printvar(crpix[0]);	printvar(crpix[1]);	printvar(i);	printvar(j);



}


//! Compute the LOS corresponding to the image pixel position
inline Cvec iijj2losobs(const float &i,const float &j,const float &xstart,const float &ystart,const float &xstep,const float &ystep) {
  float xx=xstart+j*xstep;
  float yy=ystart+i*ystep;

  return xxyy2losobs(xx,yy);
}



//! Compute ii,jj pix pos for the voxel vertices
inline void voxelvert2iijj(std::vector<std::vector<float> > &vij,std::vector<Cvec> &vert,const float &xstep,const float &ystep,const float &icntr,const float &jcntr) {

  float xx,yy;
  for(unsigned int i=0;i<vert.size();i++) {
    losobs2xxyy(vert[i],xx,yy);

    vij[0][i]=yy/ystep+icntr;
    vij[1][i]=xx/xstep+jcntr;
  }
}


//! Compute ii,jj pix pos for the voxel vertices
inline void voxelvert2iijj(float *vij,const Cvec &v,const float &xstep,const float &ystep,const float &icntr,const float &jcntr) {

  float xx,yy;
  losobs2xxyy(v,xx,yy);

  vij[0]=yy/ystep+icntr;
  vij[1]=xx/xstep+jcntr;

}


//! Model call wrapper for Dr Marque who likes Fortran
extern "C" float modelwrapper(int& modelid,
	float& x,float& y,float& z,
	float& pmparam,float& temperature);
	
//! Build a density cloud file for the frontend display
extern "C" int buildcloud(int argc, void **argv);

//! Get the electron density at a point [x,y,z] in space
extern "C" int getdensity(int argc, void **argv);

//! Get the projection in the image of a set of points in space
int getprojection(int &is,int &js,float &fovpix,
			float *obspos,float *obsang,
			float *nepos,float *neang,
			float *posin,float *posout,
			float *Hlonlat,int &nbpoints,
			bool flagclip);
			
//! Get the projection in the image of a set of points in space, IDL wrapper
extern "C" int getprojectionidl(int argc,void **argv);

//! Dump default parameters for IDL
extern "C" int getdefaultparameteridl(int argc,void **argv);

//! Check if the LOS chunk to integrate is not within or behind the sun
bool checkchunkwbsun(float rho,float r0, float r2,float vs0z,float vs2z);	

//! Adaptative Simpson's integration
void adaptabacksim(float rvar,Cvec vs,float f0ne,float f1ne,float f2ne,float losstep,Cvec vlosstep,int count,float &sum,float &sbtot,float &sbpol,CModelBase* pmod,float* pmodparam,float temperature,Cbasis nps,int &depth,int maxsubdiv,float rho);

//! Calculate the Thomson scattering geometrical factors
void getThomsonGeomFactor(const float &r,const float &rho,float &totterm,float &polterm);

//! Return the square of the number
inline float sqr(float &a) {
  return (a*a);
}

    
//! Return the size of a file
int FileSize(const char* sFileName);

//! Dump the compilation time
void dumpBuildInfo();

//! Dump compilation time for IDL
extern "C" int dumpbuildinfo(int argc,void **argv);


//! Test function for interfacing with Python
extern "C" int testParamPass(int param);


//! Compute the origin of the LOS, depending what the user specified
inline void calcqlos(const int &pofinteg,const int &frontinteg,const Cbasis &obs,const Cvec &vlosabs,const Cbasis &abs,Cbasis &poi,Cvec &qlos) {

// -- compute pos of ortho proj Sun cntr on the LOS
if (pofinteg == 0 && frontinteg == 0) qlos=orthoproj(obs.o,vlosabs,abs.o); else {
	if (pofinteg != 0) {

	// get normal vector to POI: defined as the z axis of POI coord sys
	Cvec poinormal=poi.ui.column(3);
	// compute intersection LOS - POI (see 25Aug2005 notes)
	qlos=obs.o - vlosabs * (pscal(obs.o,poinormal)/pscal(vlosabs,poinormal));
	/* Integration in the Earth POS (see 23Aug2005 otes)
	float tmpnorm=obs.o.norm();
	qlos=obs.o - vlosabs * ((tmpnorm*tmpnorm) / pscal(obs.o,vlosabs));
	*/
	}
	if (frontinteg !=0 ) {
		// ---- qlos is the observer
		qlos=obs.o;
	}
}
}




//! Integration along a line of sight
inline void losinteg(const int &pofinteg,
	const int &frontinteg,
	const Cbasis &obs,
	const Cvec &vlosabs,
	const Cbasis &abs,
	Cbasis &poi,
	Clos &los,
	const float &rho,
	CModelBase *pmod,
	float *pmodparam,
	const ModelPosition &modpos,
	const float &u,
	const float &constfactor,
	const int &flagneonly,
	float *posbtot,
	float *posbpol,
	float *posne,
    float *plosDepthIn,
    float *plosDepthOut,
    int evalDepth ) {

    
// -- Compute the origin of the LOS, depending what the user specified
Cvec qlos;
calcqlos(pofinteg,frontinteg,obs,vlosabs,abs,poi,qlos);

// -- save density in LOS profile
float *losProfile = new float[los.nbp];

// -- initialize
*posbtot=0;
*posbpol=0;
*posne=0;
*plosDepthIn=0;
*plosDepthOut=0;

// -- current position on the LOS
//float s=los.sstart;
Cvec vlosstep=vlosabs * los.ds;
Cvec vs=qlos + vlosabs * los.sstart;
Cvec vsInit;
vsInit = vs;
Cvec vrho=orthoproj(obs.o,vlosabs,abs.o);
//for(unsigned int k=0;k<los.nbp;k++,s+=los.ds,vs+=vlosstep) {
for(unsigned int k=0; k<los.nbp; k++, vs+=vlosstep) {

  // -- dist Sun cntr - LOS current pos
  float r=vs.norm();
  
  // -- integration only in front of the Sun

  if ( (r <= LIMBLIMIT) || (rho <= LIMBLIMIT && (vs-obs.o).mag() > sqrt(obs.o.mag()*obs.o.mag()-rho*rho) && pscal(vrho-obs.o,vs-obs.o) >= 0)) continue;


//  if (rho > LIMBLIMIT || (r > LIMBLIMIT && vs.v[2] > 0.)) 
    /******************************************/
    /*            MODEL CALL HERE             */
    /******************************************/
    // -- compute Ne    
    //float ner=pmod->Density((nps.u * (vs - nps.o))-nps.translation);
    float ner=pmod->Density(ChangetoDensityCoord(modpos,vs));
    losProfile[k] = ner;
    if (ner <= MIN_ELECTRON_DENSITY) continue;
     
    // total density
    *posne +=ner;
 
	// -- skip brightness calculation if requested
	if (flagneonly > 0) continue; 
		
    /******************************************/
    // -- compute geometrical functions
    float sinomega=RSUN/r;
    float sinsquareomega=sinomega*sinomega;
    float cossquareomega=1-sinsquareomega;
    float cosomega=sqrt(cossquareomega);
    
    float logterm=log((1.+sinomega)/cosomega);
    
    float a=cosomega*sinsquareomega;
    float b=-1./8.*(1.-3.*sinsquareomega-cossquareomega*((1.+3.*sinsquareomega)/sinomega)*logterm);
    
    float c=(4./3.)-cosomega-(cosomega*cossquareomega)/3.;
    float d=(1./8.)*(5.+sinsquareomega-cossquareomega*((5.-sinsquareomega)/sinomega)*logterm);
    
    // ---- sum in the pixel
    float rhooverr=rho/r;
    // the polarized brightness
    float polterm=(a+u*(b-a))*rhooverr*rhooverr;
    *posbpol +=ner*polterm;
    // the total brightness
    *posbtot +=ner*(2*(c+u*(d-c))-polterm);

  }

// multiply by the integral constant factor
*posbtot *=constfactor*los.ds;
*posbpol *=constfactor*los.ds;
*posne *=RSUN_CM*los.ds;

// -- check LOS profile and find begining and end of density
unsigned int losDensityStart, losDensityEnd;
unsigned int kk;

if (evalDepth > 0) {
    // -- search for start
    kk = 0;
    while ((losProfile[kk] <= 0.) && (kk < (los.nbp-1))) kk++;
    losDensityStart = kk;

    // -- search for end
    kk = los.nbp-1;
    while ((losProfile[kk] <= 0.) && (kk > 0)) kk--;
    losDensityEnd = kk;

    *plosDepthIn = ((vsInit + (vlosstep * losDensityStart)) - obs.o).mag();
    *plosDepthOut = ((vsInit + (vlosstep * losDensityEnd)) - obs.o).mag();
}

delete[] losProfile;
    
}


//! Check if a chunk is not within or behind the sun
inline bool checkchunkwbsun(float rho,float r0, float r2,float vs0z,float vs2z) {
    return (!((rho > LIMBLIMIT) || ((r0 > LIMBLIMIT && vs0z < 0) && (r2 > LIMBLIMIT && vs2z < 0))));
}



//! Calculate and return the Thomson scattering geometric functions
inline void getThomsonGeomFactor(const float &r,const float &rho,float &totterm,float &polterm) {
    float sinomega=RSUN/r;
    float sinsquareomega=sinomega*sinomega;
    float cossquareomega=1-sinsquareomega;
    float cosomega=sqrt(cossquareomega);

    float logterm=log((1.+sinomega)/cosomega);

    float a=cosomega*sinsquareomega;
    float b=-1./8.*(1.-3.*sinsquareomega-cossquareomega*((1.+3.*sinsquareomega)/sinomega)*logterm);

    float c=(4./3.)-cosomega-(cosomega*cossquareomega)/3.;
    float d=(1./8.)*(5.+sinsquareomega-cossquareomega*((5.-sinsquareomega)/sinomega)*logterm);

    float rhooverr=rho/r;
    // Polarized brightness term
    polterm=(a+U*(b-a))*rhooverr*rhooverr;
    // Total brightness term
    totterm=(2*(c+U*(d-c))-polterm);

}




//! Integration with adaptative step
inline void adaptabacksim(float rvar,Cvec vs,float f0ne,float f1ne,float f2ne,float losstep,Cvec vlosstep,int depthcount,float &sum,float &sbtot,float &sbpol,CModelBase* pmod,float* pmodparam,ModelPosition modpos,int &depth,int maxsubdiv,float rho) {
    float f01,f12,sum01=0.,sum12=0.,sbtot01=0.,sbtot12=0.,
                                            sbpol01=0.,sbpol12=0.;
    float tott[3],polt[3];

    // -- give up if too many recusions
    if (depthcount > maxsubdiv) {
        sum=losstep*(f0ne+4*f1ne+f2ne)/3;
        // -- get Thomson scattering geometric factors
        getThomsonGeomFactor(vs.norm(),rho,tott[0],polt[0]);
        getThomsonGeomFactor((vs+vlosstep).norm(),rho,tott[1],polt[1]);
        getThomsonGeomFactor((vs+(vlosstep*2)).norm(),rho,tott[2],polt[2]);

        sbtot=losstep*(f0ne*tott[0]+4*f1ne*tott[1]+f2ne*tott[2])/3;
        sbpol=losstep*(f0ne*polt[0]+4*f1ne*polt[1]+f2ne*polt[2])/3;

        //cout << "---- Max Depth : " << depthcount << endl;
        depth=depthcount;
        return;
    }

    // -- check if difference too big
    if (fabs(f2ne-f0ne) > rvar) {
        // -- divide the interval if too big
        // -- left
//        f01=pmod->Density(nps.u * ((vs+(vlosstep*0.5)) + nps.o));
        f01=pmod->Density(ChangetoDensityCoord(modpos,vs+(vlosstep*0.5)));

        adaptabacksim(rvar,vs,f0ne,f01,f1ne,losstep*0.5,vlosstep*0.5,depthcount+1,sum01,sbtot01,sbpol01,pmod,pmodparam,modpos,depth,maxsubdiv,rho);

        // -- divide the interval if too big
        // -- right
//        f12=pmod->Density(nps.u * ((vs+(vlosstep*1.5)) + nps.o));
        f12=pmod->Density(ChangetoDensityCoord(modpos,vs+(vlosstep*1.5)));

        adaptabacksim(rvar,vs+vlosstep,f1ne,f12,f2ne,losstep*0.5,vlosstep*0.5,depthcount+1,sum12,sbtot12,sbpol12,pmod,pmodparam,modpos,depth,maxsubdiv,rho);

        // -- sum both sides
        sum=sum01+sum12;
        sbtot=sbtot01+sbtot12;
        sbpol=sbpol01+sbpol12;

    } else {
        // ---- sum if OK
        // -- Ne
        sum=losstep*(f0ne+4*f1ne+f2ne)/3;

        // -- get Thomson scattering geometric factors
        getThomsonGeomFactor(vs.norm(),rho,tott[0],polt[0]);
        getThomsonGeomFactor((vs+vlosstep).norm(),rho,tott[1],polt[1]);
        getThomsonGeomFactor((vs+(vlosstep*2)).norm(),rho,tott[2],polt[2]);

        sbtot=losstep*(f0ne*tott[0]+4*f1ne*tott[1]+f2ne*tott[2])/3;
        sbpol=losstep*(f0ne*polt[0]+4*f1ne*polt[1]+f2ne*polt[2])/3;

        depth=depthcount;
    }
    return;
}



//! Integration along a line of sight with adaptative step
inline void losintegadaptstep(const int &pofinteg,
	const int &frontinteg,
	const Cbasis &obs,
	const Cvec &vlosabs,
	const Cbasis &abs,
	Cbasis &poi,
	Clos &los,
	const float &rho,
	CModelBase *pmod,
	float *pmodparam,
	const ModelPosition &modpos,
	const float &u,
	const float &constfactor,
	const int &flagneonly,
	const float &adapthres,
	const int &maxsubdiv,
	float *posbtot,
	float *posbpol,
	float *posne) {

Cvec qlos;
calcqlos(pofinteg,frontinteg,obs,vlosabs,abs,poi,qlos);

// -- initialize
*posbtot=0;
*posbpol=0;
*posne=0;

// -- current position on the LOS
float s=los.sstart;
Cvec vlosstep=vlosabs * los.ds;
Cvec vs=qlos + vlosabs * s;

// -- loop
bool flagfirstloop=true;
float f0ne=0,f1ne,f2ne;
for(unsigned int k=0;k<los.nbp;k++,s+=los.ds,vs+=vlosstep) {
	float r0=vs.norm();
	Cvec vs2=vs+vlosstep;
	float r2=vs2.norm();
	// -- check if chunk not within or behind the sun
	if (checkchunkwbsun(rho,r0,r2,vs.v[2],vs2.v[2])) continue;
	
	// -- if first loop then compute f0 (set a flag once done)
	if (flagfirstloop) {
		//f0ne=pmod->Density(nps.u * (vs + nps.o));
		f0ne=pmod->Density(ChangetoDensityCoord(modpos,vs));

		flagfirstloop=false;
	}
	
	// -- compute Ne : f1, f2
/*	f1ne=pmod->Density(nps.u * ((vs+(vlosstep*0.5)) + nps.o));
	f2ne=pmod->Density(nps.u * (vs2 + nps.o));*/
	f1ne=pmod->Density(ChangetoDensityCoord(modpos,vs+(vlosstep*0.5)));
	f2ne=pmod->Density(ChangetoDensityCoord(modpos,vs2));

	// -- call backtrack integrator
	float sumne02=0,stot=0,spol=0;
	int depth=0;
	adaptabacksim(adapthres,vs,f0ne,f1ne,f2ne,los.ds*0.5,vlosstep*0.5,1,sumne02,stot,spol,pmod,pmodparam,modpos,depth,maxsubdiv,rho);

	// -- add density
    	*posne +=sumne02;
	*posbtot +=stot;
	*posbpol +=spol;
	
	// -- switch last to first point
	f0ne=f2ne;
	
}
// multiply by the integral constant factor
*posbtot *=constfactor;
*posbpol *=constfactor;
*posne *=RSUN_CM;

}




//! Integration along a line of sight for UV light
inline void losintegUV(const int &pofinteg,
	const int &frontinteg,
	const Cbasis &obs,
	const Cvec &vlosabs,
	const Cbasis &abs,
	Cbasis &poi,
	Clos &los,
	const float &rho,
	CModelBase *pmod,
	float *pmodparam,
	const ModelPosition &modpos,
	const float &u,
	const float &constfactor,
	const int &flagneonly,
	float *posbtot,
	float *posbpol,
	float *posne,
	CUVEmission *uvemis,
	const int &uvinteg) {

// -- Compute the origin of the LOS, depending what the user specified
Cvec qlos;
calcqlos(pofinteg,frontinteg,obs,vlosabs,abs,poi,qlos);


// -- initialize
*posbtot=0;
*posbpol=0;
*posne=0;

// -- current position on the LOS
float s=los.sstart;
Cvec vlosstep=vlosabs * los.ds;
Cvec vs=qlos + vlosabs * s;
float temperature;

for(unsigned int k=0;k<los.nbp;k++,s+=los.ds,vs+=vlosstep) {

  // -- dist Sun cntr - LOS current pos
  float r=vs.norm();
  
  // -- integration only in front of the Sun
  if (rho > LIMBLIMIT || (r > LIMBLIMIT && vs.v[2] > 0.)) {
    /******************************************/
    /*            MODEL CALL HERE             */
    /******************************************/
    // -- compute Ne    
//    float ner=pmod->Density(nps.u * (vs + nps.o),temperature);
    float ner=pmod->Density(ChangetoDensityCoord(modpos,vs),temperature);
		float uv=uvemis->calcEmissivity(uvinteg,temperature);

    // total density
    *posne +=ner;
 
    *posbpol +=uv;
    *posbtot +=uv*ner*ner;

  }
}
// multiply by the integral constant factor
*posbtot *=los.ds;
*posbpol *=los.ds;
*posne *=RSUN_CM*los.ds;

}



//! Integration along a line of sight
inline void losintegdisttofracmax(const int &pofinteg,
	const int &frontinteg,
	const Cbasis &obs,
	const Cvec &vlosabs,
	const Cbasis &abs,
	Cbasis &poi,
	Clos &los,
	const float &rho,
	CModelBase *pmod,
	float *pmodparam,
	const ModelPosition &modpos,
	const float &u,
	const float &constfactor,
	const int &flagneonly,
	float *posbtot,
	float *posbpol,
	float *posne,
    float disttofracmax) {

// -- Compute the origin of the LOS, depending what the user specified
Cvec qlos;
calcqlos(pofinteg,frontinteg,obs,vlosabs,abs,poi,qlos);


// -- initialize: we return the half distance value in the bpol image
*posbpol=0;

// -- current position on the LOS
float bucket=0.;
Cvec vlosstep=vlosabs * los.ds;
Cvec vs=qlos + vlosabs * los.sstart;
Cvec vrho=orthoproj(obs.o,vlosabs,abs.o);
for(unsigned int k=0;k<los.nbp;k++,vs+=vlosstep) {

  // -- dist Sun cntr - LOS current pos
  float r=vs.norm();
  
    if ( (r <= LIMBLIMIT) || (rho <= LIMBLIMIT && (vs-obs.o).mag() > sqrt(obs.o.mag()*obs.o.mag()-rho*rho) && pscal(vrho-obs.o,vs-obs.o) >= 0)) continue;

  
  // -- integration only in front of the Sun
//  if (rho > LIMBLIMIT || (r > LIMBLIMIT && vs.v[2] > 0.)) {
    /******************************************/
    /*            MODEL CALL HERE             */
    /******************************************/
    // -- compute Ne    
	float temperature;
    float ner=pmod->Density(ChangetoDensityCoord(modpos,vs));
    if (ner <= 1e-1) continue;
    
    /******************************************/
    // -- compute geometrical functions
    float sinomega=RSUN/r;
    float sinsquareomega=sinomega*sinomega;
    float cossquareomega=1-sinsquareomega;
    float cosomega=sqrt(cossquareomega);
    
    float logterm=log((1.+sinomega)/cosomega);
    
    float a=cosomega*sinsquareomega;
    float b=-1./8.*(1.-3.*sinsquareomega-cossquareomega*((1.+3.*sinsquareomega)/sinomega)*logterm);
    
    float c=(4./3.)-cosomega-(cosomega*cossquareomega)/3.;
    float d=(1./8.)*(5.+sinsquareomega-cossquareomega*((5.-sinsquareomega)/sinomega)*logterm);
    
    // ---- sum in the pixel
    float rhooverr=rho/r;
    // the polarized brightness
    float polterm=(a+u*(b-a))*rhooverr*rhooverr;
    // the total brightness
   bucket +=ner*(2*(c+u*(d-c))-polterm)*constfactor*los.ds;
		if (bucket > (disttofracmax * (*posbtot))) break;
        *posbpol+=los.ds;
//  }
}
}



#endif

/*
* $Log: raytrace.h,v $
* Revision 1.16  2010-09-01 19:51:18  thernis
* Implement disttofracmax
*
* Revision 1.15  2009/04/13 21:03:57  thernis
* - Use ModelPosition class.
* - Implement extra positioning parameters for the models.
*
* Revision 1.14  2009/03/06 21:23:46  thernis
* Implement neshift
*
* Revision 1.13  2009/02/23 16:21:15  thernis
* Fix calculation of brightness in the UV renderer
*
* Revision 1.12  2009/02/09 20:51:20  thernis
* - Clean up the code
* - Change CModel::Density prototype
* - Update documentation
* - Implement multi-threading using boost thread
* - Add new models
*
* Revision 1.11  2008/09/30 18:33:29  thernis
* Replace include of losinteg.h and losintegadaptstep.h by inline functions.
*
* Revision 1.10  2008/09/23 14:08:13  thernis
* - add new testing models, - implement integration in front of the instrument
*
* Revision 1.9  2007/07/19 19:51:03  thernis
* Implement rtGetCarPosOnDisk
*
* Revision 1.8  2007/07/10 21:15:57  thernis
* Implement raytracing of a cloud of points.
*
* Revision 1.7  2007/05/14 17:19:41  thernis
* Add CVS id and log in all files
*
*/
