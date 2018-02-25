//
// C++ Implementation: rtmiscfunc
//
// $Id: rtmiscfunc.cpp,v 1.1 2009/02/09 20:46:31 thernis Exp $
//



#include "rtmiscfunc.h"


float trilininterp(float x,float y,float z,
                   float xint[2],float yint[2],float zint[2],
                   float cube[2][2][2]) {

    float t=(x-xint[0])/(xint[1]-xint[0]);
    float u=(y-yint[0])/(yint[1]-yint[0]);
    float v=(z-zint[0])/(zint[1]-zint[0]);

    float r=(1-t)*(1-u)*(1-v)*cube[0][0][0]+
            (t)  *(1-u)*(1-v)*cube[1][0][0]+
            (1-t)*(u)  *(1-v)*cube[0][1][0]+
            (t)  *(u)  *(1-v)*cube[1][1][0]+
            (1-t)*(1-u)*(v)  *cube[0][0][1]+
            (t)  *(1-u)*(v)  *cube[1][0][1]+
            (1-t)*(u)  *(v)  *cube[0][1][1]+
            (t)  *(u)  *(v)  *cube[1][1][1];

    return(r);

}

float trilininterp(const float &t,const float &u,const float &v,
										const unsigned int &xi0,
										const unsigned int &yi0,
										const unsigned int &zi0,
										float *cube,
										const unsigned int &sx,const unsigned int &sy) {

unsigned int xi1=xi0+1;
unsigned int yi1=yi0+1;
unsigned int zi1=zi0+1;
unsigned int sxsy=sx*sy;


    float r=(1-t)*(1-u)*(1-v)*cube[xi0+yi0*sx+zi0*sxsy]+
            (t)  *(1-u)*(1-v)*cube[xi1+yi0*sx+zi0*sxsy]+
            (1-t)*(u)  *(1-v)*cube[xi0+yi1*sx+zi0*sxsy]+
            (t)  *(u)  *(1-v)*cube[xi1+yi1*sx+zi0*sxsy]+
            (1-t)*(1-u)*(v)  *cube[xi0+yi0*sx+zi1*sxsy]+
            (t)  *(1-u)*(v)  *cube[xi1+yi0*sx+zi1*sxsy]+
            (1-t)*(u)  *(v)  *cube[xi0+yi1*sx+zi1*sxsy]+
            (t)  *(u)  *(v)  *cube[xi1+yi1*sx+zi1*sxsy];

    return r;

}

//! Seeks the nearest neighbor points in
//! a neutral sheet map
int wherenn(float* pnsheetmap,int slon,int slat,int srlon,int srlat,float lon,float lat,float* plon,float* plat,float* pdist,float* pval) {

    int x,y,xp,yp;
    float dist,mindist=1e10;
    int out=0;
    float xx,yy,val;
    //  float xxx=1000,yyy=1000; // debug var

    lat=lat+90;

    // -- compute axis resolution
    float rlon=((float) slon) / 360;
    float rlat=((float) slat) / 181;

    //
    int hlonpix=(int) (rlon * 180);

    int intlon=(int) (lon * rlon);
    int intlat=(int) (lat * rlat);


    // -- scan all the position within the range and compute the minimum distance if there are any points
    for(int i=intlon-srlon;i<=intlon+srlon;i++) {
        for(int j=intlat-srlat;j<=intlat+srlat;j++) {

            // -- take into account boundary conditions
            //    if scanned longitude negative then
            //    wrap to the end of the map
            if (i < 0) {
                x=i+slon;
            } else
                if (i >= slon) {
                    x=i-slon;
                } else
                    x=i;

            //   if latitude negative then
            //      if in the first longitude half then
            //         go to the oposite quadran
            if (j < 0) {
                y=j+slat;
                if (x < hlonpix) {
                    x=x+hlonpix;
                } else
                    x=x-hlonpix;
            } else
                if (j>= slat) {
                    y=j-slat;
                    if (x < hlonpix) {
                        x=x+hlonpix;
                    } else
                        x=x-hlonpix;
                } else
                    y=j;

            /*
            if ((i < 0) || (i > slon)) {
            cout << "intlon : " << intlon
            << ", intlat : " << intlat << endl;

            //  cout << "i : " << i
            << ", j : " << j
            << ", x : " << x
            << ", y : " << y << endl;
            }

            */


            // -- quick test if not out of boundaries
            if ((x < 0) || (x >= slon) || (y < 0) || (y >= slat)) {
                std::cout << "X or Y out of range" << std::endl;

                return 0;

            }


            // -- pick up the value of the map at that point
            val=*(pnsheetmap+x+slon*y);

            if (val > 0.) {
                xx=((float) i)/rlon-lon;
                yy=((float) j)/rlat-lat;

                //cout << "x : " << x
                //     << ", y : " << y << endl;

                dist=xx*xx+yy*yy;
                if (dist < (mindist)) {
                    mindist=dist;
                    *plon=i;
                    *plat=j;
                    *pval=val;
                    xp=x;
                    yp=y;
                    //xxx=xx;
                    //yyy=yy;
                    out=1;
                }
            }
        }
    }

    if (out == 1)
        *pdist=sqrt(mindist);

    return out;
}


//! Seeks the nearest neighbor points in
//!      a neutral sheet map with a smoothing to try to avoid
//!      the pixelization effect
int wherennsmoothed(float* pnsheetmap,int slon,int slat,int srlon,int srlat,float lon,float lat,float* plon,float* plat,float* pdist,float* pval) {

    int x,y,xp=0,yp=0;
    float dist,mindist=1e10;
    int out=0;
    float xx,yy,val;
    float xxx=1000,yyy=1000; // debug var

    lat=lat+90;

    // -- compute axis resolution
    float rlon=((float) slon) / 360;
    float rlat=((float) slat) / 181;

    //
    int hlonpix=(int) (rlon * 180);

    int intlon=(int) (lon * rlon);
    int intlat=(int) (lat * rlat);


    // -- scan all the position within the range and compute the minimum distance if there are any points
    for(int i=intlon-srlon;i<=intlon+srlon;i++) {
        for(int j=intlat-srlat;j<=intlat+srlat;j++) {

            // -- take into account boundary conditions
            //    if scanned longitude negative then
            //    wrap to the end of the map
            if (i < 0) {
                x=i+slon;
            } else
                if (i >= slon) {
                    x=i-slon;
                } else
                    x=i;

            //   if latitude negative then
            //      if in the first longitude half then
            //         go to the oposite quadran
            if (j < 0) {
                y=j+slat;
                if (x < hlonpix) {
                    x=x+hlonpix;
                } else
                    x=x-hlonpix;
            } else
                if (j>= slat) {
                    y=j-slat;
                    if (x < hlonpix) {
                        x=x+hlonpix;
                    } else
                        x=x-hlonpix;
                } else
                    y=j;

            /*
            if ((i < 0) || (i > slon)) {
            cout << "intlon : " << intlon
            << ", intlat : " << intlat << endl;

            //  cout << "i : " << i
            << ", j : " << j
            << ", x : " << x
            << ", y : " << y << endl;
            }

            */


            // -- quick test if not out of boundaries
            if ((x < 0) || (x >= slon) || (y < 0) || (y >= slat)) {
                std::cout << "X or Y out of range" << std::endl;

                return 0;

            }


            // -- pick up the value of the map at that point
            val=*(pnsheetmap+x+slon*y);

            if (val > 0.) {
                xx=((float) i)/rlon-lon;
                yy=((float) j)/rlat-lat;

                //cout << "x : " << x
                //     << ", y : " << y << endl;

                dist=xx*xx+yy*yy;
                if (dist < (mindist)) {
                    mindist=dist;
                    *plon=i;
                    *plat=j;
                    *pval=val;
                    xp=x;
                    yp=y;
                    xxx=xx;
                    yyy=yy;
                    out=1;
                }
            }
        }
    }

    //if (out == 1) *pdist=sqrt(mindist);

    //  return out;
    if (out == 0)
        return out;

    //if (out > 0) {
    /*
    cout << " mindist : " << mindist
         << " plon : "  << *plon
         << " plat : " << *plat
         << " xp : "  << xp
         << " yp : " << yp
         << " xxx : "  << xxx
         << " yyy : " << yyy
         << " lon : "  << lon
         << " lat : " << lat << endl;
    //}
    */

    // ---- linear interpolation to compute the smallest dist

    // -- seek for the nearest neighbor of the neutral sheet detected point
    // - scan the 8 neighbors
    //
    //  012
    //  7X3
    //  654
    float d0=0,d1=0,d2=0,d3=0,d4=0,d5=0,d6=0,d7=0;
    int xscan=xp-1,yscan=yp-1;
    int cntnz=1;
    if (xscan >=0 && yscan >=0) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon-1))/rlon-lon;
            yy=((float) (*plat-1))/rlat-lat;
            d0=xx*xx+yy*yy;
            cntnz++;
        }
    }

    xscan=xp;
    if (yscan>=0) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon))/rlon-lon;
            yy=((float) (*plat-1))/rlat-lat;
            d1=xx*xx+yy*yy;
            cntnz++;
        }
    }

    xscan=xp+1;
    if (xscan < slon && yscan>=0) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon+1))/rlon-lon;
            yy=((float) (*plat-1))/rlat-lat;
            d2=xx*xx+yy*yy;
            cntnz++;
        }
    }

    xscan=xp-1;
    yscan=yp;
    if (xscan >=0) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon-1))/rlon-lon;
            yy=((float) (*plat))/rlat-lat;
            d7=xx*xx+yy*yy;
            cntnz++;
        }
    }

    xscan=xp+1;
    if (xscan < slon) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon+1))/rlon-lon;
            yy=((float) (*plat))/rlat-lat;
            d3=xx*xx+yy*yy;
            cntnz++;
        }
    }

    xscan=xp-1;
    yscan=yp+1;
    if (xscan >=0 && yscan < slat) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon-1))/rlon-lon;
            yy=((float) (*plat+1))/rlat-lat;
            d6=xx*xx+yy*yy;
            cntnz++;
        }
    }

    xscan=xp;
    if (yscan < slat) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon))/rlon-lon;
            yy=((float) (*plat+1))/rlat-lat;
            d5=xx*xx+yy*yy;
            cntnz++;
        }
    }


    xscan=xp+1;
    if (xscan < slon && yscan < slat) {
        if (*(pnsheetmap+xscan+slon*yscan) > 0.) {
            xx=((float) (*plon+1))/rlon-lon;
            yy=((float) (*plat+1))/rlat-lat;
            d4=xx*xx+yy*yy;
            cntnz++;
        }
    }

    /*
    cout << " d0 : " << d0
         << " d1 : " << d1
         << " d2 : " << d2
         << " d3 : " << d3
         << " d4 : " << d4
         << " d5 : " << d5
         << " d6 : " << d6
         << " d7 : " << d7 << endl;
    */


    // - find the minimum dist
    mindist=(mindist+d0+d1+d2+d3+d4+d5+d6+d7)/((float) cntnz);


    if (out == 1)
        *pdist=sqrt(mindist);

    return out;
}


//! Return the pixel value at the requested position on a PFSS map
int getPosOnSSMap(float *pnsheetmap,int sang,int slat,float phi,float theta,float& dist,float& val) {

    // ---- compute position of the point on the map
    int x=(int) ((0.5+phi)*(((float) sang)/360));
    int y=(int) ((theta+90.5)*(((float) slat)/181));

    if (x < 0)
        x=0;
    else if (x > (sang-1))
        x=sang-1;
    if (y < 0)
        y=0;
    else if (y > (slat-1))
        y=slat-1;

    val=*(pnsheetmap+y*sang+x);
    dist=fabs(val);

    return 1;

}



float densyming2(const float x,const float y,const float z,const float phi0,
                 const float *prco,const float *pphico,const float *pthetaco,
                 float *pdens) {


    // rco : [rintercept,rslope,sizer]
    // phico : [phiintercept,phislope,sizephi]
    // thetaco : [thetaintercept,thetaslope,sizetheta]
    // dens : cube of density with [sizer,sizephi,sizetheta]

    float neighborcube[2][2][2];

    // -- convert x,y,z to r,theta,phi

    float r,phi,theta;

    cvcoord(x,y,z,&r,&phi,&theta);
    /*
    float r=sqrt(x*x+y*y+z*z);
    float theta=asin(z/r);
    float ratio=abs(x)/(r*cos(theta));
    float phi=0;

    if (abs((int) ratio) < 1) {
      phi=acos(ratio);
    }
    if (x < 0) phi=PI-phi;
    if (y < 0) phi=TWOPI-phi;
    //phi+=PI+phi0;
    //phi=(float) modf(2*PI,&phi);

    */

    /*
    cout << "x : " << x << endl;
    cout << "y : " << y << endl;
    cout << "z : " << z << endl;


    cout << "r : " << r << endl;
    cout << "theta : " << theta << endl;
    cout << "phi : " << phi << endl;

    cout << "rco : " 
         << *prco << " , "
         << *(prco+1) << " , " 
         << *(prco+2) << endl;
    */

    // -- find the nearest neighbors
    int mr=(int)(- *prco/ *(prco+1) + r / *(prco+1) + 1);
    if (mr >= *(prco+2) || mr <= 0)
        return(0.);

    //  cout << "mr : " << mr << endl;

    int mtheta12[2];
    int mtheta=(int)(- *pthetaco/ *(pthetaco+1) + theta / *(pthetaco+1) + 1);
    if (mtheta >= *(pthetaco+2)) {
        mtheta12[0]=(int)*(pthetaco+2)-1;
        mtheta12[1]=mtheta12[0];
    } else if (mtheta <= 0) {
        mtheta12[0]=0;
        mtheta12[1]=0;
    } else {
        mtheta12[0]=mtheta-1;
        mtheta12[1]=mtheta;
    }

    // if (mtheta >= *(pthetaco+2) || mtheta <= 0) return(0.);

    //  cout << "mtheta : " << mtheta << endl;

    int mphi12[2];
    int mphi=(int)(- *pphico/ *(pphico+1) + phi / *(pphico+1) + 1);
    if (mphi >= *(pphico+2)) {
        mphi12[0]=(int)(*(pphico+2))-1;
        mphi12[1]=0;
    } else if (mphi <= 0) {
        mphi12[0]=0;
        mphi12[1]=1;
    } else {
        mphi12[0]=mphi-1;
        mphi12[1]=mphi;
    }

    //  cout << "mphi : " << mphi << endl;


    mr--;
    mtheta--;
    mphi--;

    float *paddress;

    int xxx,yyy,zzz;

    for (int ii=0;ii<2;ii++) {
        xxx=(mr+ii);
        for (int jj=0;jj<2;jj++) {
            yyy=mphi12[jj]*((int)(*(prco+2)));
            //paddress+=(mtheta+jj)* *;
            for (int kk=0;kk<2;kk++) {
                //zzz=(mtheta+kk)*((int)(*(prco+2)) * ((int)(*(pphico+2))));
                zzz=mtheta12[kk]*((int)(*(prco+2)) * ((int)(*(pphico+2))));

                //    cout << "address : " << xxx+yyy+zzz << endl;
                paddress=pdens+(xxx+yyy+zzz);
                //    cout << "val : " << *paddress << endl;

                neighborcube[ii][jj][kk]=*paddress;

            }
        }
    }


    float rint[2]={*prco + *(prco+1)* mr,*prco + *(prco+1) * (mr+1)};
    float phiint[2]={*pphico + *(pphico+1)*mphi,*pphico + *(pphico+1)*(mphi+1)};
    float thetaint[2]={*pthetaco + *(pthetaco+1)*mtheta,
                       *pthetaco + *(pthetaco+1)*(mtheta+1)};
    /*
    cout << "rint : " << rint[0] << " , " << rint[1] << endl;
    cout << "phiint : " << phiint[0] << " , " << phiint[1] << endl;
    cout << "thetaint : " << thetaint[0] << " , " << thetaint[1] << endl;
    */

    float neinterp=trilininterp(r,phi,theta,rint,phiint,thetaint,neighborcube);

    return(neinterp*1e10);

}


// --- no interpolation: just nearest neighbor
float densymingnn(const float x,const float y,const float z,const float phi0,
                  const float *prco,const float *pphico,const float *pthetaco,
                  float *pdens) {


    // rco : [rintercept,rslope,sizer]
    // phico : [phiintercept,phislope,sizephi]
    // thetaco : [thetaintercept,thetaslope,sizetheta]
    // dens : cube of density with [sizer,sizephi,sizetheta]

    //  float neighborcube[2][2][2];

    // -- convert x,y,z to r,theta,phi

    float r,phi,theta;

    cvcoord(x,y,z,&r,&phi,&theta);

    // -- find the nearest neighbors
    int mr=(int)(- *prco/ *(prco+1) + r / *(prco+1) + 1);
    if (mr >= *(prco+2) || mr <= 0)
        return 0.;

    int mtheta=(int)(- *pthetaco/ *(pthetaco+1) + theta / *(pthetaco+1) + 0.5);
    if (mtheta >= *(pthetaco+2)) {
        mtheta=(int)*(pthetaco+2)-1;
    } else if (mtheta <= 0) {
        mtheta=0;
    }

    int mphi=(int)(- *pphico/ *(pphico+1) + phi / *(pphico+1) + 0.5);
    if (mphi >= *(pphico+2)) {
        mphi=(int)(*(pphico+2))-1;
    } else if (mphi <= 0) {
        mphi=0;
    }


    float *paddress;
    int xxx,yyy,zzz;

    xxx=(mr);
    yyy=mphi*((int)(*(prco+2)));
    zzz=mtheta*((int)(*(prco+2)) * ((int)(*(pphico+2))));
    paddress=pdens+(xxx+yyy+zzz);

    float neinterp=*paddress;

    return(neinterp*1e10);

}

