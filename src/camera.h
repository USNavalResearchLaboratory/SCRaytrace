/*! \file camera.h 
 * \brief Define image resolution, camera projection type
 *
 *  
 */


#ifndef CAMERA_H
#define CAMERA_H

#include "Cvec.h"
#include "constant.h"


//! Projection types
enum ProjType {ARC, TAN, SIN, AZP};

class Camera;
class Scene;

//! Defines the Detector characteristics of the Camera
class Detector{
public:
    Detector();
    Detector(const unsigned int &sxpix0,const unsigned int &sypix0,const float &sxmm0,const float &symm0);
    Detector(const unsigned int &sxpix0,const unsigned int &sypix0);
    Detector(const Detector &);
    const Detector &operator=(const Detector &);
    bool operator==(const Detector &) const;
    
    void setSizePix(const unsigned int &sxpix0,const unsigned int &sypix0);
    unsigned int getSizePixX() const;
    unsigned int getSizePixY() const;
    
    void setSizemm(const float &sxmm0,const float &symm0);
    float getSizemmX() const;
    float getSizemmY() const;
    
    friend class Camera;
    friend class Scene;
    
protected:
    unsigned int sxpix,sypix; //!< size of the Detector in pix
    float sxmm,symm; //!< size of the Detector in mm
};



//! Apply WCS projection to the rrr angle
inline float applyprojection(const float &rrr,const ProjType &projtype,const float &pv2_1)
{
    float rout=rrr;
    switch (projtype) {
        case ARC : break;
        case TAN : rout=tan(rrr);break;
        case SIN : rout=sin(rrr);break;
        case AZP : rout=((pv2_1+1)*sin(rrr))/(pv2_1+cos(rrr));
    }
    return rout;
}

//! Apply inverse WCS projection to the rrr angle
inline float applyinverseprojection(const float &rrr,const ProjType &projtype,const float &pv2_1)
{
    float rout=rrr;
    switch (projtype) {
        case ARC : break;
        case TAN : rout=atan(rrr);break;
        case SIN : rout=asin(rrr);break;
        case AZP : 
            float r=rrr/(pv2_1+1);
            rout=(PI/2.)-(atan2(float(1.),float(r))-asin(r*pv2_1/sqrt(r*r+1.)));
    }
    return rout;
}



//! Defines the instrument optical characteristics
class Camera{
public:
    Camera() {fovpix=0.;projtype=ARC;detector=Detector();crpix1=0;crpix2=0;pv2_1=0;pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;};
    Camera(const float &fovpix0,const ProjType &projtype0,const Detector &detector0)
    {fovpix=fovpix0;projtype=projtype0;detector=detector0;crpix1=0;crpix2=0;pv2_1=0;pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;};
    Camera(const float &fovpix0,const ProjType &projtype0,const Detector &detector0,const float &crpix10,const float &crpix20,const float &pv2_10)
    {fovpix=fovpix0;projtype=projtype0;detector=detector0;crpix1=crpix10;crpix2=crpix20;pv2_1=pv2_10;pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;};
    
    void setFovpix(const float &fovpix0) {fovpix=fovpix0;};
    float getFovpix() const {return fovpix;};
    
    void setCrpix(const float &crpix10,const float &crpix20) {crpix1=crpix10;crpix2=crpix20;};
    float getCrpix1() const {return crpix1;};
    float getCrpix2() const {return crpix2;};
    
    void setPv2_1(const float &pv2_10) {pv2_1=pv2_10;};
    float getPv2_1() const {return pv2_1;};
    
    void setPc(float pc0[4]) {pc[0]=pc0[0];pc[1]=pc0[1];pc[2]=pc0[2];pc[3]=pc0[3];};
    float getPc(const int &id) const {return pc[id];};
    
    void setProjType(const ProjType &projtype0) {projtype=projtype0;};
    void setProjType(const unsigned int &projtypecode) 
    {
        switch (projtypecode)
        {
            case 1 : projtype = ARC; break;
            case 2 : projtype = TAN; break;
            case 3 : projtype = SIN; break;
            case 4 : projtype = AZP; break;
            default : projtype = ARC;
        }
    };
    ProjType getProjType() const {return projtype;};
    
    void setDetector(const Detector &detector0) {detector=detector0;};
    Detector getDetector() const {return detector;};
    
    //! Compute the LOS of a given pixel
    inline Cvec ij2los(const float &i,const float &j)
    {
        // ---- see Calabretta & Greisen, A&A, 2002
        float ic=i-crpix1;
        float jc=j-crpix2;
        
        
        float x=pc[0]*ic+pc[1]*jc; 
        float y=pc[2]*ic+pc[3]*jc;
        
        
        float rrr=fovpix*sqrt(x*x+y*y);
        float alpha=atan2(y,x);
        
        rrr=applyinverseprojection(rrr,projtype,pv2_1);
        
        return Cvec(sin(rrr)*sin(alpha),sin(rrr)*cos(alpha),cos(rrr));
    };
    
    friend class Scene;
    
protected:
    Detector detector;      //!< Detector characteristics
    float fovpix;           //!< fovpix in rad/pix
    float crpix1,crpix2;    //!< Center of the optical axis
    ProjType projtype;      //!< projection type
    float pv2_1;            //!< projection extra parameter, useful for AZP
    float pc[4];            //!< rotation matrix
};

#endif
