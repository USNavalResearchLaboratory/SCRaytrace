/*! \file camera.h 
 * \brief Define image resolution, camera projection type
 *
 *  
 */


#ifndef CAMERA_H
#define CAMERA_H

#include <wcslib/wcslib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <map>
#include <stdexcept>


#include "Cvec.h"
#include "constant.h"


//! Projection types
enum ProjType {ARC, TAN, SIN, AZP, ZPN};

class Camera;
class Scene;

//! Defines the Detector characteristics of the Camera
class Detector{
public:
    Detector();
    Detector(const unsigned int &sxpix0,
            const unsigned int &sypix0,
            const float &sxmm0,
            const float &symm0);
    Detector(const unsigned int &sxpix0,
            const unsigned int &sypix0);
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
    // #define NPV_MAX 64
    static const int NPV_MAX=64;

    // const int NPV_MAX;
    Camera() {
        fovpix=0.1;
        projtype=ARC;
        detector=Detector();
        crpix1=64;
        crpix2=64;
        pv2_1=0;
        pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;
        setDefaultPV();
        wcs = NULL;
        setWCS();
    };
    Camera(const float &fovpix0,
        const ProjType &projtype0,
        const Detector &detector0) {
        fovpix=fovpix0;
        projtype=projtype0;
        detector=detector0;
        crpix1=0;
        crpix2=0;
        pv2_1=0;
        setDefaultPV();
        pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;
        wcs = NULL;
        setWCS();
    };
    Camera(const float &fovpix0,
        const ProjType &projtype0,
        const Detector &detector0,
        const float &crpix10,
        const float &crpix20,
        const float &pv2_10)
    {
        fovpix=fovpix0;
        projtype=projtype0;
        detector=detector0;
        crpix1=crpix10;
        crpix2=crpix20;
        pv2_1=pv2_10;
        pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;
        setDefaultPV();
        wcs = NULL;
        setWCS();
    };
    Camera(const float &fovpix0,
        const ProjType &projtype0,
        const Detector &detector0,
        const float &crpix10,
        const float &crpix20,
        const int &npv0,
        float *pv0,
        int *pv_i0,
        int *pv_m0)
    {
        fovpix=fovpix0;
        projtype=projtype0;
        detector=detector0;
        crpix1=crpix10;
        crpix2=crpix20;
        pc[0]=1;pc[1]=0;pc[2]=0;pc[3]=1;
        pvInit(npv0, pv0, pv_i0, pv_m0);
        wcs = NULL;
        setWCS();
    };
    ~Camera()
    {
        if (wcs) delete wcs;
    };


    //! Map projection type to WCS standard type
    std::map<ProjType, std::string> projtype2wcscode = {
        { ARC, "ARC" },
        { TAN, "TAN" },
        { SIN, "SIN" },
        { AZP, "AZP" },
        { ZPN, "ZPN" }
    };

private:
    //! Set the WCS projection parameters, for wcslib.
    void setWCS() 
    {
        if (wcs) {
            delete wcs;
        }
        wcs = new wcsprm;
        wcs->flag = -1;
        const int NAXIS = 2;
        wcsini(1, NAXIS, wcs);
        wcs->crpix[0] = crpix1;
        wcs->crpix[1] = crpix2;
        double *pcij;
        pcij = wcs->pc;
        *(pcij++) = pc[0]; // [0, 0] = [i, j]
        *(pcij++) = pc[1]; // [0, 1]
        *(pcij++) = pc[2]; // [1, 0]
        *(pcij++) = pc[3]; // [1, 1]
        wcs->cdelt[0] = fovpix * RADEG;
        wcs->cdelt[1] = fovpix * RADEG;

        std::string ctypelon, ctypelat;
        ctypelon.append("HPLN-");
        ctypelon.append(projtype2wcscode[projtype]);
        ctypelat.append("HPLT-");
        ctypelat.append(projtype2wcscode[projtype]);

        strcpy(wcs->ctype[0], ctypelon.c_str());
        strcpy(wcs->ctype[1], ctypelat.c_str());

        wcs->crval[0] = 0;
        wcs->crval[1] = 0;
        wcs->lonpole = 180;
        wcs->latpole = 0;
        prj.pv[1] = pv2_1;


        struct pvcard PV[NPV_MAX];
        if (pv2_1 != 0 && projtype == AZP) {
            wcs->npv = npv;

            PV[0].i = 1;
            PV[0].m = 1;
            PV[0].value = 0.0;
            wcs->pv[0] = PV[0];

            PV[1].i = 1;
            PV[1].m = 2;
            PV[1].value = 90.0;
            wcs->pv[1] = PV[1];

            PV[2].i = 1;
            PV[2].m = 3;
            PV[2].value = 180.0;
            wcs->pv[2] = PV[2];

            PV[3].i = 2;
            PV[3].m = 1;
            PV[3].value = pv2_1;
            wcs->pv[3] = PV[3];

        } else {
            wcs->npv = npv;
            int i;
            for (i = 0; i < npv; i++) {
                PV[i].i = pv_i[i];
                PV[i].m = pv_m[i];
                PV[i].value = pv[i];
                wcs->pv[i] = PV[i];
            }
        }

        if (wcsset(wcs)) {
            wcsperr(wcs, "Error");
            std::cout << "WCS no good" << std::endl;
        } else {
            std::cout << "WCS ok" << std::endl;
        }
    };

public:
    //! Print the WCS projection parameters, from wcslib.
    void printWCS() {wcsprt(wcs);};

    void setFovpix(const float &fovpix0) {
        fovpix=fovpix0;
        setWCS();
    };
    float getFovpix() const {return fovpix;};
    
    void setCrpix(const float &crpix10,
                    const float &crpix20) {
        crpix1=crpix10;crpix2=crpix20;
        setWCS();
    };
    float getCrpix1() const {return crpix1;};
    float getCrpix2() const {return crpix2;};
    
    void setPv2_1(const float &pv2_10) {
        pv2_1=pv2_10;
        setWCS();
    };
    float getPv2_1() const {return pv2_1;};
    
    void setPc(float pc0[4]) {
        pc[0]=pc0[0];
        pc[1]=pc0[1];
        pc[2]=pc0[2];
        pc[3]=pc0[3];
        setWCS();
    };
    float getPc(const int &id) const {return pc[id];};
    
    void setPv(const int &npv0,
                float *pv0,
                int *pv_i0,
                int *pv_m0) {
        pvInit(npv0, pv0, pv_i0, pv_m0);
        setWCS();
    }

    void setProjType(const ProjType &projtype0) {
        projtype=projtype0;
        setWCS();
    };
    void setProjType(const unsigned int &projtypecode) 
    {
        switch (projtypecode)
        {
            case 1 : projtype = ARC; break;
            case 2 : projtype = TAN; break;
            case 3 : projtype = SIN; break;
            case 4 : projtype = AZP; break;
            case 5 : projtype = ZPN; break;
            default : projtype = ARC;
        }
        setWCS();
    };
    ProjType getProjType() const {return projtype;};
    
    void setDetector(const Detector &detector0) {
        detector=detector0;
        setWCS();
    };
    Detector getDetector() const {return detector;};

    //! Compute the LOS of a given pixel
    inline Cvec ij2los(const float &i,const float &j) {
        return ij2loswcs(i, j);
    }


    //! Compute the LOS of a given pixel
    inline Cvec ij2losold(const float &i,const float &j)
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
    
    //! Compute the LOS for a given pixel, using wcslib
    Cvec ij2loswcs(const float &i, const float &j)
    {
        // ---- see wcslib
        const int NCOORD=1;
        double pixel[NCOORD][2];
        double imgcrd[NCOORD][2];
        double phi[NCOORD];
        double theta[NCOORD];
        double world[NCOORD][2];
        int stat[NCOORD];

        pixel[0][0] = i;
        pixel[0][1] = j;

        if (wcsp2s(wcs, 
                NCOORD,
                2, 
                pixel[0], 
                imgcrd[0], 
                phi, 
                theta, 
                world[0], 
                stat)) {
            wcsperr(wcs, "");
        }

        // wcsprt(wcs);

        // std::cout << "stat : " << stat[0] << std::endl;
        // std::cout << "pixel" << std::endl;
        // std::cout << pixel[0][0] << ", " << pixel[0][1] << std::endl;
        // std::cout << "phi, theta" << std::endl;
        // std::cout << phi[0] << ", " << theta[0] << std::endl;
        // std::cout << "world" << std::endl;
        // std::cout << world[0][0] << ", " << world[0][1] << std::endl;

        float rrr, alpha;
        rrr = (90 - theta[0]) / RADEG;
        alpha = (phi[0] - 90) / RADEG;
        // std::cout << "theta: " << theta[0] << std::endl;
        // std::cout << "phi  : " << phi[0] << std::endl;
        return Cvec(sin(rrr)*sin(alpha),
                    sin(rrr)*cos(alpha),
                    cos(rrr));
    };


    friend class Scene;
    
protected:
    Detector detector;      //!< Detector characteristics
    float fovpix;           //!< fovpix in rad/pix
    float crpix1,crpix2;    //!< Center of the optical axis
    ProjType projtype;      //!< projection type
    float pv2_1;            //!< projection extra parameter, useful for AZP
    float pc[4];            //!< rotation matrix
    struct prjprm prj;      //!< WCSLIB handler for computation of the projections
    int npv;                //!< number of pv values
    int pv_i[NPV_MAX];      //!< pv array containing the i index. Size: [npv]
    int pv_m[NPV_MAX];      //!< pv array containing the m index. Size: [npv] 
    float pv[NPV_MAX];      //!< pv array with the values. Size: [npv]
    wcsprm *wcs;            //!< wcslib wcs handler

private:
    void setDefaultPV()
    {  
        npv = 3;
        pv[0] = 0.;pv[1] = 90.;pv[2] = 180.;
        pv_i[0] = 1;pv_i[1] = 1;pv_i[2] = 1;
        pv_m[0] = 1;pv_m[1] = 2;pv_m[2] = 3;
    }

    void pvInit(const int &npv0,
                float *pv0,
                int *pv_i0,
                int *pv_m0) {
        npv=npv0;
        if (npv > NPV_MAX) {
            std::cout << "PV is too large." << std::endl;
            std::string mess = "npv out of range: max ";
            mess.append("NPV_MAX");
            throw std::out_of_range(mess);
        }

        int i;
        for (i = 0;i < npv; i++){
            pv[i] = pv0[i];
            pv_i[i] = pv_i0[i];
            pv_m[i] = pv_m0[i];
        }
    };


};

#endif
