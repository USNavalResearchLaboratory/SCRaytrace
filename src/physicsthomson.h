/** \file physicsthomson.h
 * \brief Implement Thomson scattering physics.
 */

#ifndef PHYSICSTHOMSON_H
# define PHYSICSTHOMSON_H

#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"

/** \brief Thomson scattering
 */
class PhysicsThomson : public PhysicsBase
{
protected:
    double u;                                               //!< limb darkening
    float constfactor;                                      //!< constant multiplicative factor of the Thomson scattering integral


public:
    PhysicsThomson() {
        physicsName = "Thomson Scattering";
        u = 0.58;                                           // Default limb darkening
        calcConstfactor();
      };

    ~PhysicsThomson() {};

    bool computeRadiation(const Cvec &vs, const float &r, const float &rho, float &btout, float &bpout, float &neout);

    void getConstFactors(float &btf,float &bpf,float &nef, float rho);

    void calcConstfactor() {constfactor = 1.24878E-25 / (1 - u/3) * RSUN_CM;};



    //! Compute the Thomson scattering coefficients for a given position in space
    //! \param r distance to the center of the Sun
    //! \param rho impact parameter
    //! \param btotcoeff returns the coefficient for the total brightness
    //! \param bpolcoeff returns the coefficient for the polarized brightness
    inline void getThomsonCoeff(const float &r,const float &rho,float &btotcoeff,float &bpolcoeff)
    {
        double sinomega = 1. / double(r);
        double sinsquareomega = sinomega * sinomega;
        double cossquareomega = 1 - sinsquareomega;
        double cosomega = sqrt(cossquareomega);

        double logterm = log((1. + sinomega) / cosomega);

        double a = cosomega * sinsquareomega;
        double b = -1./8. * (1. - 3. * sinsquareomega - cossquareomega * ((1. + 3. * sinsquareomega) / sinomega) * logterm);

        double c = (4./3.) - cosomega - (cosomega * cossquareomega) / 3.;
        double d = (1./8.) * (5. + sinsquareomega - cossquareomega * ((5. - sinsquareomega) / sinomega) * logterm);

        // ---- sum in the pixel
        double rhooverr = rho / r;
        // the polarized brightness
        bpolcoeff = float((a + u * (b - a)) * rhooverr * rhooverr);
        btotcoeff = float((2 * (c + u * (d - c)) - bpolcoeff));

    };

    float getLimbDarkening() {return u;};
    void setLimbDarkening(float limbdarkeningcoeff)
    {
        this->u=limbdarkeningcoeff;
        calcConstfactor();
    };

    
    void setParam(float *phyparam) {
        float limbdarkIn;
        limbdarkIn = (float)phyparam[0];
        this->setLimbDarkening(limbdarkIn);
    };
    
    void printParam() {
        printvar(this->u);
    };

};

#endif
