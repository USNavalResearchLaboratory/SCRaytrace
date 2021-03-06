/** \file physicsuv.h
 * \brief Implements UV physics.
 */

#ifndef PHYSICSUV_H
#define PHYSICSUV_H

#include "config.h"
#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"

/** \brief UV emission
*/
class PhysicsUV : public PhysicsBase
{
public:
    PhysicsUV();

    bool IsGood();
    float getyisel(const unsigned int &i,const unsigned int &j);
    float calcEmissivity(const unsigned int &kline,const float &te);

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btout,float &bpout,float &neout);
    void getConstFactors(float &btf,float &bpf,float &nef, float rho);

    //! Set the parameter of the physics model
    void setParam(float *phyparam)
    {
        wavelengthId=(unsigned int)phyparam[0];
    };

    void printParam()
    {
        printvar(wavelengthId);
    }
    
private:
    static const unsigned int NBSAMP=51;
    static const unsigned int NBLINES=4;
    float yisel[NBSAMP*NBLINES];
    float emis[NBSAMP*NBLINES];
    float ti[NBSAMP];
    float teg[NBSAMP];
    bool isgood;

    unsigned int wavelengthId;

};

#endif
