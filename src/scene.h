/** \file scene.h
 * \brief Defines the scene: Camera + Sun + Density
 */

#ifndef SCENE_H
#define SCENE_H

#include <string>
#include "camera.h"
// #include "sun.h"
#include "Cbasis.h"
#include "Cvec.h"
#include "Clos.h"
#include "CModelBase.h"
#include "rtmiscfunc.h"
#include "ModelPosition.h"
#include "physicsbase.h"
#include <cstdio>

//! Defines the Scene: Camera + Sun + Density
class Scene
{
  
  protected:
    Cbasis abs,obs; //!< absolute and observer coordinate system
    CModelBase *pmod; //!< points to the density model
    PhysicsBase *pphy; //!< points to the model of physics
    float *btot,*bpol,*netot; //!< points to the total brightness, polarized brightness and total density images
    float fracmax;      //!< fraction of the total brightness for computation of the depth of the contribution
    float *pintegrand;  //!< save all the integration points in this structure if requested
    
  private:
    bool *pisrunning;
    bool neonly,quiet,frontinteg;
    unsigned int nbpix,chunksize,lastchunkremain,nbchunk;
    bool runfracmax;    //!< flag that is used to compute or not the depth contribution
    bool runDumpInteg;  //!< set this flag for saving the integration points in the pintegand array
    
  public:
    Camera camera;
//     Sun csun;
    Clos los;
    ModelPosition modelposition; //!< model positioning

    friend class PhysicsBase;
    friend class PhysicsThomson;
    friend class PhysicsUV;
    friend class PhysicsMie;
    friend class PhysicsIsotropic;
    friend class PhysicsVSF;
    friend class PhysicsVSFVaryDist;
    friend class PhysicsVariableVSF;

    friend class CModelBase;
    friend class CModel73;
    friend class CModel80;
    friend class CModel81;
    friend class CModel82;
    friend class CModel83;
    friend class CModel84;
    
    Scene();
    
    ~Scene()
    {
        delete pmod;
        delete pphy;
    }
    
    
    //! Compute the image with multi-threading for each pixel
    void computeImagebyRay(float *btot,float *bpol,float *netot,const unsigned int nbthread);
    
    //! Compute the image with multi-threading by chunk of the image
    void computeImagebyChunk(float *btot,float *bpol,float *netot,const unsigned int nbthread,const unsigned int nbchunk);
    
    
    void setabs(Cbasis abs)
    {
        this->abs=abs;
    };
    
    
    void setobs(Cbasis obs)
    {
        this->obs=obs;
        if(quiet != 1) std::cout << "obs : " << obs << std::endl;
    };

    
    //! Set the density model to be used
    //! \param modelid the id of the model
    //! \param *pmodparam points to the vector of parameters
    void setDensityModel(int modelid,float *pmodparam)
    {
        delete pmod;
        pmod = modelselect(modelid);
        pmod->initParam(pmodparam);
        pmod->setParentScene(this);
    };

    
    void setNeonly(int neonly)
    {
        this->neonly=(bool)neonly;
    };
    
    
    void setQuiet(int quiet)
    {
        this->quiet=(bool)quiet;
    };
    
    
    void setFrontInteg(int frontinteg)
    {
        this->frontinteg=(bool)frontinteg;
    };
    
    
    void setPhysics(PhysicsType phytype);
    
    
    string getPhysics()
    {
        return pphy->getPhysics();
    };
    
    
    void setPhysicsParam(float *phyparam) {pphy->setParam(phyparam);};
    void printPhysicsParam() {pphy->printParam();};
    void setFracMax(float fracmax);
    void setRunFracMaxOn()  { runfracmax = true ;};
    void setRunFracMaxOff() { runfracmax = false;};
    void setRunDumpIntegOn(float *pintegrand) { this->pintegrand = pintegrand;
                                                runDumpInteg = true;};
                                                
    void setRunDumpIntegOff() { runDumpInteg = false;};
    void setDumpIntegrand(bool flagrunDumpInteg, float *pintegrand);

private:

//! Generic LOS integration: the physics has to be initialized beforehand
    inline void losinteg(const unsigned int &i,const unsigned int &j)
    {
      
        Cvec vlosobs=camera.ij2los(float(i),float(j));
        Cvec vlosabs=obs.ui * vlosobs;
        
        float btout,bpout,neout;
        bool flagnull=0;

        // -- compute rho : dist LOS - Sun cntr: impact parameter
        float rho=psldist(obs.o,vlosabs,abs.o);

        // -- origin of the LOS: closest distance to sun by default
        Cvec qlos;
        if (frontinteg)
        {
            qlos=obs.o;
        }
        else
        {
            qlos=orthoproj(obs.o,vlosabs,abs.o);
        }

        // -- Pass the LOS unit vector to the model. This is needed for some models.
//         pphy->initDensityModel(vlosabs);
        
        // -- init LOS integration parameters
        Cvec vlosstep=vlosabs * los.ds;
        Cvec vs=qlos + vlosabs * los.sstart;
        
        unsigned int pos=i+j*camera.detector.sxpix;

        // if(i < 130 && i > 125 && j < 130 && j > 125){ //(DEBUG)
        //     printf("Pre: (%d, %d): BTOT: %f; WLOS: (%f, %f, %f); VS: (%f, %f, %f); obs.o: (%f, %f, %f)\n", i, j, btot[pos]*pow(10, 10), vlosabs[0], vlosabs[1], vlosabs[2], vs[0], vs[1], vs[2], obs.o[0], obs.o[1], obs.o[2]); //(DEBUG)
        // } //(DEBUG)

        for (unsigned int k=0;k<los.nbp;k++,vs+=vlosstep)
        {
            // -- dist Sun cntr - LOS current pos
            float r = vs.norm();

            // -- no integration within and behind the Sun
            // if ( (r <= 1.001) || (rho <= 1.001 && (vs-obs.o).mag() > sqrt(obs.o.mag()*obs.o.mag()-rho*rho))) continue;

            // -- Compute radiation of the volume element depending on chosen physics
            //    Also call the density in here.
            flagnull = pphy->computeRadiation(vs,r,rho,btout,bpout,neout);

            if (flagnull) continue;

            netot[pos] += neout;
            if (neonly) continue;

            btot[pos] += btout;
            bpol[pos] += bpout;

            // if(i == 128 && j == 128 && btout != 0){ //(DEBUG)
                // printf("Intensity = (%f) at (%d, %d) with VS: (%f, %f, %f)\n", btout*pow(10, 10), i, j, vs[0], vs[1], vs[2]); //(DEBUG)
            // } //(DEBUG)
        }

        // ---- multiply by the integral constant factor, depending on the physics
        float btf, bpf, nef;
        pphy->getConstFactors(btf, bpf, nef, rho);
        btot[pos]  *= btf;
        bpol[pos]  *= bpf;
        netot[pos] *= nef;
        // printf("Post: (%d, %d): BTOT: %f; WLOS: (%f, %f, %f); VS: (%f, %f, %f); obs.o: (%f, %f, %f)\n", i, j, btot[pos]*pow(10, 10), vlosabs[0], vlosabs[1], vlosabs[2], vs[0], vs[1], vs[2], obs.o[0], obs.o[1], obs.o[2]); //(DEBUG)
    }

    //! Compute distance to the fraction of total brightness
    inline void losintegtofrac(const unsigned int &i,const unsigned int &j)
    {
      
        Cvec vlosobs=camera.ij2los(float(i),float(j));
        Cvec vlosabs=obs.ui * vlosobs;
        
        float btout,bpout,neout;
        bool flagnull=0;

        // -- compute rho : dist LOS - Sun cntr: impact parameter
        float rho=psldist(obs.o,vlosabs,abs.o);

        // -- origin of the LOS: observer by default
        Cvec qlos;
        if (frontinteg)
        {
            qlos=obs.o;
        }
        else
        {
            qlos=orthoproj(obs.o,vlosabs,abs.o);
        }

        // -- init LOS integration parameters
        Cvec vlosstep=vlosabs * los.ds;
        Cvec vs=qlos + vlosabs * los.sstart;
        
        unsigned int pos=i+j*camera.detector.sxpix;

        // -- get the integral constant factor
        float btf,bpf,nef;
        pphy->getConstFactors(btf, bpf, nef, rho);
        
        // -- define a bucket where we will add up the contribution up until the requested fraction of the total
        float bucket = 0.;
        float fracBtot = btot[pos] * fracmax;
        bpol[pos] = 0.;
        
        for (unsigned int k=0; k<los.nbp; k++,vs+=vlosstep)
        {
            // -- dist Sun cntr - LOS current pos
            float r=vs.norm();

            // -- no integration within and behind the Sun
            // if ( (r <= 1.001) || (rho <= 1.001 && (vs-obs.o).mag() > sqrt(obs.o.mag()*obs.o.mag()-rho*rho))) continue;

            // -- Compute radiation of the volume element depending on chosen physics
            //    Also call the density in here.
            flagnull=pphy->computeRadiation(vs,r,rho,btout,bpout,neout);

            if (flagnull) continue;

            // -- compute only total brightness
            bucket += btout * btf;
            
            // -- and add up the integration distance in the polarized brigthness image
            bpol[pos] += los.ds;

            // -- compare to the sum of total brightness and exit if above the threshold
            if (bucket >= fracBtot) break;
            
        }

    }
    
    
    //! Same as losinteg but save all the integration points in an array
    inline void losintegDumpIntegrand(const unsigned int &i,const unsigned int &j)
    {
      
        Cvec vlosobs=camera.ij2los(float(i),float(j));
        Cvec vlosabs=obs.ui * vlosobs;
        
        float btout,bpout,neout;
        bool flagnull=0;

        // -- compute rho : dist LOS - Sun cntr: impact parameter
        float rho=psldist(obs.o,vlosabs,abs.o);

        // -- get integrals constant factors
        float btf,bpf,nef;
        pphy->getConstFactors(btf, bpf, nef, rho);

        // -- origin of the LOS: observer by default
        Cvec qlos;
        if (frontinteg)
        {
            qlos=obs.o;
        }
        else
        {
            qlos=orthoproj(obs.o,vlosabs,abs.o);
        }

        // -- init LOS integration parameters
        Cvec vlosstep=vlosabs * los.ds;
        Cvec vs=qlos + vlosabs * los.sstart;
        
        // -- compute pixel position in lexicographic ordered output arrays
        unsigned int pos=i+j*camera.detector.sxpix;

        // -- init integral point index for the integrand output array
        unsigned long posIntegrand;
                
        for (unsigned int k=0;k<los.nbp;k++,vs+=vlosstep)
        {
            // -- dist Sun cntr - LOS current pos
            float r=vs.norm();

            // -- no integration within and behind the Sun
            // if ( (r <= 1.001) || (rho <= 1.001 && (vs-obs.o).mag() > sqrt(obs.o.mag()*obs.o.mag()-rho*rho))) continue;

            // -- Compute radiation of the volume element depending on chosen physics
            //    Also call the density in here.
            flagnull=pphy->computeRadiation(vs,r,rho,btout,bpout,neout);
            
            // -- compute index position in the integrand array
            posIntegrand = pos + camera.detector.sxpix * camera.detector.sypix * k;
            
            pintegrand[posIntegrand] = btout * btf;
//             pintegrand[posIntegrand] = bpout;
            
            if (flagnull) continue;

            netot[pos]+=neout;
            if (neonly) continue;

            btot[pos] +=btout;
            bpol[pos] +=bpout;
        }

        // ---- multiply by the integral constant factor, depending on the physics
        btot[pos]  *=btf;
        bpol[pos]  *=bpf;
        netot[pos] *=nef;
    }

    
    
    //! Multi-threaded Thomson scattering integration along a LOS, for the "by Ray" integration
    void losintegray(const unsigned int &i,const unsigned int &j,const unsigned int &threadid);
    //! Multi-threaded Thomson scattering integration along a LOS, for the "by Chunk" integration
    void losintegchunk(const unsigned int &i,const unsigned int &threadid);

};





#endif
