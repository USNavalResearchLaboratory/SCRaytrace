// $Id$

#ifndef PHYSICSMIE_H
#define PHYSICSMIE_H

#include <fstream>
#include "physicsbase.h"
#include "Cvec.h"
#include "constant.h"
#include "miehps.h"

/**
Line of sight integration for Mie scattering
*/
class PhysicsMie : public PhysicsBase
{
public:
    PhysicsMie() {physicsName="Mie Scattering";};

    ~PhysicsMie() { delete mie;};

    bool computeRadiation(const Cvec &vs,const float &r,const float &rho,float &btot,float &bpol,float &neout);

    void getConstFactors(float &btf,float &bpf,float &nef);

    
    //! Set the parameter of the physics model
    void setParam(float *phyparam)
    {
        particleRadius = (float)phyparam[0];
        wavelength = (float)phyparam[1];
        refractiveIndexRealPart = (float)phyparam[2];
        refractiveIndexImaginaryPart = (float)phyparam[3];
        
        // -- wave number
        k = 2 * PI / wavelength;
        
        // -- init mie scattering phase functions
        x_sizeparam = 2 * PI * particleRadius / wavelength;
        m = complex<double>(refractiveIndexRealPart,refractiveIndexImaginaryPart);

        cout << "x parameter : " << x_sizeparam << endl;
        
        mie = new Mie(x_sizeparam,m);
        mie->calc();
        
        // -- compute the lookup tables
        float angle;
        for (unsigned int i=0;i<SIZELOOKUPTABLE;i++)
        {
          angle = float(i) * PI / float(SIZELOOKUPTABLE-1);
          mie->calcS(angle);
          
          // -- See "Light Scattering by Small Particles", H.C. Van de Hulst, Dover Publication N.Y., 1957, 1981, ISBN 0-486-64228-3, p. 35
          // I_perp = \frac{i_1}{k^2 r^2} I_0  : irradiance for the perpendicular polarization
          // I_para = \frac{i_2}{k^2 r^2} I_0  : irradiance for the parallel polarization
          // I = \frac{\frac{1}{2}(i_1 + i_2)}{k^2 r^2} I_0 : natural light
          // Q = \frac{\frac{1}{2}(i_2 - i_1)}{k^2 r^2} I_0 : polarized light (see Stokes parameters I Q U V)
          i1table[i] = pow(abs(mie->getS1()),2);
          i2table[i] = pow(abs(mie->getS2()),2);
          
      
        }
        
        // -- Save tables in a file
        fstream file;
        file.open("i1table.dat",ios_base::out | ios_base::binary);
        file.write((char*)i1table, sizeof(float) * SIZELOOKUPTABLE);
        file.close();
        
        file.open("i2table.dat",ios_base::out | ios_base::binary);
        file.write((char*)i2table, sizeof(float) * SIZELOOKUPTABLE);
        file.close();
        
        
    };

    void printParam()
    {
        printvar(particleRadius);
        printvar(wavelength);
        printvar(m);
   }

  private:
    float particleRadius;   // Radius of the spherical particle.
    float wavelength;       // Wavelength, in the same units than the particle radius.
    float refractiveIndexRealPart;
    float refractiveIndexImaginaryPart;
    float x_sizeparam;      // Size parameter defined as 2 * pi * particleRadius / wavelength.
    complex<double> m;      // Index of refraction of the particule.
    Mie *mie;               // The Mie scattering class, using miehps implementation.
    float k;                // wave number
    
    complex<double> s1;
    complex<double> s2;
        
    static const unsigned int SIZELOOKUPTABLE=36001;  // Size of the lookup table
    
    float i1table[SIZELOOKUPTABLE];
    float i2table[SIZELOOKUPTABLE];

    
};

#endif
