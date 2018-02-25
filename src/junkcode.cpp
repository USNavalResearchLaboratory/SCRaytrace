//
// C++ Implementation: junkcode
//
// Description: 
//
//
// Author:  <>, (C) 2011
//
// Copyright: See COPYING file that comes with this distribution
//
//



//! Thomson scattering integration along a LOS
    inline void losintegThomson(const unsigned int &i,const unsigned int &j)
    {

        Cvec vlosobs=camera.ij2los(float(i),float(j));
        Cvec vlosabs=obs.ui * vlosobs;

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

        unsigned int pos=i+j*camera.ccd.sxpix;

        for (unsigned int k=0;k<los.nbp;k++,vs+=vlosstep)
        {
            // -- dist Sun cntr - LOS current pos
            float r=vs.norm();

            // -- no integration within and behind the Sun
            if ( (r <= 1.001) || (rho <= 1.001 && (vs-obs.o).mag() > sqrt(obs.o.mag()*obs.o.mag()-rho*rho))) continue;

            // -- Call density here
            float ner=pmod->Density(ChangetoDensityCoord(modelposition,vs));

            if (ner <= 1e-1) continue;

            netot[pos]+=ner;
            if (neonly) continue;

            float btotcoeff,bpolcoeff;
            csun.getThomsonCoeff(r,rho,btotcoeff,bpolcoeff);

            btot[pos] +=ner*btotcoeff;
            bpol[pos] +=ner*bpolcoeff;

        }

        // multiply by the integral constant factor
        btot[pos] *=csun.constfactor*los.ds;
        bpol[pos] *=csun.constfactor*los.ds;
        netot[pos] *=RSUN_CM*los.ds;

    };

