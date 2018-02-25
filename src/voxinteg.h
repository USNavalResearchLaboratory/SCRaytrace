// $Id: voxinteg.h,v 1.2 2007-05-14 17:19:41 thernis Exp $

{
    // ---- compute the vertices pos on the image
    voxelvert2iijj(vij,voxel.vertobs,xstep,ystep,icntr,jcntr);

    // ---- find the extrema of pix pos in the image
    vector<float>::const_iterator it= min_element(vij[0].begin(),vij[0].end());
    float iminf=*it;
    int imin=(int) (*it +0.5); // smallest element in I
    it = min_element(vij[1].begin(),vij[1].end());
    float jminf=*it;
    int jmin=(int) (*it +0.5); // smallest element in J
    it = max_element(vij[0].begin(),vij[0].end());
    float imaxf=*it;
    int imax=(int) (*it +1); // largest element in I
    it = max_element(vij[1].begin(),vij[1].end());
    float jmaxf=*it;
    int jmax=(int) (*it +1); // largest element in J

    // -- check if not out of image bounds
    if (imin>=*is || jmin>=*js || imax<0 || jmax<0)
        continue;
    if (imin < 0) {
        imin=0;
        iminf=-0.5;
    }
    if (imax >= *is) {
        imax=*is-1;
        imaxf=float(imax)+0.5;
    }
    if (jmin < 0) {
        jmin=0;
        jminf=-0.5;
    }
    if (jmax >= *js) {
        jmax=*js-1;
        jmaxf=float(jmax)+0.5;
    }


    /*
    // ---- Test if voxel under resolution of 1 pix
    //      -0.5 for shift to the edge of the pix: 0,0 is the cntr of the pix
    if ( imaxf-floor(iminf+0.5)+0.5<1 && jmaxf-floor(jminf+0.5)+0.5<1) {
        countunder++;
        // ---- just use the vox center pos and total Ne
        int ip=int(iminf+0.5);
        int jp=int(jminf+0.5);

        if (ip < 0 || ip >= *is || jp < 0 || jp >= *js)
            continue;

        int impos=jp * *is + ip;

        // -- compute Ne
        netot[impos]+=voxel.dens*voxel.volume;

        // -- loop if only Ne is requested
        if (*flagneonly)
            continue;

        // -- impact parameter of the voxel center
        Cvec Vvoxcenterabs=obs.ui*voxel.voxcenterobs;
        float rho=psldist(obs.o,Vvoxcenterabs,abs.o);
        // -- skip integration if LOS within the occulter
        if (rho < *poccrad)
            continue;

        // -- dist vox cntr - sun cntr
        float r=(Vvoxcenterabs+obs.o).mag();

        // -- integration only in front of the Sun
        if (!(rho > LIMBLIMIT || (r > LIMBLIMIT && Vvoxcenterabs[2] < 0)))
            continue;
        float lengththru=0;
        voxel.GetThruLengthCenter(lengththru,0);

        // -- Compute Thomson scattering geometric factors and then the brightness
        float totterm,polterm;
        getThomsonGeomFactor(r,
                             rho,
                             totterm,polterm);

        btot[impos]+=constfactor*lengththru*totterm*voxel.dens;
        bpol[impos]+=constfactor*lengththru*polterm*voxel.dens;
        continue;
    }
    
    */
    
    /*
    // -- Test if voxel under resolution.
    //    If it falls into that case it's that the pix is in
    //    between two pix.
    //    Compute barycenter of vox projection and spread the density in a
    //    maxrespix boxcar.
    //      -0.5 for shift to the edge of the pix: 0,0 is the cntr of the pix
    //    for the moment use a 2x2 boxcar
    if ( false && imaxf-iminf<=1 && jmaxf-jminf<=1) {
        countunder2++;

        // ---- just use the vox center pos and total Ne
        // -- compute projection of the vox center on the image
        float vij[2];
        voxelvert2iijj(vij,voxel.voxcenterobs,xstep,ystep,icntr,jcntr);

        int ip=int(vij[0]);
        int jp=int(vij[1]);
        int impos=jp * *is + ip;
        if (impos <  maxnbpix) {
            float xw=(vij[0]-float(ip));
            float yw=(vij[1]-float(jp));

            // -- compute Ne
            netot[impos]+=voxel.dens*((1-xw)*(1-yw))*voxel.volume;
            if ((imin+1) < *is)
                netot[impos+1]+=voxel.dens*(xw*(1-yw))*voxel.volume;
            if (imin+1 < *is && jmin+1 < *js)
                netot[impos+1+ *is]+=voxel.dens*(xw*yw)*voxel.volume;
            if (jmin+1 < *js)
                netot[impos+ *is]+=voxel.dens*((1-xw)*yw)*voxel.volume;
        }

        // -- loop if only Ne is requested
        if (*flagneonly)
            continue;

        // -- impact parameter of the voxel center
        Cvec Vvoxcenterabs=obs.ui*voxel.voxcenterobs;
        float rho=psldist(obs.o,Vvoxcenterabs,abs.o);
        // -- skip integration if LOS within the occulter
        if (rho < *poccrad)
            continue;

        // -- dist vox cntr - sun cntr
        float r=(Vvoxcenterabs+obs.o).mag();

        // -- integration only in front of the Sun
        if (!(rho > LIMBLIMIT || (r > LIMBLIMIT && Vvoxcenterabs[2] < 0)))
            continue;
        float lengththru=0;
        voxel.GetThruLengthCenter(lengththru,0);

        // -- Compute Thomson scattering geometric factors and then the brightness
        float totterm,polterm;
        getThomsonGeomFactor(r,
                             rho,
                             totterm,polterm);

        if (impos <  maxnbpix) {
            float xw=(vij[0]-float(ip));
            float yw=(vij[1]-float(jp));

            // -- compute btot and bpol
            btot[impos]+=constfactor*lengththru*totterm*voxel.dens*((1-xw)*(1-yw));
            bpol[impos]+=constfactor*lengththru*polterm*voxel.dens*((1-xw)*(1-yw));
            rrr[impos]=lengththru;
            if ((imin+1) < *is) {
                btot[impos+1]+=constfactor*lengththru*totterm*voxel.dens*(xw*(1-yw));
                bpol[impos+1]+=constfactor*lengththru*polterm*voxel.dens*(xw*(1-yw));
            }
            if (imin+1 < *is && jmin+1 < *js) {
                btot[impos+1+ *is]+=constfactor*lengththru*totterm*voxel.dens*(xw*yw);
                bpol[impos+1+ *is]+=constfactor*lengththru*polterm*voxel.dens*(xw*yw);
            }
            if (jmin+1 < *js) {
                btot[impos+ *is]+=constfactor*lengththru*totterm*voxel.dens*((1-xw)*yw);
                bpol[impos+ *is]+=constfactor*lengththru*polterm*voxel.dens*((1-xw)*yw);
            }
        }
        continue;
    }
*/
    
    
    
    
    
    
    // -- dist obs to voxcenter
    float distobs2vox=voxel.voxcenterobs.mag();
    // -- surface area corresponding to the pixel square cone at the distance of the observer
    float voxpixedge=2*distobs2vox*tan(*fovpix/2); //*RSUN_CM
    float voxpixsurf=voxpixedge*voxpixedge;

    // ---- Scan the LOSes going through the voxel
    Cvec Vcntrvoxel;
    float lengththru;
    Cvec Vsuncntrvoxel=voxel.base.u*(voxel.base.o-abs.o); // pos of the sun cntr in the voxel coord sys
    countbig++;

    float totterm,polterm; // thomson scattering geometrical factors
    for (int jj=jmin;jj<=jmax;jj++) {
        int jpos=jj * *is;
        for (int ii=imin;ii<=imax;ii++) {
            int impos=ii + jpos;
            if (vlosimobs[impos].IsNull()) {
                vlosimobs[impos]=iijj2losobs((float) ii,(float) jj,xstart,ystart,xstep,ystep);
            }

            int outvoxelflag=voxel.GetThruLength(vlosimobs[impos],Vcntrvoxel,lengththru,0);

            // -- loop if out LOS out of the voxel
            if (outvoxelflag)
                continue;
            // -- compute Ne
            netot[impos]+=voxel.dens*lengththru*voxpixsurf;


            // -- loop if only Ne is requested
            if (*flagneonly)
                continue;

            if (vlosimabs[impos].IsNull()) {
                vlosimabs[impos]=obs.ui * vlosimobs[impos];
                rho[impos]=psldist(obs.o,vlosimabs[impos],abs.o);
            }
            // -- skip integration if LOS within the occulter
            if (rho[impos] < *poccrad)
                continue;

            // -- dist vox cntr - sun cntr
            float r=(Vsuncntrvoxel+Vcntrvoxel).mag();

            // -- integration only in front of the Sun
            if (!(rho[impos] > LIMBLIMIT || (r > LIMBLIMIT && vlosimabs[impos][2] < 0)))
                continue;

            // -- Compute Thomson scattering geometric factors and then the brightness
            getThomsonGeomFactor(r,
                                 rho[impos],
                                 totterm,polterm);

            btot[impos]+=constfactor*lengththru*totterm*voxel.dens;
            bpol[impos]+=constfactor*lengththru*polterm*voxel.dens;
        }
    }
}



/*
* $Log: voxinteg.h,v $
* Revision 1.2  2007-05-14 17:19:41  thernis
* Add CVS id and log in all files
*
*/
