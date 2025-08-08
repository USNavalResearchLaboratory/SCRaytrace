"""
scraytrace
==========

.. module:: scraytrace
   :synopsis: Python interface for the solar corona ray-tracing tools.

Python interface for the solar corona ray-tracing tools.

.. todo::
   * Dynamically find the compiled libraries based on the architechture
   
"""


import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import wcs
from matplotlib import pyplot as plt
from astropy.constants import iau2015
import pandas as pd

import mathutil as mu
from fcor_kl import fcor_kl

import scraytrace as scrt


if __name__ == "__main__" :

    plt.close('all')

#    runname = 'COR2'
#    runname = 'wispr'
#    runname = 'secchi'
#    runname = 'disk'
    runname = 'saito'


    if runname == 'COR2':

        # COR2 image
        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/cor2/20120802/20120802_125400_d4c2A.fts')
        minmaxfits = (2e3, 1e4)

        mu.dispFits(hdu[0].data, hdu[0].header, log=True, minmax=minmaxfits, cmap='flag')
        
        w = wcs.WCS(hdu[0].header)
        
            
        sim1 = scrt.scraytrace(losnbp=100)
        sim1.setParamFromWCS(hdu[0].header, [256, 256])
        sim1.raytrace()
    #    sim1.dispim()
        
        sim1.printRTParam()
        sim1.createWCSHead()
        print(sim1.wcsrt)
    
        for k in sim1.fitsrt: print(k, sim1.fitsrt[k])
        
        print(w)
        
        print(wcs.WCS(sim1.fitsrt))    
        
        
        mu.dispFits(sim1.btot, sim1.fitsrt, log=True, minmax=(5e-15, 2e-11), units='B0', cmap='flag')




    elif runname == 'wispr':

        modelid = 78
    
        minmaxmod = {}
        # Read and display WISPR fits image
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20181101/psp_L2_wispr_20181101T231550_1221.fits')
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190326/psp_L2_wispr_20190326T155032_1221.fits')
        minmaxfits = (1e-13, 5e-12)
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190323/psp_L2_wispr_20190323T155032_1221.fits')
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190405/psp_L2_wispr_20190405T150730_1221.fits')
#        minmaxfits = (1e-12, 5e-11)
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190405/psp_L2_wispr_20190405T152258_2222.fits')
#        minmaxfits = (5e-13, 2e-11)
        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L1/20190405/psp_L1_wispr_20190405T234142_1221.fits')
        minmaxfits = (1e3, 3e4)

        minmaxmod[74] = (1e-12, 7e-11)
        minmaxmod[78] = (1e-15, 2e-13)

        # COR2 A
#        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/cor2/20120802/20120802_125400_d4c2A.fts')
#        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/cor2/20120530/20120530_235400_d4c2A.fts')
#        minmaxfits = (2e3, 1e4)

        # HI1 A
#        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/hi_1/20181101/20181101_124901_s4h1A.fts')
#        minmaxfits = (2e3, 1e6)
        
        # HI2 A
#        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/hi_2/20181101/20181101_120921_s4h2A.fts')
#        minmaxfits = (2e3, 1e6)
#        mimmaxmod78 = (4e-17, 4e-15)
        
#        hdu[0].header['HAEX_OBS'] = 0.        
#        hdu[0].header['HAEY_OBS'] = -205.64 * iau2015.R_sun.to_value(u.meter)       
#        hdu[0].header['HAEZ_OBS'] = 0.        
        
        
        title = "{0}, dist: {1:.1f} Rsun".format(hdu[0].header['INSTRUME'], 
                                             hdu[0].header['DSUN_OBS'] / iau2015.R_sun.to_value(u.meter))        
        
        mu.dispFits(hdu[0].data, hdu[0].header, log=True, 
                    minmax=minmaxfits, cmap='flag',
                    title=title,
                    windowXYPos=(50, 50))
        
        w = wcs.WCS(hdu[0].header)
        
        
        # -- orientation of the F corona in the HAE
        # Values are from Guillermo's paper ApJ 2017
        #-- The inclination beyond 12 de elongation (from ~1au) of the F-corona 
        # plane of symmetry in the HI-1 FOV *** according to our paper *** is 
        # about 3.8 - 4 deg, rather close to that of Venus' inclination. 
        #
        #-- The ascending node for the same range of elongations is about 80 deg. 
        #
        #-- At inner elongation these values start changing rather abruptly.
        
        
#        nerotmat = np.dot(mu.rotmat(np.deg2rad(80 + 180), 'x'), 
#                          mu.rotmat(np.deg2rad(4), 'y'))
#        
#        
#        q = quaternion.from_rotation_matrix(nerotmat)
#        print(q)
#        
#        rx = np.arctan2(-nerotmat[1, 2], nerotmat[2, 2])
#        ry = np.arcsin(nerotmat[0, 2])
#        rz = np.arctan2(-nerotmat[0, 1], nerotmat[0, 0])
        
#        neang = np.array([rx, ry, rz])        


#        neang = np.deg2rad([80, 0, 3.8])  # Parameters work on COR2 and HI
        neang = np.deg2rad([80, 0, 6])
        print("neang : \n", np.rad2deg(neang))
        
        sim1 = scrt.scraytrace(losnbp=200, modelid=modelid, neang=neang, physics=5)
#        sim1 = scraytrace(losnbp=200, modelid=79, neang=neang,
#                          modparam=[1, 0.2, 50.], physics=5)
                          
#        imsize = [hdu[0].header['NAXIS1'], hdu[0].header['NAXIS2']]
        imsize = [256, 256]
        sim1.setParamFromWCS(hdu[0].header,
                             imsize,
                             frame='HAE')
        sim1.raytrace()
        
        sim1.printRTParam()
        sim1.createWCSHead()
#        print(sim1.wcsrt)
    
#        for k in sim1.fitsrt: print(k, sim1.fitsrt[k])
        
#        print(w)
        
#        print(wcs.WCS(sim1.fitsrt))    
        
        minmaxmod[79] = (1e-13, 1e-8)
#        mu.dispFits(sim1.btot, sim1.fitsrt, log=True, minmax=(5e-15, 1e-10),
        mu.dispFits(sim1.btot, sim1.fitsrt, log=True, 
                    minmax=minmaxmod[modelid],
                    units='B0', cmap='flag',
                    windowXYPos=(1100, 50))

#        mu.dispim(sim1.btot, log=True, minmax=minmaxmod[79],
#                    units='B0', cmap='gray',
#                    windowXYPos=(1100, 50))


        print("neang : \n", np.rad2deg(neang))
        print("Obs lon : {0}deg, lat : {1}deg, dist : {2}Rsun".format(np.rad2deg(sim1.obslonlat[0]), np.rad2deg(sim1.obslonlat[1]), sim1.obslonlat[2]))
    
        for i in ('X', 'Y', 'Z'):
            print("HAE{0}_OBS : {1} Rsun".format(i, hdu[0].header['HAE{0}_OBS'.format(i)] / iau2015.R_sun.to_value(u.meter)))    
    
            
    
        lon_HAE = np.arctan2(hdu[0].header['HAEY_OBS'], hdu[0].header['HAEX_OBS'])
        lat_HAE = np.arcsin(hdu[0].header['HAEZ_OBS'] / hdu[0].header['DSUN_OBS'])
    
        print('lon_HAE : {0} deg'.format(np.rad2deg(lon_HAE)))    
        print('lat_HAE : {0} deg'.format(np.rad2deg(lat_HAE)))    
        print('Dsun_OBS : {0} Rsun'.format(hdu[0].header['DSUN_OBS'] / iau2015.R_sun.to_value(u.meter)))    
    
    
    # ------------------------------------------------------------------------
    elif runname == 'secchi':
    
        minmaxmod = {}
        # Read and display WISPR fits image
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20181101/psp_L2_wispr_20181101T231550_1221.fits')
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190326/psp_L2_wispr_20190326T155032_1221.fits')
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190323/psp_L2_wispr_20190323T155032_1221.fits')
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190405/psp_L2_wispr_20190405T150730_1221.fits')
#        minmaxfits = (1e-12, 5e-11)
#        hdu = fits.open('/net/ares/pub/wispr/fm/ql/Images/fits/L2/20190405/psp_L2_wispr_20190405T152258_2222.fits')
#        minmaxfits = (5e-13, 2e-11)
        minmaxmod[78] = (5e-15, 1e-10)
        
        # COR2 A
#        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/cor2/20120802/20120802_125400_d4c2A.fts')
#        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/cor2/20120530/20120530_235400_d4c2A.fts')
#        minmaxfits = (2e3, 1e4)

        # HI1 A
#        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/hi_1/20181101/20181101_124901_s4h1A.fts')
#        minmaxfits = (2e3, 1e6)
        
        # HI2 A
        hdu = fits.open('/net/cronus/opt/secchi/lz/L0/a/img/hi_2/20181101/20181101_120921_s4h2A.fts')
        minmaxfits = (2e3, 1e6)
#        mimmaxmod78 = (4e-17, 4e-15)
        
#        hdu[0].header['HAEX_OBS'] = 0.        
#        hdu[0].header['HAEY_OBS'] = -205.64 * iau2015.R_sun.to_value(u.meter)       
#        hdu[0].header['HAEZ_OBS'] = 0.        
        
        
        title = "{0}, dist: {1:.1f} Rsun".format(hdu[0].header['INSTRUME'], 
                                             hdu[0].header['DSUN_OBS'] / iau2015.R_sun.to_value(u.meter))        
        
        mu.dispFits(hdu[0].data, hdu[0].header, log=True, 
                    minmax=minmaxfits, cmap='flag',
                    title=title,
                    windowXYPos=(50, 50))
        
        w = wcs.WCS(hdu[0].header)
        
        
        # -- orientation of the F corona in the HAE
        # Values are from Guillermo's paper ApJ 2017
        #-- The inclination beyond 12 de elongation (from ~1au) of the F-corona 
        # plane of symmetry in the HI-1 FOV *** according to our paper *** is 
        # about 3.8 - 4 deg, rather close to that of Venus' inclination. 
        #
        #-- The ascending node for the same range of elongations is about 80 deg. 
        #
        #-- At inner elongation these values start changing rather abruptly.
        
        
#        nerotmat = np.dot(mu.rotmat(np.deg2rad(80 + 180), 'x'), 
#                          mu.rotmat(np.deg2rad(4), 'y'))
#        
#        
#        q = quaternion.from_rotation_matrix(nerotmat)
#        print(q)
#        
#        rx = np.arctan2(-nerotmat[1, 2], nerotmat[2, 2])
#        ry = np.arcsin(nerotmat[0, 2])
#        rz = np.arctan2(-nerotmat[0, 1], nerotmat[0, 0])
        
#        neang = np.array([rx, ry, rz])        
        neang = np.deg2rad([80, 0, 3.8])        
        print("neang : \n", np.rad2deg(neang))
        
#        sim1 = scraytrace(losnbp=100, modelid=78, neang=neang)
        sim1 = scrt.scraytrace(losnbp=200, modelid=79, neang=neang,
                          modparam=[1, 0.5, 170.], physics=5)
                          
#        imsize = [hdu[0].header['NAXIS1'], hdu[0].header['NAXIS2']]
        imsize = [256, 256]
        sim1.setParamFromWCS(hdu[0].header,
                             imsize,
                             frame='HAE')
        sim1.raytrace()
        
        sim1.printRTParam()
        sim1.createWCSHead()
#        print(sim1.wcsrt)
    
#        for k in sim1.fitsrt: print(k, sim1.fitsrt[k])
        
#        print(w)
        
#        print(wcs.WCS(sim1.fitsrt))    
        
        minmaxmod[79] = (1e-13, 1e-8)
#        mu.dispFits(sim1.btot, sim1.fitsrt, log=True, minmax=(5e-15, 1e-10),
        mu.dispFits(sim1.btot, sim1.fitsrt, log=True, minmax=minmaxmod[79],
                    units='B0', cmap='gray',
                    windowXYPos=(1100, 50))

#        mu.dispim(sim1.btot, log=True, minmax=minmaxmod[79],
#                    units='B0', cmap='gray',
#                    windowXYPos=(1100, 50))


        print("neang : \n", np.rad2deg(neang))
        print("Obs lon : {0}deg, lat : {1}deg, dist : {2}Rsun".format(np.rad2deg(sim1.obslonlat[0]), np.rad2deg(sim1.obslonlat[1]), sim1.obslonlat[2]))
    
        for i in ('X', 'Y', 'Z'):
            print("HAE{0}_OBS : {1} Rsun".format(i, hdu[0].header['HAE{0}_OBS'.format(i)] / iau2015.R_sun.to_value(u.meter)))    
    
            
    
        lon_HAE = np.arctan2(hdu[0].header['HAEY_OBS'], hdu[0].header['HAEX_OBS'])
        lat_HAE = np.arcsin(hdu[0].header['HAEZ_OBS'] / hdu[0].header['DSUN_OBS'])
    
        print('lon_HAE : {0} deg'.format(np.rad2deg(lon_HAE)))    
        print('lat_HAE : {0} deg'.format(np.rad2deg(lat_HAE)))    
        print('Dsun_OBS : {0} Rsun'.format(hdu[0].header['DSUN_OBS'] / iau2015.R_sun.to_value(u.meter)))    
    
    # ------------------------------------------------------------------------
    elif runname == 'disk':
        neang = np.deg2rad([80, 30, 0])
        
        sim1 = scrt.scraytrace(losnbp=100, modelid=79, neang=neang,
                          modparam=[1e5, 1., 100.], physics=0,
                          imsize=[256, 256],
                          obspos=[0,0,-215],
                          fovpix=np.deg2rad(30.)/256)
#        sim1.setParamFromWCS(hdu[0].header, [256, 256], frame='HAE')
        sim1.raytrace()
        mu.dispim(sim1.btot, log=True, minmax=(1e-13, 8e-9),
                    units='B0', cmap='gray',
                    windowXYPos=(1100, 50))

#------------------------------------------------------------------------
    elif runname == 'saito':
        save = False
        savePath = './'
        
        neang = np.deg2rad([0, 0, 0])
        
        profNbPix = 512
        fovRadius_Rsun = 30.
        sunAngRadius_deg = 0.27 
        fovpix = np.deg2rad(fovRadius_Rsun * 2. * sunAngRadius_deg) / profNbPix
        losrange = [0, 215. * 2]
        
        df = pd.DataFrame()
        df['modelid'] = [19, 20, 21]
        df['modelName'] = ['K, Equatorial Hole', 'K, Polar Regions (Hole)', 'K, Background (Equator)']
        df['shortName'] = ['K-EquatHole', 'K-PolarRegionsHole', 'K-BackgroundEquator']
                
        sim = {}
        
        for i, row in df.iterrows():
            s = scrt.scraytrace(losnbp=1000, 
                                modelid=row.modelid,
                                losrange=losrange,
                                neang=neang,
                                physics=0,
                                imsize=[1, profNbPix],
                                obspos=[0,0,-215],
                                fovpix=fovpix)
            s.raytrace()
            sim[row.modelid] = s.btot
            
        
        
        r = np.linspace(-fovRadius_Rsun, fovRadius_Rsun, profNbPix)
        m = r > 3
        
        fco = fcor_kl(r[m])

        # --- save data to Excel
        dfOut = pd.DataFrame({'Rsun': pd.Series(r[m]),
                              'Fcorona': pd.Series(fco),
                              df.loc[df.modelid == 19, 'shortName'].to_numpy()[0]: pd.Series(np.reshape(sim[19][m], fco.shape[0])),
                              df.loc[df.modelid == 20, 'shortName'].to_numpy()[0]: pd.Series(np.reshape(sim[20][m], fco.shape[0])),
                              df.loc[df.modelid == 21, 'shortName'].to_numpy()[0]: pd.Series(np.reshape(sim[21][m], fco.shape[0]))})
        if save:
            dfOut.to_excel(savePath + 'K-F-Models_01.xlsx')
            
        figsize = None
        fig, ax = plt.subplots(figsize=figsize)
        for i, row in df.iterrows():
            p = ax.semilogy(r[m], sim[row.modelid][m], label=row.modelName)
        
        p = ax.semilogy(r[m], fco, label='F-Corona')
        
        ax.set_xlabel('Radial Distance [Rsun]')
        ax.set_ylabel('Mean Solar Brightness [Bsun]')
        ax.set_xlim((0, fovRadius_Rsun))
        ax.grid(linestyle='--', which='minor')
        ax.legend()

        plt.show()
        plt.tight_layout()
        
        if save:
            fig.savefig(savePath + 'K-F-Models_01.png')
    
    
    else:
        print('No runname for ', runname)
    