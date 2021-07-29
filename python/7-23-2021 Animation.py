#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

This code simulates an expanding CME. 


**Important** This code only runs with the version of mathutil I uploaded to github because I had to
change code in mathutil.py (lines 215-217) for this to work. 


"""




# import unittest
import numpy as np
import copy
from matplotlib import pyplot as plt
from astropy.io import fits
import pathlib
import scraytrace as sc

import mathutil as mu

# ---- Define tests here


def test_CME():
    """Raytrace and display a GCS CM model"""
    
    plt.close('all')
    
    imsize = [256, 256]
    modelid = 54
    # model parameters: see models51to60.cpp, CModel54::dumpDefaultParamForIDL
    # 
    # ("rb","2.55","dist to bottom of structure","Rsun"));    
    # ("alpha","0.52","angle between Oy and leg axis","rad"));    
    # ("h","6.","leg height","Rsun"));    
    # ("kappa","0.4","aspect ratio","")); 
    # ("nemin","1e6","Ne","cm^-3"));  
    # ("thick","0.1","Thickness of the skeleton axis.","Rsun"));
    # ("neaxis","0.","Electron density in the skeleton or axis","cm^-3"));
    # ("stiffness","0.","NOT USED",""));  
    # ("skinsigmain","0.1","inner sigma","")
        
    # ("skinsigmafr","0.1","front sigma",""));
    alpha = 0.52  #0.52
     #this is the time dependence
    kappa = 0.4 #0.4
    leadEdgeHeightInit_Rsun = 2. # starting height of the CME
    leadingEdgeEndHeight_Rsun = 70 # final height of the CME

    sequenceNumberofImages = 20 # set the number of images of the movie sequence

    if sequenceNumberofImages > 1:
        leadEdgeSpeedPerImage = (leadingEdgeEndHeight_Rsun - leadEdgeHeightInit_Rsun) / (sequenceNumberofImages - 1) # "speed" of the CME, from one image to the other
    else:
        leadEdgeSpeedPerImage = leadingEdgeEndHeight_Rsun - leadEdgeHeightInit_Rsun

    
    # modparam = [1.1, alpha, h, kappa, 1e6, 0.2, 0, 0.2, 0.2] #values correspond to description above. First value can be thought of as "time" - distance to sun - length of leg
    physics = 0
    DisttoSun_Rsun = 215. #distance to sun

    fullFOV_deg = 10. # Field of view, edge to edge of the image, defined in deg



    #Define variables we use in sc.scraytrace
    obsang = [0.,0,0]
    obspos = [0., 0., -DisttoSun_Rsun]
    fovpix = np.deg2rad(fullFOV_deg / imsize[0])
    neang = np.deg2rad([170, 5, 48]) # CME 3D direction


    obsang_head = str('[') + str(obsang[0]) + ',' + str(obsang[1]) + ',' + str(obsang[2]) + str(']')
    obspos_head = str('[') + str(obspos[0]) + ',' + str(obspos[1]) + ',' + str(obspos[2]) + str(']')
    neang_head = str('[') + str(neang[0]) + ',' + str(neang[1]) + ',' + str(neang[2]) + str(']')


    movieName = 'Movie05'


    """


    Below this section is where I implement the animation loop.
    It begins by running the animation with the current value for h, and then after it runs, it
    chooses a new value of h for the next iteration.
    I chose an arbitrary 10 times to run it. I also chose to halve h through
    each iteration because, otherwise, the animation is too fast. """


    # -- define display palette
    palette = copy.copy(plt.cm.gray)
    palette.set_under(color='black')
    palette.set_bad(color='black')


    # -- set output path
    foldername = movieName
    rootPath = pathlib.Path('/Users/BrandonBonifacio/SCRaytrace_real/python/Output/test01')
    pathlib.Path(rootPath).joinpath(foldername).mkdir(parents=True, exist_ok=True)

    
    for i in range(sequenceNumberofImages):

        # leadingEdgeHeight = sc.model54_calcLeadingEdgeHeight(h, kappa, alpha) #This calculates the leading edge height
        
        
        # -- compute the position of the CME leading edge for the current image ID
        leadingEdgeHeight = leadEdgeHeightInit_Rsun + leadEdgeSpeedPerImage * i
        
        # -- compute the h parameter from the leadingEdgeHeight and other CME params
        h = sc.model54_calcLegHeight(leadingEdgeHeight, kappa, alpha)
        


        #Below is where we set a new modparam and h for the next iteration
        # Good practice is to have that at the beginning of the loop, not the end
            
        #modparam = [1.1, alpha, h, kappa, 1e6, 0.2, 0, 0.2, 0.2]
        modparam = [1.1, alpha, h, kappa, 1e6, 0.2, 0, 0.2, 0.2]


       # obsang = [0.,0,0]
        #obspos = [0., 0., -DisttoSun_Rsun]
       # fovpix = np.deg2rad(fullFOV_deg / imsize[0])
        

        rt = sc.scraytrace(imsize=imsize,
                        frontinteg=True,
                        losrange=[215-100, 215+50], #slice of space where we perform integration. If CME is outside region, will not integrate
                        losnbp=2000, #How many points we sample for a slice. The more the better, but slower
                        modelid=modelid, 
                        modparam=modparam, 
                        physics=physics, 
                        fovpix=fovpix, #FOV
                        projtypecode=1,
                        obsang=obsang,
                        obspos=obspos, 
                        neang=neang, #How we orient CME in space
                        nbthreads=16, 
                        nbchunks=16,
                        phyparam=[0.58])

        # -- Run the raytrace

        

        
        rt.raytrace()
        print('here')
        

        # -- display image
        # plt.close('all')
        fig, ax = rt.dispim(cmap=palette, minmax=(1e-13, 1e-9))
        
        print("Leading edge height {0} Rsun".format(leadingEdgeHeight))

                #Below is where we set up the file saving part of the code

        # for j in range(len(modparam)):
        #     if j != 2:
        #         foldername += str(modparam[j]) + "_"

        # filename = '#' + str(i) + '_' +str(h) + '.fits'
        filename = '{0}_im{1:03d}.fits'.format(movieName, i)

        hdu = fits.PrimaryHDU()
        hdul = fits.HDUList([hdu])
 
        
 
        # You need to assign the data to the fits file, else you were saving an empty file
        hdul[0].data = rt.btot
        

        # Then you can add keywords to the fits header. This saves the values in the file so that we know how it was generated later on
        hdul[0].header['movname'] = (movieName, 'Movie name')
        hdul[0].header['modelid'] = (modelid, 'Model ID')
        headparam = modparam + [obsang_head, obspos_head, neang_head, fovpix, leadingEdgeHeight]
        for iii in range(len(headparam)):
            hdul[0].header['mdpar{0:02d}'.format(iii)] = (headparam[iii], 'Model param. {0}'.format(iii))
            
        # TODO:
            # add obsang, obspos, neang, fovpix, leadingEdgeHeight to the FITS header
        #!!!! Important to remember for the fits header keywords: they are limited to 8 characters
        
        
        
        fullFilename = pathlib.Path(rootPath).joinpath(foldername).joinpath(filename)
        hdul.writeto(fullFilename, overwrite=True)
        



if __name__ == "__main__" :
    test_CME()


    # check fits file
    print()
    print('Checking fits file...')
    fits_image_filename = pathlib.Path('/Users/BrandonBonifacio/SCRaytrace_real/python/Output/test01/Movie05/Movie05_im002.fits')
    with fits.open(fits_image_filename) as hdul:
        print(repr(hdul[0].header))
        fig, ax = mu.dispim(hdul[0].data, log=True, minmax=(1e-11, 1e-8))
        
        



"""


Next Task:

Save the output into .fits file
Name the images so that we can read them later

Instead of 10, do 20 images - increase integration range

For a small movie, very important to save parameters we used for modelling - autogenerated .txt file in output folder

or .yaml file

Export modparam for each iteration

Do 4 or 5 sequences of CMEs with different parameters:

- neang  -[30,220,0] (towards earth), more or less all towards Earth [0,180,90]
- simulated image for coronagraphs
-

Creating artificial sequences of CME images




"""


        


    
