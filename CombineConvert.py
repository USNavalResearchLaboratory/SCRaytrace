
"""
Written by Brandon Bonifacio
Naval Research Laboratory
August 2021

This code assumes that there already exists two fits files modeling a CME with the same parameters at the same time.
One of these fits files is generated with modelid = 54, and the other is generated with modelid = 73. 

Using these two fits files, this code:

- Overlays the images onto eachother to create one fits file
- Resizes the images to H:2048, V:1920 pixels using opencv resize
- Multiplies images by vignetting function
- Convert the images from Bsun to photon per pixel
- Add the photon noise
- Convert the images to DN/s (digital number per second)
- Convert the images to DN, exposure time is 5s.
- Add the column bias
- Save the images


"""

import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import FITS_tools
import numpy as np

plt.style.use(astropy_mpl_style)

#This first section retrieves Model_54 and displays it

model_54 = get_pkg_data_filename('Model_54_Standard_im001.fits')

fits.info(model_54)

model_54_data = fits.getdata(model_54, ext=0)

print(model_54_data.shape)

plt.figure()
plt.imshow(model_54_data, cmap='gray')
plt.colorbar()
plt.show()


#This second section retrieves Model_73 and displays it
model_73 = get_pkg_data_filename('Model_F_Standard_im001.fits')

fits.info(model_73)

model_73_data = fits.getdata(model_73, ext=0)

print(model_73_data.shape)

plt.figure()
plt.imshow(model_73_data, cmap='gray')
plt.colorbar()
plt.show()


#This section overlays Model 54 and Model 73 and shows displays the combined image

combined_data = model_54_data + model_73_data

print(combined_data.shape)

plt.figure()
plt.imshow(combined_data, cmap='gray')
plt.colorbar()
plt.show()


#This section converts the fits file to H:2048, V:1920 and displays it. 

data = np.zeros((1920, 2048), dtype=np.float64)
hdu = fits.PrimaryHDU(data=data)
hdu.writeto('convert.fits', overwrite=True)

convert = get_pkg_data_filename('convert.fits')

fits.info(convert)

convert_data = fits.getdata(convert, ext=0)

print(convert_data.shape)

plt.figure()
plt.imshow(convert_data, cmap='gray')
plt.colorbar()
plt.show()

for i in range(len(convert_data)):
    for j in range(len(convert_data[i])):
        convert_data[i][j] = combined_data[round(i/7.5)-1][round(j/8)-1]

print(convert_data.shape)

plt.figure()
plt.imshow(convert_data, cmap='gray')
plt.colorbar()
plt.show()


#This section multiplies our converted image by the vignetting function. 

vignet = get_pkg_data_filename('Vig-CCOR1_041421.fits')

fits.info(vignet)

vignet_data = fits.getdata(vignet, ext=0)

print(vignet_data.shape)

plt.figure()
plt.imshow(vignet_data, cmap='gray')
plt.colorbar()
plt.show()

for i in range(len(convert_data)):
    for j in range(len(convert_data[i])):
        convert_data[i][j] = convert_data[i][j]*vignet_data[i][j]

print(convert_data.shape)

plt.figure()
plt.imshow(convert_data, cmap='gray')
plt.colorbar()
plt.show()









