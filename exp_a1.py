# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 20:37:39 2021

@author: maria
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
hdulist = fits.open("Mosaic.fits/mosaic.fits")
image_data = fits.getdata("Mosaic.fits/mosaic.fits", ext=0)
#plt.imshow(image_data) #this displays the image
data=hdulist[0].data
hist_data1=data.ravel()
hdr=hdulist[0].header
plt.figure()
plt.hist(hist_data1, bins=data.max())

plt.figure()
plt.imshow(data, cmap='gray')

