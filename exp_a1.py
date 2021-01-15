# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 20:37:39 2021
@author: maria
"""
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit
hdulist = fits.open("Mosaic.fits/mosaic.fits")
image_data = fits.getdata("Mosaic.fits/mosaic.fits", ext=0)
#plt.imshow(image_data) #this displays the image
#plt.imshow(data, cmap='gray')
hdr=hdulist[0].header
data=hdulist[0].data
hist_data=data.ravel()
hist_data_sorted=sorted(hist_data)
background_data=[]
for i in range(len(hist_data_sorted)):
	if hist_data_sorted[i]<6000:
		background_data.append(hist_data_sorted[i])		
def gaussian_fit(x,a,mu,sig):
    gaussian = a*sp.exp(-(x-mu)**2/(2*sig**2))
    return gaussian
yhistogram, xhistogram = np.histogram(hist_data, bins=10000)
popt, pcov = sp.optimize.curve_fit(gaussian_fit, xhistogram[1:10001], yhistogram, [8000000, 3420, 20])
plt.hist(background_data,bins=1000)
i = np.linspace(3000, 4000, 1001)
plt.plot(i, gaussian_fit(i, *popt))

mask = np.zeros((np.shape(data)[0],np.shape(data)[1]))
def starPhotons(x_centre, y_centre, radius):
	photon_vals = []
	for i in range(x_centre-radius, x_centre+radius):
		for j in range(y_centre-radius,y_centre+radius):
			d = np.sqrt((i-x_centre)**2+(j-y_centre)**2)
			if d<radius and mask[i][j]==0: #how to ensure we dont redo ones?
				photon_vals.append(data[i][j])
				mask[i][j] = 1 #this way the value of the mask is changed only if it's used
	photon_count = np.sum(photon_vals)
	photon_no = len(photon_vals)
	return(photon_count, photon_no)


