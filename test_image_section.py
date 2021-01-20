# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 09:57:16 2021

@author: maria
"""

from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import curve_fit

#read data from the file
hdulist = fits.open("Mosaic.fits/mosaic.fits")#plt.imshow(data, cmap='gray')
hdr=hdulist[0].header #meta data 
data=hdulist[0].data #astro image data
hist_data=data.ravel() #made data file a 1d array of pixel values 
hist_data_sorted=sorted(hist_data) #sort data points into ascending order
background_data=[] #empty list - store the relevant background data
hdulist.close()
datatest = data[880:1059,1915:2222] #[y,x]
plt.imshow(np.log(datatest))
testimagex = np.shape(datatest)[0]
testimagey = np.shape(datatest)[1]
masktest = np.zeros((testimagex,testimagey))

datatest1d=datatest.ravel()
datatest1d_sorted = sorted(datatest1d)
#sums values of pixels inside 1st aperture, pixels of galaxy
def galaxyPhotons(x_centre,y_centre, radius):
	photon_vals=[];
	for i in range(x_centre-radius,x_centre+radius):
		for j in range(y_centre-radius,y_centre+radius):
			if i<testimagex and i>=0 and j<testimagey and j>=0: #only consider pixels within the image!
				d=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
				if d<radius: 
					photon_vals.append(datatest[i,j])
				masktest[i,j]=1
	return(np.sum(photon_vals))

#averages values of pixels between 1st and 2nd aperture, local bakcground mean
def localBackground(x_centre,y_centre, initialRadius, secondaryRadius):
	bg_photon_vals=[];
	for i in range(x_centre-secondaryRadius,x_centre+secondaryRadius):#outer x
		for j in range(y_centre-secondaryRadius,y_centre+secondaryRadius):#outer y
			if i<testimagex and i>=0 and j<testimagey and j>=0: #only consider pixels within the image!
				d2=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
				if d2<secondaryRadius and d2>initialRadius: #only consider pixels between the two apertures
					bg_photon_vals.append(datatest[i,j])
	return(np.sum(bg_photon_vals)/len(bg_photon_vals))

def pixelPos(data, value):
	indices = np.where(datatest==value) #gives tuples of coordinates
	listIndices= list(zip(indices[1], indices[0])) # gives pairs of coordinates in form (x,y)
	return listIndices

	
initialRadius = 12
secondaryRadius = 15
galaxyCounts = []
for i in range(len(datatest1d_sorted)): 
	tempPixVal = datatest1d_sorted[-(1+i)] #find highest pixel value
	if tempPixVal > np.mean(datatest)+np.std(datatest):
		tempPixPos = pixelPos(datatest, tempPixVal) #find position of highest pixel value
		for j in range(len(tempPixPos)): #potentially multiple locations with same pixel value
			loc = tempPixPos[j]
			if masktest[loc[1],loc[0]] == 0: #if pixel is available to use 
				galaxyBrightness=galaxyPhotons(loc[1],loc[0],initialRadius)
				galaxyBackground = localBackground(loc[1],loc[0],initialRadius, secondaryRadius)
				galaxyCounts.append(galaxyBrightness-galaxyBackground)
				print(tempPixVal)
plt.figure()
plt.imshow(masktest)			
