# -*- coding: utf-8 -*-
"""
Created on Mon Jan  14 20:37:39 2021

@author: User
"""
#import relevant libraries
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


#stores ascending pixel values upto the value of 6000 
#this excludes intrumental artifacts
# for i in range(len(hist_data_sorted)):
# 	if hist_data_sorted[i]<6000:
# 		background_data.append(hist_data_sorted[i])

#defines the Gaussian function 
#x = position 
#a = peak height 
#mu = avg 
#sig = std dev
def gaussian_fit(x,a,mu,sig):
    gaussian = a*np.exp(-(x-mu)**2/(2*sig**2))
    return gaussian

#bin frequencies and bins from histogram data (10000 artificial bins)
yhistogram, xhistogram = np.histogram(hist_data, bins=10000)
#popt has optimized variable values 
#fit the scipy function to the data
popt, pcov = sp.optimize.curve_fit(gaussian_fit, xhistogram[1:10001], yhistogram, [8000000, 3420, 20])
#plot the histogram 
plt.ylabel("Bin frequency")
plt.xlabel("Pixel value")
plt.hist(background_data,bins=1000)
i = np.linspace(3000, 4000, 1001)
#plot the gaussian on top
plt.plot(i, gaussian_fit(i, *popt))
plt.show()


#N is the number of noisy data points below the estimated mean
N = sum(yhistogram[0: np.abs(xhistogram - popt[1]).argmin()])
muindex = np.abs(xhistogram - popt[1]).argmin()
i = 0
#finds the value of xhistogram that is N datapoints above the mean 
while N > 0 : 
    N -= yhistogram[muindex - 1]
    i +=1
    
print('the pixel val N points away from the mean', xhistogram[muindex+i])
print('pixel val 1 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ popt[2])).argmin()])
print('pixel val 2 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ 2*popt[2])).argmin()])
print('pixel val 2 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ 3*popt[2])).argmin()])

mask = np.zeros((np.shape(data)[0],np.shape(data)[1]))

for i in range (0, len(data[0])):
    for j in range in range(0, len(data[1])):
        #masking the borders
        if i <= 200 or j <= 200 or i >= 2370 or j>= 4411: 
            mask[i][j] = 1
        #masking pixel values below 2 std devs from the mean 
        if mask[i][j] <= popt[1] + 2*popt[2]:
            mask[i][j] = 1
        #mask the central vertical bloom 
        if 1425 <= i <= 1447:
            mask[i][j] = 1
        #unmask the central brightest star
        if 1230 <= i <= 1652 and 3000 <= j <= 3425:
            mask[i][j] = 0 
        #masking 3 horizontal blooms (with rectangles)
        if 1193 <= i <= 1655 and 418 <= j <= 477:
            mask[i][j] = 1
        if 1020 <= i <= 1706 and 310 <= j <= 378:
            mask[i][j] = 1
        if 1392 <= i <= 1477 and 200 <= j <= 283:
            mask[i][j] = 1
        #masking small blooming patch on the central LHS
        if 770 <= i <= 784:
            if 3370<= j <=3420 :
                mask[i][j] = 1
            if 3200<= j <= 3279:
                mask[i][j] = 1
        
        
        
        
        
def starPhotons(x_centre, y_centre, radius):
    photon_vals = []
    for i in range(x_centre-radius, x_centre+radius):
        for j in range(y_centre-radius,y_centre+radius):
            d = np.sqrt((i-x_centre)**2+(j-y_centre)**2)
            if d<radius and mask[i][j]==0:
                photon_vals.append(data[i][j])
                mask[i][j] = 1
    photon_count = np.sum(photon_vals)
    photon_no = len(photon_vals)
    return(photon_count, photon_no)


