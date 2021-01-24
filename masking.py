# -*- coding: utf-8 -*-
"""
Created on Fri Jan 22 16:42:35 2021

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
hdulist.close()

#number of rows
imagey = np.shape(data)[0]
#number of columns
imagex = np.shape(data)[1]

popt = np.loadtxt('image_parameters.txt')
#consider 1s don't consider 0s
mask = np.ones((imagey,imagex)) #initiate a mask image 

#masking of foreground pixels including stars, bleeding and noisy edges
for i in range (imagey):
    for j in range(imagex):
        #masking the borders
        if i <= 200 or j <= 200 or j >= 2370 or i>= 4411: 
            mask[i,j] = 0
        # masking pixel values below 2 std devs from the mean 
        if data[i][j] <= popt[1] + 2*popt[2]:
            mask[i][j] = 0
        #mask oversaturated pixels
        if data[i][j] >= 50000:
            mask[i][j] = 0

#masking central vertical bloom
mask[:,1425:1447] = 0
#mask horizontal blooming
mask[415:480,1195:1655] = 0
mask[310:380,1020:1705] = 0
mask[200:285,1390:1480] = 0
#masking misc blooming:
mask[3370:7420,770:784] = 0
mask[3200:3279,770:784] = 0
#masking the stars (identified by eye)
mask[2900:3525,1130:1752] = 0
mask[2228:2369,858:978] = 0
mask[2700:2843,916:1031] = 0
mask[3197:3423,712:830] = 0
mask[3713:3813,2097:2180] = 0
mask[3267:3340,2222:2304] = 0
mask[2268:2346,2091:2171] = 0
mask[1371:1464,2043:2140] = 0
mask[1746:1814,1375:1446] = 0
mask[4015:4061,1429:1490] = 0
mask[2250:2306,2280:2324] = 0
mask[3820:3875,2249:2302] = 0
mask[3819:3879,2242:2308] = 0
mask[2277:2330,425:472] = 0
mask[1465:1516,609:667] = 0
mask[540:605,1734:1813] = 0
mask[4066:4131,525:587] = 0
plt.figure()
plt.imshow(mask)
#plt.figure()
#plt.imshow(np.log10(data))
#print('we have this many data points: ', imagex*imagey - np.count_nonzero(mask))
np.savetxt('mask.txt', mask)
