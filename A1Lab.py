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

imagex = shape(data)[0]
imagey = shape(data)[1]
#stores ascending pixel values upto the value of 6000 
#this excludes intrumental artifacts
#for i in range(len(hist_data_sorted)):
#	if hist_data_sorted[i]<6000:
#		background_data.append(hist_data_sorted[i])

#defines the Gaussian function 
#x = position 
#a = peak height 
#mu = avg 
#sig = std dev
def gaussian_fit(x,a,mu,sig):
    gaussian = a*np.exp(-(x-mu)**2/(2*sig**2))
    return gaussian

#bin frequencies and bins from histogram data (one bin for each pixel value)
yhistogram, xhistogram = np.histogram(hist_data, bins=np.max(hist_data))
#popt has optimized variable values 
#fit the scipy function to the data
popt, pcov = sp.optimize.curve_fit(gaussian_fit, xhistogram[0:-1], yhistogram, [8000000, 3420, 20])
#plot the histogram 
#plt.ylabel("Bin frequency")
#plt.xlabel("Pixel value")
#plt.hist(hist_data,bins=np.max(hist_data))
#i = np.linspace(3000, 4000, 1001)
#plot the gaussian on top
#plt.plot(i, gaussian_fit(i, *popt))
#plt.show()


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
meanBackground = popt[1]
sigmaBackground = popt[2]*1000 # this can be altered accordingly - used to determine if galaxy significant
#0 value means pixel can be considered
#1 value means pixel should be ignored
mask = np.zeros((np.shape(data)[1],np.shape(data)[0])) #initiate a mask image (need to swap x and y round)

def galaxyPhotons(x_centre,y_centre, radius):
	photon_vals=[];
	for i in range(x_centre-radius,x_centre+radius):
		for j in range(y_centre-radius,y_centre+radius):
			if i<imagex and i>=0 and j<imagey and j>=0: #only consider pixels within the image!
				#(i,j)
				d=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
				if d<radius: 
					photon_vals.append(data[j][i])
			mask[i][j]=1
	return(np.sum(photon_vals))

def pixelPos(data, value):
	indices = np.where(data==value) #gives tuples of coordinates
	listIndices= list(zip(indices[1], indices[0])) # gives pairs of coordinates in form (x,y)
	return listIndices


initialRadius = 12 #this can be varied as an extension to account for differently sized galaxies
secondaryRadius=16 #this number needs to be better identified, but for now choose arbitrarily
galaxyBrightnesses = []

for i in range(len(hist_data_sorted)): 
	print(i/len(hist_data_sorted)*100)
	tempPixVal = hist_data_sorted[-(1+i)] #find highest pixel value
	if tempPixVal > meanBackground+sigmaBackground: #if the pixel is significantly brighter than background mean, continue
		tempPixPos = pixelPos(data, tempPixVal) #find position of highest pixel value
		for j in range(len(tempPixPos)): #potentially multiple locations with same pixel value
			loc = tempPixPos[j]
			if mask[loc[1]][loc[0]] == 0: #if pixel is available to use 
				galaxyBrightnesses.append(galaxyPhotons(loc[1],loc[0],initialRadius))
#it would be more efficient to assess mask image earlier, but unsure how to implement this
	else:
		break 
		
