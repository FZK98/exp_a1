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
reduced_data = hdulist[0].data
hist_data=data.ravel() #made data file a 1d array of pixel values 
hist_data_sorted=sorted(hist_data) #sort data points into ascending order
background_data=[] #empty list - store the relevant background data
hdulist.close()

imagey = np.shape(data)[0]
imagex = np.shape(data)[1]

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
    
#print('the pixel val N points away from the mean', xhistogram[muindex+i])
#print('pixel val 1 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ popt[2])).argmin()])
#print('pixel val 2 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ 2*popt[2])).argmin()])
#print('pixel val 2 sig from the mean', xhistogram[np.abs(xhistogram - (popt[1]+ 3*popt[2])).argmin()])
meanBackground = popt[1]
sigmaBackground = popt[2] # this can be altered accordingly - used to determine if galaxy significant

#0 value means pixel can be considered
#1 value means pixel should be ignored
mask = np.zeros((np.shape(data)[0],np.shape(data)[1])) #initiate a mask image (need to swap x and y round)

#masking of foreground pixels including stars, bleeding and noisy edges
for i in range (imagey):
    for j in range(imagex):
        #masking the borders
        if i <= 200 or j <= 200 or j >= 2370 or i>= 4411: 
            mask[i,j] = 1
        # masking pixel values below 2 std devs from the mean 
        if data[i][j] <= popt[1] + 2*popt[2]:
            mask[i][j] = 1
#masking central vertical bloom
mask[:,1425:1447] = 1
#mask horizontal blooming
mask[415:480,1195:1655] = 1
mask[310:380,1020:1705] = 1
mask[200:285,1390:1480] = 1
#masking misc blooming:
mask[3370:7420,770:784] = 1
mask[3200:3279,770:784] = 1
#masking the stars (identified by eye)
mask[1230:1652,3000:3425] = 1
mask[2228:2369,858:978] = 1
mask[2700:2843,916:1031] = 1
mask[3197:3423,712:830] = 1
mask[3713:3813,2097:2180] = 1
mask[3267:3340,2222:2304] = 1
mask[2940:3465,1182:1690] = 1
mask[2268:2346,2091:2171] = 1
mask[1371:1464,2043:2140] = 1
mask[1746:1814,1375:1446] = 1
mask[4015:4061,1429:1490] = 1
mask[2250:2306,2280:2324] = 1
mask[3820:3875,2249:2302] = 1
mask[3819:3879,2242:2308] = 1
mask[2277:2330,425:472] = 1
mask[1465:1516,609:667] = 1
mask[540:605,1734:1813] = 1
mask[4066:4131,525:587] = 1
plt.figure()
plt.imshow(mask)
#plt.figure()
#plt.imshow(np.log10(data))
#print('we have removed this many data points: ', np.count_nonzero(mask))

#sums values of pixels inside 1st aperture, pixels of galaxy
def galaxyPhotons(x_centre,y_centre, radius):
	photon_vals=[];
	for i in range(x_centre-radius,x_centre+radius):
		for j in range(y_centre-radius,y_centre+radius):
			if i<imagex and i>=0 and j<imagey and j>=0: #only consider pixels within the image
				d=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
				if d<radius: 
					photon_vals.append(data[i,j])
				mask[i,j]=1
	return(np.sum(photon_vals)) #returns total photon counts for the galaxy
	
#averages values of pixels between 1st and 2nd aperture, local bakcground mean	
def localBackground(x_centre,y_centre, initialRadius, secondaryRadius):
	bg_photon_vals=[];
	for i in range(x_centre-secondaryRadius,x_centre+secondaryRadius):#outer x
		for j in range(y_centre-secondaryRadius,y_centre+secondaryRadius):#outer y
			if i<imagex and i>=0 and j<imagey and j>=0: #only consider pixels within the image
				d2=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
				if d2<secondaryRadius and d2>initialRadius: #only consider pixels in between the two apertures
					bg_photon_vals.append(data[i,j])
	return(np.sum(bg_photon_vals)/len(bg_photon_vals)) #returns mean local background

def pixelPos(data, value):
	indices = np.where(data==value) #gives tuples of coordinates
	listIndices= list(zip(indices[1], indices[0])) # gives pairs of coordinates in form (x,y)
	return listIndices


initialRadius = 6 #this can be varied as an extension to account for differently sized galaxies
secondaryRadius=10 #this number needs to be better identified, but for now choose arbitrarily
galaxyCounts = []

#reduces the image data set to only include unmasked values
for j in range(imagey):
    for i in range (imagex):
        if mask[j][i] != 0 :
            reduced_data[j][i] = 0


#allData[i] = list(filter(lambda a: a !=0, allData[i]))
reduced_hist_data=reduced_data.ravel() #made data file a 1d reduced array of pixel values
reduced_hist_data = list(filter(lambda a:a !=0, reduced_hist_data)) 
reduced_hist_data_sorted=sorted(reduced_hist_data) #sorts into ascending order
            
    

for i in range(len(reduced_hist_data_sorted)):  
 	print(i/len(reduced_hist_data_sorted)*100)
 	tempPixVal = reduced_hist_data_sorted[-(1+i)] #find next highest pixel value in image
 	if tempPixVal > meanBackground+sigmaBackground: #if the pixel is significantly brighter than background mean, continue
         tempPixPos = pixelPos(data, tempPixVal) #find position of this pixel value
         for j in range(len(tempPixPos)): #potentially multiple locations with same pixel value, so iterate through them
             loc = tempPixPos[j]
             galaxyBrightness=galaxyPhotons(loc[1],loc[0],initialRadius)
             galaxyBackground = localBackground(loc[1],loc[0],initialRadius, secondaryRadius)
             galaxyCounts.append(galaxyBrightness-galaxyBackground)
