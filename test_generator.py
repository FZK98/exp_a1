# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:14:24 2021

@author: maria
"""
no_objects = 100

object_loc = []

image_size = 500
import numpy as np
import matplotlib.pyplot as plt

basearray = np.zeros((image_size,image_size)) #test image
masktest = np.zeros((image_size,image_size))

for i in range(image_size):
	for j in range(image_size):
		basearray[i][j] = 100+np.random.randint(0,100) #assign background vals between 100-200
		
for i in range(no_objects):
	object_loc.append(((np.random.randint(0,image_size)),(np.random.randint(0,image_size))))
	i += 1

for i in range(no_objects):
	loctemp=object_loc[i]
	basearray[loctemp[0]][loctemp[1]] = 1000+np.random.randint(0,500) #assign galaxies value 1000-1500
	basearray[loctemp[0]+1][loctemp[1]] = 1000
	basearray[loctemp[0]][loctemp[1]+1]= 1000
	basearray[loctemp[0]+1][loctemp[1]+1]= 1000
	basearray[loctemp[0]-1][loctemp[1]]= 1000
	basearray[loctemp[0]][loctemp[1]-1]= 1000
	basearray[loctemp[0]-1][loctemp[1]-1]= 1000
	basearray[loctemp[0]+1][loctemp[1]-1]= 1000
	basearray[loctemp[0]-1][loctemp[1]+1]= 1000 #extend galaxies to 3x3 with outer vals 1000
	#may get indexing errors as the randomly chosen locations will be at edges - run multiple times until works
plt.figure()
plt.imshow(basearray)
base1d = basearray.ravel() #make test image 1D
base1d_sorted = sorted(base1d) #sort values in test image in ascending order


def galaxyPhotons(x_centre,y_centre, radius):
	photon_vals=[];
	for i in range(x_centre-radius,x_centre+radius):
		for j in range(y_centre-radius,y_centre+radius):
			if i<500 and i>=0 and j<500 and j>=0: #only consider pixels within the image!
				#(i,j)
				d=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
				if d<radius: 
					photon_vals.append(basearray[j][i])
			masktest[i][j]=1
	return(np.sum(photon_vals))

def localBackground(x_centre,y_centre, initialRadius, secondaryRadius):
	bg_photon_vals=[];
	for i in range(x_centre-secondaryRadius,x_centre+secondaryRadius):
		for j in range(y_centre-secondaryRadius,y_centre+secondaryRadius):
			if i<500 and i>=0 and j<500 and j>=0: #only consider pixels within the image!
				d2=np.sqrt((i-x_centre)**2+(j-y_centre)**2)
				if d2<secondaryRadius: 
					bg_photon_vals.append(basearray[j][i])
			masktest[i][j]=1
	return(np.sum(bg_photon_vals)/len(bg_photon_vals))

def pixelPos(data, value):
	indices = np.where(basearray==value) #gives tuples of coordinates
	listIndices= list(zip(indices[1], indices[0])) # gives pairs of coordinates in form (x,y)
	return listIndices

	
initialRadius = 3
secondaryRadius = 6
galaxyBrightnesses = []
for i in range(len(base1d_sorted)): 
	tempPixVal = base1d_sorted[-(1+i)] #find highest pixel value
	#print(tempPixVal)
	if tempPixVal > 900:
		tempPixPos = pixelPos(basearray, tempPixVal) #find position of highest pixel value
		#print(tempPixPos, tempPixVal)
		for j in range(len(tempPixPos)): #potentially multiple locations with same pixel value
			loc = tempPixPos[j]
			if masktest[loc[1]][loc[0]] == 0: #if pixel is available to use 
				galaxyBrightnesses.append(galaxyPhotons(loc[1],loc[0],initialRadius))
				#masktest[loc[1]][loc[0]] = 1
				print(tempPixVal)
plt.figure()
plt.imshow(masktest)			
