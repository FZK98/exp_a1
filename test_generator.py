# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 16:14:24 2021

@author: maria
"""
no_objects = 50
object_loc = []
image_size = 500
import numpy as np
import matplotlib.pyplot as plt
basearray = np.zeros((image_size,image_size))

for i in range(image_size):
	for j in range(image_size):
		basearray[i][j] = np.random.randint(0,100)
		
for i in range(no_objects):
	object_loc.append(((np.random.randint(0,image_size)),(np.random.randint(0,image_size))))
	i += 1

for i in range(no_objects):
	loctemp=object_loc[i]
	basearray[loctemp[0]][loctemp[1]] = 1000
	
def sersic(r,b):
	return 1000*np.exp(-b*(r**0.25))

plt.imshow(basearray)

def specific_Pixel_Pos(data, value):
	return np.where(data==value) #gives tuple of xy coordinates
def max_Pixel_Val(data):
	return np.max(data)
def max_Pixel_Pos(data):
	return np.where(data==max_Pixel_Val(data)) #gives tuple of xy coordinates 
lowerp = np.percentile(data, 1)
upperp = np.percentile(data,99)

#testing the number of pixels with a certain value to troubleshoot the histogram
#shape1 = shape(data)		
#counter=0		 
#for i in range(shape1[0]):
#	for j in range(shape1[1]):
#		if data[i][j] == 3421:
#			counter+=1			
#print(counter)