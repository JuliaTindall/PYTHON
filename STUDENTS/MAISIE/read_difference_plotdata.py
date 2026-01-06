#!/usr/bin/env python3

# created by Julia Tindall on 19/01/2024

# This program will use iris to read in two netcdf files
# it will find the difference between the files
# it will plot the difference
# 
#  

import matplotlib.pyplot as plt
import numpy as np
import iris
import iris.quickplot as qplt

# read first file  (this is Near Surface Temperature with CO2=280ppm)

pi280_cube = iris.load_cube('E280.NearSurfaceTemperature.allmean.nc')
#print(pi280_cube) - this statement would show you what was in the file


# read second file (this is NearSurfaceTemperature with CO2=400ppm)
pi400_cube = iris.load_cube('E400.NearSurfaceTemperature.allmean.nc')

# find the difference between the two cubes

difference_cube = pi400_cube - pi280_cube
#print(difference_cube)

# plot the difference using iris.quickplot (which works a lot like matplotlib)

qplt.contourf(difference_cube, cmap='RdBu_r',
              levels=np.arange(-3.0, 3.25, 0.25),extend='both')
# all these things will make the plot look nicer
plt.title('difference between E400 and E280')
plt.gca().coastlines()

#plt.show() # sends the plot to the screen
plt.savefig('difference.png')  # puts the plot to a file
