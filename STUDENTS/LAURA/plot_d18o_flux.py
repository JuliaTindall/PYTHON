#!/usr/bin/env python3

# created by Julia Tindall on 19/01/2024

# This program will plot the ratio of isotope flux correction to 
# p-e flux correction
# 
#  

import matplotlib.pyplot as plt
import numpy as np
import iris
import iris.quickplot as qplt
import sys

# read julia's file and plot d18o ratio of fresthwater flux
cubes = iris.load('/nfs/hera1/earjcti/um/xpsia/restart/xpsiao#da000000010c1+.nc')
for cube in cubes:
    if cube.var_name == 'field672':
        j_pmine_cube = cube
    if cube.long_name == 'Stash code = 217':
        j_pmine_iso_cube = cube
    if cube.var_name == 'outflow':
        j_river_cube = cube
    if cube.var_name == 'outflow_1':
        j_river_iso_cube = cube


d18o_ratio = iris.util.squeeze(j_pmine_iso_cube / j_pmine_cube)
d18o_river_j = iris.util.squeeze(j_river_iso_cube / j_river_cube)


#plt.subplot(2,2,1)
#qplt.contourf(d18o_ratio)
#print('julia mean is',(d18o_ratio.collapsed(['latitude','longitude'],iris.analy#sis.MEAN).data - 2005.2E-6) / 2005.2E-9,' permil')
#plt.title('julia xpsia - d18o flux corr')


# read laura's file and plot d18o ratio of fresthwater flux
# and river outflux
cubes = iris.load('xprato#da00000570671+.nc')
for cube in cubes:
    print(cube.var_name, cube.long_name)
    if cube.var_name == 'field672':
        l_pmine_cube = cube
    if cube.var_name == 'outflow':
        l_river_cube = cube
    if cube.long_name == 'Stash code = 217':
        l_pmine_iso_cube = cube
    if cube.var_name == 'outflow_1':
        l_river_iso_cube = cube


d18o_ratio_l = iris.util.squeeze(l_pmine_iso_cube / l_pmine_cube)

print(l_river_iso_cube,l_river_cube)
d18o_river_l = iris.util.squeeze(l_river_iso_cube / l_river_cube)

#plt.subplot(2,2,2)
#qplt.contourf(d18o_ratio_l)
#print('laura mean is',(d18o_ratio_l.collapsed(['latitude','longitude'],iris.ana#lysis.MEAN).data - 2005.2E-6) / 2005.2E-9,' permil')
#plt.title('laura xptmb -  d18o flux corr')
#plt.close()

plt.subplot(2,1,1)
vals = np.arange(-100,120,20)
qplt.contourf((d18o_river_l - 2005.2E-6)/2005.2E-9, levels=vals, extend='both',cmap='RdBu_r')
print('laura mean river is',(d18o_river_l.collapsed(['latitude','longitude'],iris.analysis.MEAN).data - 2005.2E-6) / 2005.2E-9,' permil')
plt.title('laura xptmb -  d18o river outflow')

plt.subplot(2,1,2)
qplt.contourf((d18o_river_j - 2005.2E-6)/2005.2E-9, levels=vals, extend='both',cmap='RdBu_r')
print('julia mean river is',(d18o_river_j.collapsed(['latitude','longitude'],iris.analysis.MEAN).data - 2005.2E-6) / 2005.2E-9,' permil')
plt.title('julia xpsia -  d18o river outflow')


plt.show()
sys.exit(0)

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
