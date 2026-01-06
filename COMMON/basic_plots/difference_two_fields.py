#NAME
#    difference_two_fields.py
#PURPOSE 
#
#  This program will plot a difference of two fields to the screen

# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
import iris.quickplot as qplt
from iris.experimental.equalise_cubes import equalise_attributes
import sys
#from netCDF4 import Dataset, MFDataset
from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
 
    
#=================================================================
# MAIN PROGRAM STARTS HERE

field = '_p_1'
#field = '_Temperature'
expt1='xoorb'
expt2='xogzl'

file1 = '/nfs/hera1/earjcti/um/' + expt1 + '/database_averages/'+expt1+'_Annual_Average_a@pd'+field+'.nc'
file2 = '/nfs/hera1/earjcti/um/' + expt2 + '/database_averages/'+ expt2+ '_Annual_Average_a@pd'+field+'.nc'

if field == '_p_1':
    FIELDNAME = 'PSTAR AFTER TIMESTEP'
    levels = np.arange(-3000,3500,250)
if field == '_Temperature':
    FIELDNAME = 'TEMPERATURE AT 1.5M'
    levels = np.arange(-3.,3.2,0.2)


cube1 = iris.util.squeeze(iris.load_cube(file1,FIELDNAME))
cube1.long_name = expt1 + field
cube2 = iris.util.squeeze(iris.load_cube(file2,FIELDNAME))
cube2.long_name = expt2 + field

cubediff = cube2 - cube1
cubelist = [cube1,cube2]
qplt.contourf(cubediff, levels, extend='both',cmap='RdBu_r')
plt.gca().coastlines()
plt.title(FIELDNAME + ': xogzl - xoorb')
if field == '_p_1':
    plt.savefig('pressure_diff.png')
    plt.savefig('pressure_diff.pdf')
    iris.save(cubelist,'pressures.nc')
if field == '_Temperature':
    plt.savefig('temp_diff.png')
    plt.savefig('temp_diff.pdf')
    iris.save(cubelist,'temperatures.nc')
plt.close()

#



