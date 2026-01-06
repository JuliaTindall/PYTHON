#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 18 14:13:50 2019

#@author: earjcti
#
# This program will extract the fields needed by Zhongshi Zhang for the
# PlioMIP MOC paper.  It will extract the monthly averages from
# steve hunters processed data 
#
#
#

import os
import numpy as np
import scipy as sp
#import cf
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import iris.plot as iplt
import sys

#####################################
def extract_fields(filename,fieldnames,fileoutstart,filefieldnames):

    # load required cubes
    cubes=iris.load_cubes(filename,fieldnames)
    
    for i in range(0,len(fieldnames)):
        fileout=fileoutstart+filefieldnames[i]+'.nc'
        if filefieldnames[i]=='SSS':
            subcube=cubes[i]
            cube_extract= subcube.extract(iris.Constraint(z2=5))
            iris.save(cube_extract,fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
        else:
            iris.save(cubes[i],fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)


def avg_MOC(fileMOC,fileout,fieldname):
    cubes=iris.load(fileMOC,fieldname)
  
    print(np.shape(cubes[0]))
    print(cubes)
    print(cubes[0].data)
    avgcube=cubes[0].data
    for i in range(1,len(cubes)):
        subcube=cubes[i]
        subcube_1d=subcube[0,:,50]
        iplt.plot(subcube_1d,color='r')
        avgcube=avgcube+cubes[i].data
        cubedata=cubes[i].data
       
    avgcube=avgcube/len(cubes)
    # put average cube in subcube area
    subcube.data=avgcube
    subcube_1d=subcube[0,:,50]
    iplt.plot(subcube_1d,color='blue')
    iris.save(subcube,fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
    
    

    
   
   

##########################################################
# main program

#####################################################
# this is extracting fields from a mean file

#fieldnames=['OCN TOP-LEVEL TEMPERATURE K',
#            'SALINITY (OCEAN) (PSU)',
#            'AICE : ICE CONCENTRATION',
#            'SALINITY (OCEAN) (PSU)',]
#filefieldnames=['SST','SSS','iceconc','salinity']
#filename='C:\\Users\\julia\\OneDrive\\DATA\\HadCM3_DATA\\Eoi400_2400-2499_Monthly.nc'
#fileoutstart='C:\\Users\\julia\\OneDrive\\DATA\\HadCM3_DATA\\EOI400_2400_2499_Monthly_'

#filename='/nfs/a103/palaeo_share/PlioMIP2/processed/Eoi400_2400-2499_Monthly.nc'
#fileoutstart='/nfs/hera1/earjcti/um/tenvj/pk2/EOI400_2400_2499_Monthly_'
#extract_fields(filename,fieldnames,fileoutstart,filefieldnames)

#print('here')
#filename='/nfs/a103/palaeo_share/PlioMIP2/processed/Preind_E280_2900-2999_Monthly.nc'
#fileoutstart='/nfs/hera1/earjcti/um/tenvo/pk2/E280_2900_2999_Monthly_'
#extract_fields(filename,fieldnames,fileoutstart,filefieldnames)


####################################################
# this is averaging MOC from the MOC scripts
Basin=['Atlantic','Indian','Pacific','Global']

for i in range(0,len(Basin)):
    fieldname='Meridional Overturning Stream Function ('+Basin[i]+')'

#    fileMOC='/nfs/hera1/earjcti/um/tenvo/pk2/tenvoo@pgt*c1.nc'
#    fileout='/nfs/hera1/earjcti/um/tenvo/pk2/E280_avg_'+Basin[i]+'_MOC.nc'
#    avg_MOC(fileMOC,fileout,fieldname)

#    fileMOC='/nfs/hera1/earjcti/um/tenvj/pk2/tenvjo@pgo*c1.nc'
#    fileout='/nfs/hera1/earjcti/um/tenvj/pk2/EOI400_avg_'+Basin[i]+'_MOC.nc'
#    avg_MOC(fileMOC,fileout,fieldname)

# E400
    fileMOC='/nfs/hera1/earjcti/um/tenvq/pk2/tenvqo@pgt*c1.nc'
    fileout='/nfs/hera1/earjcti/um/tenvq/pk2/E400_avg_'+Basin[i]+'_MOC.nc'
    avg_MOC(fileMOC,fileout,fieldname)

#sys.exit(0)
