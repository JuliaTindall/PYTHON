#!/usr/bin/env python2.7
#NAME
# create climatological SST
#PURPOSE
#  
# we want to regrid the E280/EOI400 SST onto the atmospheric grid for forcing the atmosphere model
#
#
# Julia June 2025



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import cartopy.crs as ccrs
import cf_units
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess




def fixSST(cubeout,cubemask,dirin):
    """
    this will make sure the LSM on the cube matches the mask ancil file
    """

    #print(cubeout)
    #print(cubeout.data)
    #print(cubeout.data.mask)

    #1. replace masked temperatures with 1000.
    cubeout.data = np.ma.where(cubeout.data < -900., 1000., cubeout.data)

    #1. replace the mask on cubeout
    print(np.shape(cubeout.data.mask))
    for mon in range(0,12):
        cubeout.data.mask[mon,:,:,:] = cubemask.data

    # 2. now the unset points that need setting are set to 1000.
    #    for SST we will get these from the surface temperature

    surffile = dirin + 'xqbwc_Monthly_Average_#pd_temp_1_3900_4000.nc'
    surfcube = iris.load_cube(surffile)
    
    cubedata=cubeout.data
    nt,nz,ny,nx = np.shape(cubedata)
    for t in range(0,nt):
        for k in range(0,nz):
            for j in range(0,ny):
                for i in range(0,nx):
                    if cubedata[t,k,j,i] > 500.:
                        valuse = max([surfcube.data[t,k,j,i] - 273.15,-1.8])
                        cubedata[t,k,j,i]=valuse
    
    return cubeout


def fixice(cubeout,cubemask,dirin):
    """
    this will check the ice is set properly where the LSM has changed
    """

    # if the atmospheric temp < -1.8degC we will set to 100% seaice
    
    
    surffile = dirin + 'xqbwc_Monthly_Average_#pd_temp_1_3900_4000.nc'
    surfcube = iris.load_cube(surffile)
    
    cubedata=cubeout.data
    nt,nz,ny,nx = np.shape(cubedata)
    for t in range(0,nt):
        for k in range(0,nz):
            for j in range(0,ny):
                for i in range(0,nx):
                    if cubedata[t,k,j,i] > -9000.:
                        if (surfcube.data[t,k,j,i]-273.15 < -1.8
                            and cubedata[t,k,j,i]==0.0):
                            cubedata[t,k,j,i]=1.0
    
    return cubeout
\
################################
# main program

time = 'E280_EPLSM'   # E280 / E280_LPLSM / E280_EPLSM
                      #  option also to specify E280 temperautres on a
                      # LP or EP LSM
field = 'SST'    # SST, SeaIceConc
rootdir =  '/nfs/hera1/earjcti/'
#rootdir = '/home/earjcti/'

print(time)

if time == 'E280' or time == 'E280_LPLSM' or time == 'E280_EPLSM':
    dirin = rootdir + 'um/xqbwc/database_averages/'
    filemid = 'xqbwc_Monthly_Average_#pf_'+field+'_3900_4000'
    filein = dirin + filemid

if time == 'E280':
    maskfile = rootdir + 'ancil/preind2/qrparm.mask.nc'
if time == 'E280_LPLSM':
    maskfile = rootdir + 'ancil/P4_enh/P4_enh_qrparm.mask.nc'
if time == 'E280_EPLSM':
    maskfile = rootdir + 'ancil/EP/EP_mask.nc'

cubein = iris.load_cube(filein + '.nc')
print(cubein.coord('time'))
print(cubein.metadata)
#cubein.coord('time').units = 'days since 0001-01-01 00:00'
#cubein.coord('time').calendar = '360_day'
cubein.coord('time').units = cf_units.Unit('days since 0001-01-01 00:00',
                                           calendar='360_day')
print(cubein.coord('time'))
#sys.exit(0)

# get LSM
cubemask = iris.load_cube(maskfile,'LAND MASK (LOGICAL: LAND=TRUE)')

print(cubemask)

cubeout = cubein.regrid(cubemask,iris.analysis.Nearest())
cubeout.data = np.ma.where(cubemask.data == 1, -99999,cubeout.data)
cubeout.data.mask = cubemask.data


# if we are changing the LSM to pliocene then we need to fill in the gaps

if time == 'E280_LPLSM' or time == 'E280_EPLSM':
    if field == 'SST':
        cubeout = fixSST(cubeout,cubemask,dirin)
    if field == 'SeaIceConc':
        cubeout = fixice(cubeout,cubemask,dirin)

print(cubeout)

iris.save(cubeout,filemid + '_' + time + '.nc',
          netcdf_format = "NETCDF3_CLASSIC",
          fill_value = -99999)
