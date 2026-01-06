#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created 16 May 2023
#
#@author: earjcti
#"""
#
# Niels noted that for his site the summer warming was more than the
# winter warming.
# He wondered whether this was a global signal or just a local signal.
# This program is to see what the models show
#
# Ammernded September 2023:
# This is based on summer_warming_vs_winter_warming but instead
# of using the multimodel mean we will plot a subset of the models
# these are:
# CCSM4-UoT,  CCSM4_Utr,  COSMOS,  EC-Earth3.3, GISS2.1G,  HadCM3,
# IPSLCM6A,  MIROC and NorESM1-F
#


import cartopy as cart
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import numpy as np
import sys

       

def get_data():
    """
    gets the multimodel mean data from the subset of the models
    """
    filestart = '/nfs/hera1/earjcti/regridded/'
    for i, model in enumerate(MODELNAMES):
        print(i,model)
        cube_e280 = iris.load_cube(filestart + model + '/E280.' 
                                   + field + '.mean_month.nc')
        cube_eoi400 = iris.load_cube(filestart + model + '/EOI400.' 
                                   + field + '.mean_month.nc')
        anom_cube = cube_e280.copy(data=cube_eoi400.data - cube_e280.data)
        iris.util.squeeze(anom_cube)

        if i == 0:
            print(np.shape(anom_cube.data))
            nt, ny,nx = np.shape(anom_cube.data)
            dataarr = np.zeros((len(MODELNAMES), nt, ny, nx))
            
        dataarr[i,:,:,:] = anom_cube.data
            
    dataarr2 = np.where(dataarr > 1.0E10, np.nan, dataarr)
    print(np.shape(dataarr2))
    meandata = np.nanmean(dataarr2,axis=0)
    meancube = anom_cube.copy(data=meandata)
    if field == 'SST':
        meancube.coord('time').points = np.arange(1,13,1)
        meancube.coord('time').attributes = None
        meancube.coord('time').units = None

    return meancube


   
#################
# MAIN PROGRAM
################


# read in multimodel mean monthly SST data (EOI400-E280)

field = 'SST' # NearSurfaceTemperature, SST, TotalPrecipitation

MODELNAMES = ['HadCM3','CCSM4-UoT',  'CCSM4-Utr',  'COSMOS',  'EC-Earth3.3',
              'GISS2.1G',   'IPSLCM6A',  'MIROC4m','NorESM1-F']

if field == 'TotalPrecipiation':
    cmapname = 'bwr'
else:
    cmapname = 'rainbow'

anom_cube = get_data() # get the multimodel mean for MODELNAMES

print(anom_cube)
print(anom_cube.coord('time').points)
print(anom_cube.coord('time'))

if field == 'SST':
    print('here')
    jan_cube = anom_cube.extract(iris.Constraint(time=1))
    print(jan_cube)
    jul_cube = anom_cube.extract(iris.Constraint(time=7))
else:
    jan_cube = anom_cube.extract(iris.Constraint(month=1))
    jul_cube = anom_cube.extract(iris.Constraint(month=7))

summer_min_winter_cube = jul_cube - jan_cube
summer_min_winter_cube.long_name='(mPWP_July - mPWP_Jan) - (PI_July - PI_Jan)'
print(summer_min_winter_cube)

vals=np.arange(-5.0,6.0,1.0)


ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-80,35,20,90],crs=ccrs.PlateCarree())
qplt.contourf(summer_min_winter_cube,cmap='RdBu_r',levels=vals,extend='both')
plt.title(' ')
#plt.gca().coastlines()
if field == 'SST': 
    ax.add_feature(cart.feature.LAND, zorder=1,facecolor='white',edgecolor='black')
else:
    ax.coastlines()
plt.savefig('summer_minus_winter_'+field+'_NA_subset.png')                      
plt.savefig('summer_minus_winter_'+field+'_NA_subset.eps')

sys.exit(0)                             

levels=np.arange(0,11.0,1.0)
if field == 'TotalPrecipitation':
    levels = np.arange(-2.2, 2.4, 0.4)
print(levels)
ax=plt.subplot(1,2,1,projection=ccrs.Orthographic(45,90))
#ax.set_extent([-180,180,45,90])
qplt.contourf(jan_cube,levels=levels,extend='both',cmap='rainbow')
plt.title('January')
plt.gca().coastlines()
ax=plt.subplot(1,2,2,projection=ccrs.Orthographic(45,90))
qplt.contourf(jul_cube,levels=levels,extend='both')
plt.title('July')
plt.gca().coastlines()
plt.savefig('summer_and_winter_orth_'+field+'_subset.png')
plt.savefig('summer_and_winter_orth_'+field+'_subset.eps')

ax=plt.subplot(1,2,1,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
qplt.contourf(jan_cube,levels=levels,extend='both',cmap=cmapname)
plt.title('January (mPWP-PI)')
plt.gca().coastlines()
ax=plt.subplot(1,2,2,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
qplt.contourf(jul_cube,levels=levels,extend='both',cmap=cmapname)
plt.title('July (mPWP - PI)')
plt.gca().coastlines()
plt.savefig('summer_and_winter_NA_'+field+'_subset.png')
plt.savefig('summer_and_winter_NA_'+field+'_subset.eps')


ax=plt.subplot(1,2,1,projection=ccrs.PlateCarree(central_longitude=0.0))
#ax.set_extent([-100,20,20,90],crs=ccrs.PlateCarree())
qplt.contourf(jan_cube,levels=levels,extend='both',cmap=cmapname)
plt.title('January (mPWP-PI)')
plt.gca().coastlines()
ax=plt.subplot(1,2,2,projection=ccrs.PlateCarree(central_longitude=0.0))
#ax.set_extent([-100,20,20,90],crs=ccrs.PlateCarree())
qplt.contourf(jul_cube,levels=levels,extend='both',cmap=cmapname)
plt.title('July (mPWP - PI)')
plt.gca().coastlines()
plt.savefig('summer_and_winter_'+field+'_subset.png')
plt.savefig('summer_and_winter_'+field+'_subset.eps')
