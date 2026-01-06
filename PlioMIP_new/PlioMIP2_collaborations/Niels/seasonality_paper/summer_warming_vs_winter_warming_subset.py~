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
#
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

       


#################
# MAIN PROGRAM
################


# read in multimodel mean monthly SST data (EOI400-E280)

field = 'NearSurfaceTemperature' # NearSurfaceTemperature, SST, TotalPrecipitation

if field == 'TotalPrecipiation':
    cmapname = 'bwr'
else:
    cmapname = 'rainbow'


anom_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/'+field+'_multimodelmean_month.nc', field + 'plio - pi')

jan_cube = anom_cube.extract(iris.Constraint(time=1))
jul_cube = anom_cube.extract(iris.Constraint(time=7))

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
plt.savefig('summer_minus_winter_'+field+'_NA.png')                             
plt.savefig('summer_minus_winter_'+field+'_NA.eps')

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
plt.savefig('summer_and_winter_orth_'+field+'.png')
plt.savefig('summer_and_winter_orth_'+field+'.eps')

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
plt.savefig('summer_and_winter_NA_'+field+'.png')
plt.savefig('summer_and_winter_NA_'+field+'.eps')


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
plt.savefig('summer_and_winter_'+field+'.png')
plt.savefig('summer_and_winter_'+field+'.eps')
