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

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import numpy as np

       


#################
# MAIN PROGRAM
################


# read in multimodel mean monthly SST data (EOI400-E280)

mPWP_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/SST_multimodelmean_month.nc','SSTmean_plio')
PI_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/SST_multimodelmean_month.nc','SSTmean_pi')


janmpwp_cube = mPWP_cube.extract(iris.Constraint(time=1))
julmpwp_cube = mPWP_cube.extract(iris.Constraint(time=7))

janpi_cube = PI_cube.extract(iris.Constraint(time=1))
julpi_cube = PI_cube.extract(iris.Constraint(time=7))

pi_contrast_cube = julpi_cube - janpi_cube
mpwp_contrast_cube = julmpwp_cube - janmpwp_cube


levels = np.arange(0,11,1)

ax=plt.subplot(1,2,1,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-100,20,20,90],crs=ccrs.PlateCarree())
qplt.contourf(mpwp_contrast_cube,levels=levels,extend='both')
plt.title('mPWP (July - January)')
plt.gca().coastlines()
ax=plt.subplot(1,2,2,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-100,20,20,90],crs=ccrs.PlateCarree())
qplt.contourf(pi_contrast_cube,levels=levels,extend='both')
plt.title('PI (July - January)')
plt.gca().coastlines()
plt.savefig('seasonal_contrast_NA.png')
plt.savefig('seasonal_contrast_NA.eps')

levels = np.arange(4,12.5,0.5)
ax=plt.subplot(1,2,1,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-20,20,40,70],crs=ccrs.PlateCarree())
qplt.contourf(mpwp_contrast_cube,levels=levels,extend='both')
plt.title('mPWP (July - January)')
plt.gca().coastlines()
ax=plt.subplot(1,2,2,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-20,20,40,70],crs=ccrs.PlateCarree())
qplt.contourf(pi_contrast_cube,levels=levels,extend='both')
plt.title('PI (July - January)')
plt.gca().coastlines()
plt.savefig('seasonal_contrast_NA_east.png')
plt.savefig('seasonal_contrast_NA_east.eps')


levels = np.arange(3,8.5,0.5)
ax=plt.subplot(1,2,1,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-100,-60,20,60],crs=ccrs.PlateCarree())
qplt.contourf(mpwp_contrast_cube,levels=levels,extend='both')
plt.title('mPWP (July - January)')
plt.gca().coastlines()
ax=plt.subplot(1,2,2,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-100,-60,20,60],crs=ccrs.PlateCarree())
qplt.contourf(pi_contrast_cube,levels=levels,extend='both')
plt.title('PI (July - January)')
plt.gca().coastlines()
plt.savefig('seasonal_contrast_NA_west.png')
plt.savefig('seasonal_contrast_NA_west.eps')
