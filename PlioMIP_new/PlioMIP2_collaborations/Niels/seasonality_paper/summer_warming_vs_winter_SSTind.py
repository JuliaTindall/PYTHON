#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created 16 May 2023
#
#@author: earjcti
#"""
#
# Niels noted that for his site the summer warming was more than the
# winter warming.  We have been looking to see if this is related to clouds.
#
# I have plotted a multimodel mean from clouds but we don't have all the models
#  This will do a partial MMM plot for SST

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import iris
import iris.quickplot as qplt
from iris.cube import CubeList
import iris.plot as iplt
import numpy as np
import sys

def get_MMM_cube(fileend):
    """
    gets the multimodel mean from fileend
    """

    filestart = '/nfs/hera1/earjcti/regridded/'
    allcubes = CubeList([])

    for i, model in enumerate(MMM_modelnames):
        filename = filestart + model + fileend
        cube = iris.load_cube(filename)
            
        if i == 0:
            nt, ny, nx = np.shape(cube.data)
            alldata = np.zeros((len(MMM_modelnames), nt, ny ,nx))
        alldata[i,:,:,:] = cube.data

    alldatamean = np.mean(alldata, axis=0)
    MMM_cube = cube.copy(data=alldatamean)
    
    print(MMM_cube.data)

    return MMM_cube

#################
# MAIN PROGRAM
################



field = 'NearSurfaceTemperature'
model='MMMsubset'

MMM_modelnames = ['CESM2','HadCM3','IPSLCM6A','MIROC4m'] # to be used in the MMM

label = {'SST' : 'degC', 'NearSurfaceTemperature' : 'degC',
         'TotalPrecipitation' : 'mm/day'}

if model == 'MMMsubset':
    eoi400_cube = get_MMM_cube('/EOI400.' + field + '.mean_month.nc')
    e280_cube = get_MMM_cube('/E280.' + field + '.mean_month.nc')
else:
    eoi400_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/'+model+'/EOI400.' + field + '.mean_month.nc')
    e280_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/'+model+'/E280.' + field + '.mean_month.nc')


janpi_cube = e280_cube.extract(iris.Constraint(month=1))
julpi_cube = e280_cube.extract(iris.Constraint(month=7))
janplio_cube = eoi400_cube.extract(iris.Constraint(month=1))
julplio_cube = eoi400_cube.extract(iris.Constraint(month=7))

janpi_cube.data = np.where(janpi_cube.data > 1.0e10, np.nan, janpi_cube.data)
julpi_cube.data = np.where(julpi_cube.data > 1.0e10, np.nan, julpi_cube.data)
janplio_cube.data = np.where(janplio_cube.data > 1.0e10, np.nan, janplio_cube.data)
julplio_cube.data = np.where(julplio_cube.data > 1.0e10, np.nan, julplio_cube.data)


levels=np.arange(-10., 40., 5.)

if field == 'TotalPrecipitation':
    levels=levels/10.
ax=plt.subplot(3,2,1,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
axplot=qplt.contourf(janplio_cube,levels=levels,extend='both')
axplot.cmap.set_under('white')
plt.title('January (MP):' + model)
plt.gca().coastlines()

ax=plt.subplot(3,2,2,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
axplot=qplt.contourf(julplio_cube,levels=levels,extend='both')
axplot.cmap.set_under('white')
plt.title('July (MP):' + model)
plt.gca().coastlines()

ax=plt.subplot(3,2,3,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
axplot=qplt.contourf(janpi_cube,levels=levels,extend='both')
axplot.cmap.set_under('white')
plt.title('January (PI):' + model)
plt.gca().coastlines()

ax=plt.subplot(3,2,4,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
axplot=qplt.contourf(julpi_cube,levels=levels,extend='both')
plt.title('July (pi):' + model)
plt.gca().coastlines()

levels=np.arange(-7, 9, 2)
if field == 'TotalPrecipitation':
    levels=levels/10.

ax=plt.subplot(3,2,5,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
axplot=qplt.contourf(janplio_cube - janpi_cube,levels=levels,extend='both',cmap='RdBu_r')
plt.title('January (mpwp-pi):' + model)
plt.gca().coastlines()

ax=plt.subplot(3,2,6,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-90,20,20,70],crs=ccrs.PlateCarree())
axplot=qplt.contourf(julplio_cube - julpi_cube,levels=levels,extend='both',cmap='RdBu_r')
plt.title('July (mpwp-pi):' + model)
plt.gca().coastlines()

plt.tight_layout()
plt.savefig('summer_and_winter_' + field + '_' + model + '.png')
plt.savefig('summer_and_winter_' + field + '_' + model + '.eps')

plt.close()
### plot jul-jan plio - jul-jan pi
anom_anom_cube = (julplio_cube - janplio_cube) - (julpi_cube - janpi_cube)
#levels=np.arange(-22, 26, 4)
cs = iplt.contourf(anom_anom_cube,levels=levels,cmap='bwr',extend='both')
cbar = plt.colorbar(cs,orientation='horizontal',label=label.get(field),ticks=levels)


plt.gca().coastlines()
plt.title(field + ':' +model+' (mPWP_jul - MPWP_jan) - (PI_jul - PI_jan)')
plt.savefig('anom_' + field + '_'+ model + '.png')
plt.savefig('anom_' + field + '_'+ model + '.eps') 
plt.close()

### plot jul_plio-jul_pi  and jan_plio - jan_pi on one plot
julcube = (julplio_cube - julpi_cube)
jancube = (janplio_cube - janpi_cube)
levels=np.arange(1, 11, 1)
if field == 'TotalPrecipitation':
    levels=levels/10.

plt.subplot(1,2,2)
cs = iplt.contourf(julcube,levels=levels,cmap='gist_heat_r',extend='both')
cs.cmap.set_under('white')
cbar = plt.colorbar(cs,orientation='horizontal',label=label.get(field),ticks=levels)
plt.title('jul ' +field+' anom mpwp -pi')
plt.gca().coastlines()
plt.subplot(1,2,1)
cs = iplt.contourf(jancube,levels=levels,cmap='gist_heat_r',extend='both')
cs.cmap.set_under('white')
cbar = plt.colorbar(cs,orientation='horizontal',label=label.get(field),ticks=levels)
plt.title('jan '+field+' anom mpwp -pi')
plt.gca().coastlines()

plt.savefig('jul_jan_anom_' + field + '_'+ model + '.eps') 
plt.savefig('jul_jan_anom_' + field + '_' +  model + '.png')

plt.close()


# plot over north atlantic only

levels=np.arange(-7, 9, 2)
if field == 'TotalPrecipitation':
    levels=levels/10.

lat_constraint = iris.Constraint(latitude = lambda cell:  20 < cell < 70.0)  

lon_slice  =  anom_anom_cube.intersection(longitude=(-90., 20.))       
anom_anom_region_cube = lon_slice.extract(lat_constraint) 

cs = iplt.contourf(anom_anom_region_cube,levels=levels,cmap='bwr',extend='both')
cbar = plt.colorbar(cs,orientation='horizontal',label=label.get(field),ticks=levels)


plt.gca().coastlines()
plt.savefig('anom_anom_'+ field + '_' + model + '_NA.png')
plt.savefig('anom_anom_'+ field + '_' + model + '_NA.eps') 
