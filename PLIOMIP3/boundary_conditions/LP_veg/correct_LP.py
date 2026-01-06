#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will correct the LP metadata to keep USGS happy

"""

import numpy as np
#import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import netCDF4
from datetime import date
#import cmocean to use topo colormap
import sys

   

##############################################
def get_veg():
    """
    reads in the vegetation and lsm for the LP experiment.
    Makes sure there are no Nans in the vegetation experiment
    also makes sure that the LSM is correct
    """

    lsmfile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc')
    vegfile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_mbiome_v1.0.nc')

    vegcube = iris.load_cube(vegfile,'P3veg_P4Lsm')
    vegdata = vegcube.data
    lsmcube = iris.load_cube(lsmfile)
    lsmdata = lsmcube.data

   
    # look for Nan and correct (copied from LP_Mod_EAIS experiments
    ny,nx=np.shape(vegdata)
    for j in range(0,ny):
        for i in range(0, nx):
            if np.isfinite(vegdata[j,i]):
                pass
            else:
                if vegdata.mask[j,i] == False:
                    print('unknown veg code',i,j,vegcube.coord('longitude').points[i],vegcube.coord('latitude').points[j],vegcube.data[j,i],vegcube.data.mask[j,i],vegdata[j,i],vegdata.mask[j,i])
                    # need to try and correct
                   # print(i,j,np.shape(vegcube.data),np.shape(vegdata.mask))
                   # print('up above',vegcube.data[j-1,i],vegcube.data.mask[j-1,i],vegdata[j-1,i],vegdata.mask[j-1,i])
                   # print('down below',vegcube.data[j+1,i],vegcube.data.mask[j+1,i],vegdata[j+1,i],vegdata.mask[j+1,i])
                   # print('to the left',vegcube.data[j,i-1],vegcube.data.mask[j+1,i],vegdata[j,i-1],vegdata.mask[j,i-1])
                    #if i < 359:
                    #    print('to the right',vegcube.data[j,i+1],vegcube.data.mask[j,i+1],vegdata[j,i+1],vegdata.mask[j,i+1])
                    #print(' ')
                    # we have found that one of the problems is fiji
                    # set to tropical forest (biome0)
                    if i ==359 and j == 73:
                        vegdata[j,i]=1
                    # the other problem is near the Bering Strait
                    # set to 7
                    if j==161:
                        vegdata[j,i]=7

     # look to see where the LSM doesn't match and correct.
    # where lsm cube doesn't mach newvegcube
    # remember lsm=1: land 0: ocean

    veg_01mask = np.where(vegdata.mask==True, 0.0, 1.0)
    for j in range(0,ny):
        for i in range(0, nx):
            #if j==151 and i==10:
            #    print(lsmcube.data[j,i], veg_01mask[j,i],vegdata[j,i])
            if ((lsmdata[j,i] != veg_01mask[j,i]) or
                (lsmdata[j,i] == 1.0 and vegdata[j,i] == 0.0)):
            
                print('found one',i,j,vegdata[j,i],lsmdata[j,i],veg_01mask[j,i])
                if lsmcube.data[j,i] == 0.0: # ocean
                    print('setting land to ocean')
                    vegdata[j,i] = 0.0
                    vegdata.mask[j,i] = True
                else:
                    # we are going to set it to the biome at j-1
                    # tests show this is appropriate
                    vegdata[j,i] = vegdata[j-1,i]
                    vegdata.mask[j,i] = vegdata.mask[j-1,i]
                    print('mask diff',j,i,
                          vegcube.coord('longitude').points[i],
                          vegcube.coord('latitude').points[j],
                          lsmcube.data[j,i], 
                          veg_01mask[j,i],vegdata[j,i])
                    print('surrounddata',vegdata[j-1,i], vegdata[j+1,i],
                          vegdata[j,i-1])
                    if i < 359:
                        print('surrounddata2',vegdata[j,i+1])
                    else:
                        print('surrounddata2',vegdata[j,0])

   
    newvegcube = vegcube.copy(data=vegdata)
   
    return newvegcube

###############################################
# main program

# topography
files = ['icemask','lake','LSM','mbiome','soil','topo']
datetoday = date.today()
datestr = str(datetoday)
  
field = {'icemask' : 'icesheets',
         'lake':'lakes',
         'LSM':'land sea mask',
         'mbiome' : 'vegetation',
         'soil': 'soils',
         'topo' : 'topography'}

for f in files:
    cube = iris.load_cube('/nfs/hera1/earjcti/PlioMIP3_Boundary_conds_preJuly2023/LP/LP_'+f+'_v1.0.nc')
    # update attributes
    cube.attributes.update({'Netcdf_author' : 'Tindall, J. C.'})  
    cube.attributes.update({'Email' : 'J.C.Tindall@leeds.ac.uk'})  
    cube.attributes.update({'Netcdf_Creation_Date' : datestr}) 

    print(f)
    if (f == 'soil' or f == 'lake' or f == 'mbiome' or f == 'icemask'):
        infofield =  'This is the PRISM4 Pliocene ' + field.get(f) + ' using the land sea mask as presented in Dowsett et al (2016)'
    else:
        infofield =  'This is the PRISM4 Pliocene ' + field.get(f) + ' as presented in Dowsett et al (2016)'

    cube.attributes.update({'Information' : infofield})
  

    try:
        del cube.attributes["Code Used"]
    except:
        print('no code used')

    try:
        del cube.attributes["Information2"]
    except:
        print('no information2')

    try:
        del cube.attributes["Code_author"]
    except:
        print('no code author')


    if f=='icemask':
        icedata = cube.data
        newicedata = np.where(icedata > 1000., 0.0, icedata)
        cube=cube.copy(data=newicedata)
  
        cube.attributes["Details"] = "Ocean=0, Land=1, Ice=2 - ice is as PRISM4 (Dowsett et al 2016)" 

    if f=='lake':
        cube.attributes.update({'Other_key':'0=ocean, 1000-1100=Land indicator (1000) + Lake_percentage of grid cell: 1000=0%, 1100=100%'})
        data = cube.data
        newdata = np.where(cube.data > 10000, 0, cube.data)
        cube = cube.copy(data=newdata)

    if f=='soil':
        cube.attributes.update({'Other_Key':'0 = Ocean, 28 = Ice'})
        data = cube.data
        newdata = np.where(cube.data > 10000, 0, cube.data)
        cube = cube.copy(data=newdata)

    if f=='mbiome':
        cube.attributes.update({'Veg_Key_0':'0 = Ocean'})
        data = cube.data
        newdata = np.where(cube.data > 1000, 0, cube.data)
        cube = cube.copy(data=newdata)

    if f=='LSM':
        del cube.attributes["invalid_units"]
       
    print(cube)

    # save 

    iris.save(cube,'LP_'+f+'_v1.0.nc')
       

