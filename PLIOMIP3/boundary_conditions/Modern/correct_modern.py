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

   


###############################################
# main program

# topography
files = ['LSM','topo','mbiome','soil']
datetoday = date.today()
datestr = str(datetoday)
  

for f in files:
    cube = iris.load_cube('/nfs/hera1/earjcti/PlioMIP3_Boundary_conds_preJuly2023/Modern_std/Modern_std_'+f+'_v1.0.nc')
    # update attributes
    cube.attributes.update({'Netcdf_author' : 'Tindall, J. C.'})  
    cube.attributes.update({'Email' : 'J.C.Tindall@leeds.ac.uk'})  
    cube.attributes.update({'Netcdf_Creation_Date' : datestr}) 
    try:
        del cube.attributes["Code Used"]
    except:
        print('no code used')


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
        cube.attributes.update({'Other_Key':'0=ocean, 1000-1100=Land indicator (1000) + Lake_percentage of grid cell: 1000=0%, 1100=100%'})
        data = cube.data
        newdata = np.where(cube.data > 10000, 0, cube.data)
        cube = cube.copy(data=newdata)

    if f=='soil':
        otherkey = cube.attributes.get("Other_Key")
        cube.attributes.update({"Other_Key": otherkey + ' 0:ocean'})
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

    if f=='soil':
        iris.save(cube,'Modern_std_soil_lake_v1.0.nc')
    else:
        iris.save(cube,'Modern_std_'+f+'_v1.0.nc')
       

