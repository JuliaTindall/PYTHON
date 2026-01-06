#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created July 2023 by Julia

This program will Check the Early Pliocene boundary conditions
"""

import numpy as np
import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import netCDF4
from datetime import date
#import cmocean to use topo colormap
import sys

def check_lsm(file_to_check, details, mask_ind):
    """
    check all the files created in this section have the same land sea mask
    details is what will be written out to show what section we are in
    mask_ind = 'm' if sea is masked out
    mask_ind = 0 if sea is set to zero
    mask_ind = 't' if sea is < 1 land > 1
    """     
    print(' ')
    print(details) 
    cubeLSM = iris.load_cube('EP_LSM_v1.0.nc')  # 0 sea 1 land
    cubecheck = iris.load_cube(file_to_check)
    if mask_ind == 'm':
        filelsm = np.where(cubecheck.data.mask == True, 0.0, 1.0)
    elif mask_ind == '0':
        filelsm = np.where(cubecheck.data == 0, 0.0, 1.0)
    elif mask_ind == 't':
        filelsm = np.where(cubecheck.data <= 0.0, 0.0, 1.0)
    else:
        print('wrong mask indicator')
        sys.exit(0)
      
        
    diff = cubeLSM.data - filelsm
    #print(cubecheck.data)
    #plt.subplot(211)
    #qplt.contourf(cubeLSM)
    #plt.subplot(212)
    #qplt.contourf(cubeLSM.copy(data=filelsm))
    #plt.show()
    #sys.exit(0)

    differences = np.sum(np.where(diff != 0.0, 1.0, 0.0))
    print('number of different gridboxes =' + np.str(differences))
    
    if differences > 0.0:
        if details == "veg":
            print('veg')
        else:
            print('julia',details)
        diffcube = cubeLSM.copy(data=diff)
        for j, lat in enumerate(diffcube.coord('latitude').points):
            for i, lon in enumerate(diffcube.coord('longitude').points):
                if diffcube.data[j,i] != 0.0:
                    print(lat, lon, j, i, diffcube.data[j,i], cubeLSM.data[j,i], filelsm[j,i],cubecheck.data[j,i])



# check all the files have the same land sea mask
check_lsm('EP_topo_v1.0.nc','topography', 't')
check_lsm('EP_icemask_v1.0.nc','ice', 'm')
check_lsm('EP_mbiome_v1.0.nc','veg', 'm')
check_lsm('EP_soil_v1.0.nc','soil', 'm')
check_lsm('EP_lake_v1.0.nc','lake', 'm')
