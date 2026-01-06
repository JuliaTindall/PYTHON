#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will correct the EOI400 vegetation for LP vegetation.
It will get rid of the Nans and make sure the LSM matches

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


vegcube = get_veg()
iris.save(vegcube,'LP_mbiome_v1.0.nc',
          fill_value=0.0)


