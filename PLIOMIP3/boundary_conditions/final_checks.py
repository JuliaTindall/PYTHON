#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created July 2023 by Julia

This program will check the final boundary conditions.
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
    cubeLSM = iris.load_cube(EXPDIR + EXPFILE+'_LSM_v1.0.nc')  # 0 sea 1 land
    #cubeLSM= iris.load_cube('/nfs/hera1/earjcti/PlioMIP3_Boundary_conds/LP/LP_LSM_v1.0.nc')

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


def check_metadata():
    """
    prints out all the cubes to check the metadata
    """
    cube = iris.load_cube(EXPDIR + EXPFILE+'_LSM_v1.0.nc')
    print(cube)
    print(cube.attributes)
    input("Press Enter to continue...")

    cube = iris.load_cube(EXPDIR + EXPFILE+'_topo_v1.0.nc')
    print(cube)
    print(cube.attributes)
    input("Press Enter to continue...")

    if EXP == 'Modern_std':
        cube = iris.load_cube(EXPDIR + EXPFILE+'_soil_lake_v1.0.nc')
        print(cube)
        print(cube.attributes)
        input("Press Enter to continue...")

    cube=iris.load_cube(EXPDIR + EXPFILE+'_icemask_v1.0.nc')
    print(cube)
    print(cube.attributes)
    input("Press Enter to continue...")

    cube=iris.load_cube(EXPDIR + EXPFILE+'_mbiome_v1.0.nc')
    print(cube)
    print(cube.attributes)
    input("Press Enter to continue...")

    cube=iris.load_cube(EXPDIR + EXPFILE+'_lake_v1.0.nc')
    print(cube)
    print(cube.attributes)
  
    input("Press Enter to continue...")

    cube = iris.load_cube(EXPDIR + EXPFILE+'_soil_v1.0.nc')
    print(cube)
    print(cube.attributes)
  
    input("Press Enter to continue...")


def check_ice():
    """
    check if we have ice in the ice mask then we also have ice in 
    the vegetation mask and the soil mask
    """     
    cubeice = iris.load_cube(EXPDIR + EXPFILE+'_icemask_v1.0.nc')  # 2ice
    cubeveg = iris.load_cube(EXPDIR + EXPFILE+'_mbiome_v1.0.nc')  # 28 ice
    cubesoil = iris.load_cube(EXPDIR + EXPFILE+'_soil_v1.0.nc')  # 28ice
    
    icemask = np.where(cubeice.data == 2, 1.0, 0.0)
    veg_icemask = np.where(cubeveg.data == 28.0, 1.0, 0.0)
    soil_icemask = np.where(cubesoil.data == 28.0, 1.0, 0.0)

    print('ice veg equal',np.array_equal(icemask, veg_icemask))
    print('ice soil equal',np.array_equal(icemask, soil_icemask))
    
   

#####################################################
EXP = 'EP' # ep LP_pi-EAIS  PI  LP LP_min_LSM_change Modern_std
EXPFILE = EXP
if EXP == 'LP_min_LSM_change':
    EXPFILE = 'LP_MinLSM'
if EXP == 'PI':
    EXPFILE = 'Modern_std'
EXPDIR = '/nfs/hera1/earjcti/PlioMIP3_Boundary_conds/' + EXP + '/'

check_lsm_ind = 'y'
check_metadata_ind = 'y'
check_ice_ind = 'y'

# check all the files have the same land sea mask
if check_lsm_ind == 'y':
    check_lsm(EXPDIR + EXPFILE+'_topo_v1.0.nc','topography', 't')
    check_lsm(EXPDIR + EXPFILE+'_mbiome_v1.0.nc','veg', '0')
    if EXP == 'Modern_std':
        check_lsm(EXPDIR + EXPFILE+'_soil_lake_v1.0.nc','soil', '0')
    else:
        check_lsm(EXPDIR + EXPFILE+'_soil_v1.0.nc','soil', '0')
        check_lsm(EXPDIR + EXPFILE+'_icemask_v1.0.nc','ice', '0') 
        check_lsm(EXPDIR + EXPFILE+'_lake_v1.0.nc','lake', '0')


if check_metadata_ind == 'y':
    check_metadata()

if check_ice_ind == 'y':
    check_ice()
