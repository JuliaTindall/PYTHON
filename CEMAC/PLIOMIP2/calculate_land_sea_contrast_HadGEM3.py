#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Created on Fri Sep 18 10:42:28 2020

IPCC were not happy with the way that we had done the land sea contrast in the paper.
This program will calculate it based on the individual models land sea mask.
This superceeds extract_ipcc_data.py

@author: julia
"""


import numpy as np
import pandas as pd
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import iris.analysis.cartography
import iris.coord_categorisation
import cf_units as unit
from iris.experimental.equalise_cubes import equalise_attributes
import netCDF4
import sys
#import os


###########################
def get_land_sea_mask():
    """
    the land mask is where the land_frac = 100% in both pliocene & pi
    the sea mask is where the sea_frac = 100% in both pliocene & pi
    returns land_mask and sea_mask as a cube
    """
    f = Dataset(LSM, "r")
    lsm_data = f.variables['lsm'][:]
    lats = f.variables['latitude'][:] 
    lons = f.variables['longitude'][:]
    f.close()
  
    lsm_data = np.squeeze(lsm_data)
    land_mask = lsm_data
    sea_mask = (lsm_data -1.0) * 1.0

    return land_mask, sea_mask, lats, lons

    
def get_nsat_data():
    """
    get the average temperature from the pliocene and the preindustrial
    """

    f = Dataset(FILENAME_PLIO, "r")
    plio_data_all = f.variables[FIELDNAME][:]
    lats = f.variables['latitude'][:] 
    lons = f.variables['longitude'][:]
    f.close()

    plio_data = (np.mean(plio_data_all,axis=0)-273.15)

    f = Dataset(FILENAME_PI, "r")
    pi_data_all = f.variables[FIELDNAME][:]
    lats2 = f.variables['latitude'][:] 
    lons2 = f.variables['longitude'][:]
    f.close()

    pi_data = (np.mean(pi_data_all,axis=0)-273.15)


    for i, lat in enumerate(lats):
        if lat != lats2[i]:
            print('lats dont match', lat, i, lats2[i])
            sys.exit(0)

    for i, lon in enumerate(lons):
        if lon != lons2[i]:
           print('lons dont match', lon, i, lons2[i])
           sys.exit(0)

  
    return plio_data, pi_data, lats, lons


def get_global_avg(land_mask, sea_mask, dataarr, lats, lons):
    """
    gets global average temperature, and also global avg for
    the land and the ocean
    """
  
    grid_areas_lat = np.zeros(len(lats))
    for j, lat in enumerate(lats):
        grid_areas_lat[j] = np.cos(2. * np.pi * lat / 360.)

    print(np.shape(land_mask))
    global_mean = 0.
    global_mean_weights = 0.
    global_mean_land = 0.
    global_mean_land_weights = 0.
    global_mean_sea = 0.
    global_mean_sea_weights = 0.

    for j in range(0, len(lats)):
        for i in range(0, len(lons)):
            global_mean = global_mean + (dataarr[j, i] * grid_areas_lat[j])
            global_mean_weights = global_mean_weights + grid_areas_lat[j]
            if land_mask[j, i] == 1.0:
                global_mean_land = (global_mean_land + 
                                   (dataarr[j, i] * grid_areas_lat[j]))
                global_mean_land_weights = (global_mean_land_weights + 
                                            grid_areas_lat[j])
            else:
                global_mean_sea = (global_mean_sea + 
                                   (dataarr[j, i] * grid_areas_lat[j]))
                global_mean_sea_weights = (global_mean_sea_weights + 
                                            grid_areas_lat[j])

    global_mean = global_mean / global_mean_weights
    global_mean_land = global_mean_land / global_mean_land_weights
    global_mean_sea = global_mean_sea / global_mean_sea_weights

   
    return global_mean, global_mean_land, global_mean_sea, grid_areas_lat


def get_regional_landsea(rmax, rmin, land_mask, dataarr, lats, lons,
                         grid_areas):

    """
    gets the mean temperature for latitude bands for average and for
    land and sea  
    """

    grid_areas_use = grid_areas * 1.0
    for j, lat in enumerate(lats):
        if rmin > lat or rmax < lat:
            grid_areas_use[j] = 0.0
    #grid_areas_land = grid_areas_use * land_cube.data
    #grid_areas_sea = grid_areas_use * sea_cube.data
    
    reg_mean = 0.
    reg_mean_weights = 0.
    reg_mean_land = 0.
    reg_mean_land_weights = 0.
    reg_mean_sea = 0.
    reg_mean_sea_weights = 0.

    for j in range(0, len(lats)):
        if grid_areas_use[j] != 0.0:
            for i in range(0, len(lons)):
                reg_mean = reg_mean + (dataarr[j, i] * grid_areas_use[j])
                reg_mean_weights = reg_mean_weights + grid_areas_use[j]
                if land_mask[j, i] == 1.0:
                    reg_mean_land = (reg_mean_land + 
                                     (dataarr[j, i] * grid_areas_use[j]))
                    reg_mean_land_weights = (reg_mean_land_weights + 
                                             grid_areas_use[j])
                else:
                    reg_mean_sea = (reg_mean_sea + 
                                    (dataarr[j, i] * grid_areas_use[j]))
                    reg_mean_sea_weights = (reg_mean_sea_weights + 
                                            grid_areas_use[j])
    print(rmax, rmin)
    print(rmax, rmin, reg_mean / reg_mean_weights)
    print(reg_mean, reg_mean_land, reg_mean_sea)
    print(reg_mean_weights, reg_mean_land_weights, reg_mean_sea_weights)
  
    reg_mean = reg_mean / reg_mean_weights
    reg_mean_land = reg_mean_land / reg_mean_land_weights
    reg_mean_sea = reg_mean_sea / reg_mean_sea_weights

    return [reg_mean, reg_mean_land, reg_mean_sea]
           
                                                       
def write_to_spreadsheet(avg_T_anom, avg_T_landanom, avg_T_seaanom,
                         land_sea_anom, regionmax, regionmin):
    """
    write the information to a pandas dataframe
    """
    print(avg_T_anom, MODELNAME)

    data = [['Global', avg_T_anom], ['Global (over land)', avg_T_landanom],
            ['Global (over sea)', avg_T_seaanom]]

    for i, rmax in enumerate(regionmax):
        if rmax > 0:
            latrange = np.str(np.around(rmax)) + 'N'
        else:
            latrange = np.str(np.around(np.abs(rmax))) + 'S'
        if regionmin[i] > 0.:
            latrange = latrange + np.str(np.around(regionmin[i])) + 'N'
        else:
            latrange = latrange + np.str(np.around(np.abs(regionmin[i]))) + 'S'
        
        data.append(['glob_' + latrange, land_sea_anom[0,i]])
        data.append(['land_' + latrange, land_sea_anom[1,i]])
        data.append(['sea_' + latrange, land_sea_anom[2,i]])
        

    df = pd.DataFrame(data, columns = ['Simulated temperature', MODELNAME])
 

    # save dataframe as a excel file
    filename = ('/nfs/hera1/earjcti/PLIOMIP2/IPCC/' + MODELNAME + 
          'mPWP_CMIP6_land_sea.csv')
    #df.to_excel(filename)
    df.to_csv(filename)


##############################################
def get_land_sea_contrast():
    """
    get the land sea contrast for this model
    """

    print('moodelname is', MODELNAME)
    print('filename is', FILENAME_PLIO)
    print('lsm is', LSM)

    # get land and sea mask
    land_mask, sea_mask, lsm_lats, lsm_lons = get_land_sea_mask()

    # get temporally averaged nsat data
    print('getting temporally averaged nsat data')
    plio_data, pi_data, lats, lons = get_nsat_data()

    # check grid
    for i, lat in enumerate(lsm_lats):
        if lat != lats[i]:
           print('lsm lat doesnt match', i, lat, lats[i])
           sys.exit(0)
    for i, lon in enumerate(lsm_lons):
        if lon !=lons[i]:
           print('lsm lon doesnt match', i, lon, lons[i])
    

    # get land and sea temperatures
    print('getting land sea temperatures')
    (avg_T_plio, avg_land_T_plio, 
     avg_sea_T_plio, grid_areas) = get_global_avg(land_mask, sea_mask, 
                                                 plio_data, lats, lons)
                                    
    (avg_T_pi, avg_land_T_pi,
     avg_sea_T_pi, grid_areas) = get_global_avg(land_mask, sea_mask,
                                                pi_data, lats, lons)

    avg_T_anom = avg_T_plio - avg_T_pi
    avg_T_landanom = avg_land_T_plio - avg_land_T_pi
    avg_T_seaanom = avg_sea_T_plio - avg_sea_T_pi



    regionmax = [90.0, 60.0, 30.0, 0.0, -30.0, -60.0]
    regionmin = [60.0, 30.0, 0.0, -30.0, -60.0, -90.0]
    land_sea_region_pi = np.zeros((3, len(regionmax)))
    land_sea_region_plio = np.zeros((3, len(regionmax)))


    for i, rmax in enumerate(regionmax):
        land_sea_region_pi[:, i] = get_regional_landsea(rmax, regionmin[i],
                                                        land_mask, pi_data,
                                                        lats, lons, 
                                                        grid_areas)
        land_sea_region_plio[:, i] = get_regional_landsea(rmax, regionmin[i],
                                                        land_mask, plio_data,
                                                        lats, lons,
                                                        grid_areas)

       
    
 # write to a spreadsheet
    land_sea_anom = land_sea_region_plio - land_sea_region_pi

    write_to_spreadsheet(avg_T_anom, avg_T_landanom, avg_T_seaanom,
                         land_sea_anom, regionmax, regionmin)
    


    

   


##########################################################
# main program

FILENAME  =  ' '
LINUX_WIN  =  'l'
MODELNAME  = "HadGEM3"
FIELDNAMEIN = ['tas']


START = '/nfs/hera1/pliomip2/data/HadGEM3_new/'
FILENAME_PI = START + 'climatologies/E280/atmos/clims_hadgem3_pi_airtemp_final.nc'
FILENAME_PLIO = START + 'climatologies/Eoi400/atmos/clims_hadgem3_pliocene_airtemp_final.nc'
FIELDNAME = 'temp'
LSM = START + 'hadgem3.mask.nc'
FIELDLSM = 'land_binary_mask'#

if LINUX_WIN  == 'l':
    FILESTART = '/nfs/hera1/pliomip2/data/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'

get_land_sea_contrast()

#sys.exit(0)
