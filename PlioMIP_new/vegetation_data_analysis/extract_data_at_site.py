#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on October 2021
# 
#  This program will extract the temperature cycle at a given site
#1
#
#
#import os
import numpy as np
import pandas as pd
#import scipy as sp
#import cf
import iris
#import iris.util
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import netCDF4
from mpl_toolkits.basemap import Basemap, shiftgrid
#from netCDF4 import Dataset, MFDataset
#import iris.analysis.cartography
#import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
#import cf_units as unit
#from iris.experimental.equalise_cubes import equalise_attributes
import cartopy.crs as ccrs
#import matplotlib.ticker as mticker
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from mpl_toolkits.basemap import Basemap

import sys

 
def get_land_obs(SITE_REQUIRED):
    """
    reads in the spredsheet from ulrich and returns 
    the field, latitude and longitude from the required site. 
    """

    if SITE_REQUIRED == 'Meighen Island':
        fielddata = 10.0
        lat = 80.0
        lon = -99.0

    if SITE_REQUIRED == 'Beaver Pond':
        fielddata = 1.0
        lat = 79.0
        lon = -82.0
   
    if SITE_REQUIRED == 'Flyes Leaf Bed':
        fielddata = 1.0
        lat = 79.0
        lon = -83.0

    if SITE_REQUIRED == 'Lake Elgygytgyn':
        fielddata = 16.0
        lat = 67.0
        lon = -172.0
    
    if SITE_REQUIRED == 'Alpes-Maritimes N':
        lat = 44.0
        lon = -7.19
        fielddata = 1.0
   
    
    if SITE_REQUIRED == 'Alpes-Maritimes S':
        lat = 43.5
        lon = -7.0
        fielddata = 1.0
   
   
    return  lat, lon





  
def main():
    """
    get temperatures from mmm
    get temperatures from individual models
    
    write out annual temperatures
    write out monthly temperatures
    write out seasonal temperatures
    """
    
    proxy_lat, proxy_lon = get_land_obs(SITE_REQUIRED)
    if proxy_lon < 0: proxy_lon = proxy_lon + 360.
    
    field = 'NearSurfaceTemperaturemean_plio'
    mmm_cube = iris.load_cube(NSAT_MMM_FILE,field)
  
    lat_ix = (np.abs(mmm_cube.coord('latitude').points - proxy_lat)).argmin()
    lon_ix = (np.abs(mmm_cube.coord('longitude').points - proxy_lon)).argmin()
    
    mmmtemp = mmm_cube.data[:, lat_ix, lon_ix]
    mmmtemp_avg = np.mean(mmmtemp)
    mmmtemp_djf = (mmmtemp[0] + mmmtemp[1] + mmmtemp[11]) / 3.0
    mmmtemp_mam = (mmmtemp[2] + mmmtemp[3] + mmmtemp[4]) / 3.0
    mmmtemp_jja = (mmmtemp[5] + mmmtemp[6] + mmmtemp[7]) / 3.0
    mmmtemp_son = (mmmtemp[8] + mmmtemp[9] + mmmtemp[10]) / 3.0
   
    print(proxy_lat, proxy_lon, 
          mmm_cube.coord('latitude').points[lat_ix],
          mmm_cube.coord('longitude').points[lon_ix])

    print('mean annual temp = ',mmmtemp_avg)
    print('mean annual djf = ',mmmtemp_djf)
    print('mean annual mam = ',mmmtemp_mam)
    print('mean annual jja = ',mmmtemp_jja)
    print('mean annual son = ',mmmtemp_son)

    allmodel_summer = []
    allmodel_mean = []
    for model in MODELNAMES:
        filename = (FILESTART + 'regridded100/' + model +  
                    '/EOI400.NearSurfaceTemperature.mean_month.nc')
        model_cube = iris.load_cube(filename)
  
        lat_ix = (np.abs(model_cube.coord('latitude').points - proxy_lat)).argmin()
        lon_ix = (np.abs(model_cube.coord('longitude').points - proxy_lon)).argmin()
        modeltemp = model_cube.data[:, lat_ix, lon_ix]
  
    
        modeltemp_jja = (modeltemp[5] + modeltemp[6] + modeltemp[7]) / 3.0
        modeltemp_mean = np.mean(modeltemp)
        allmodel_summer.append(modeltemp_jja)
        allmodel_mean.append(modeltemp_mean)
      #  print('jja temp for ',model,' = ',modeltemp_jja)
        print('ann temp for ',model,' = ',modeltemp_mean)
  
   # print('these models median jja',np.median(np.asarray(allmodel_summer)))
   # print('these models 20th percentile jja',np.percentile(np.asarray(allmodel_summer),20))
   # print('these models 80th percentile jja',np.percentile(np.asarray(allmodel_summer),90))
  
    print('these models median ann',np.median(np.asarray(allmodel_mean)))
    print('these models 20th percentile ann',np.percentile(np.asarray(allmodel_mean),20))
    print('these models 80th percentile ann',np.percentile(np.asarray(allmodel_mean),80))

 #########################################################
# main program


LINUX_WIN = 'l'
EXPTNAME = 'EOI400'

LINUX_WIN = 'l'
FILESTART = '/nfs/hera1/earjcti/'
DATABASE = '/nfs/hera1/pliomip2/data/'

MODELNAMES = [
               'HadGEM3', 'CESM2',
              'IPSLCM6A', 
              'COSMOS', 
              'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
              'MIROC4m', 'IPSLCM5A2', 'HadCM3',
              'GISS2.1G', 
               'CCSM4', 
              'CCSM4-Utr', 'CCSM4-UoT', 
              'NorESM-L',  'NorESM1-F'
            ,  'MRI2.3'
              ]


NSAT_MMM_FILE = (FILESTART + 'regridded100/' + 
                 'NearSurfaceTemperature_multimodelmean_month.nc')


SITE_REQUIRED= 'Beaver Pond'

main()
