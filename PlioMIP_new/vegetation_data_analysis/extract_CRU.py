#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on September 2020
# 
#  This program will extract the cru temperatures for all of the sites  
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
   
    dfs = pd.read_excel(LAND_DATAFILE)
    found = 'n'
    
    colreq = 9
    
    
    for index, row in dfs.iterrows():
        if row[0] == SITE_REQUIRED:
            found = 'y'
            fielddata = row[colreq]
            lat = row[2]
            lon = row[3]
            pass
            
    
    if found == 'n' and SITE_REQUIRED != 'Meighen Island':
        print('couldnot find site in Ulrichs file')
   #     sys.exit(0)
        
    
   
    return fielddata, lat, lon





  
def main():
    """
    calling structure
    a) put cru data into an iris cube
    b) get lat and long from from Ulrichs spreadsheet
    c) get cru temperature at nearest lat and long
    c) printout
    """
    
    cru_cube = iris.load_cube(FILENAME, 'near-surface temperature')

 
    for site in SITES_REQUIRED:
        #    # get data from Ulrichs spreadsheet
    
        proxy_temperature, proxy_lat, proxy_lon = get_land_obs(site)
    
       
        lat_ix = (np.abs(cru_cube.coord('latitude').points - proxy_lat)).argmin()
        lon_ix = (np.abs(cru_cube.coord('longitude').points - proxy_lon)).argmin()
    
        crutemp = cru_cube.data[lat_ix, lon_ix]
   
        print(site, proxy_lat, proxy_lon, 
              cru_cube.coord('latitude').points[lat_ix],
              cru_cube.coord('longitude').points[lon_ix], crutemp)
 #########################################################
# main program


LINUX_WIN = 'l'
EXPTNAME = 'EOI400'

if LINUX_WIN == 'l':
    FILESTART = '/nfs/hera1/earjcti/' 
    FILENAME = (FILESTART + 'regridded/CRUTEMP/' + 
                'E280.NearSurfaceTemperature.allmean.nc')
    LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')
    PLOTTYPE = '_' + EXPTNAME + '.eps'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
    LAND_DATAFILE = (FILESTART + '/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')
    PLOTTYPE = '_' + EXPTNAME + '.png'
 


SITES_REQUIRED = ['Lake Baikal', 'Lake Elgygytgyn', 'Meighen Island',
                  'Beaver Pond',
                  'Flyes Leaf Bed', 'Lost Chicken Mine',
                  'James Bay Lowland',
                  'Pula Maar',
                  'Alpes-Maritimes',
                  'Tarragona',
                  'Rio Maior',
                  'Yallalie, Perth',
                  'Alpes-Maritimes N', 'Alpes-Maritimes S']

main()
