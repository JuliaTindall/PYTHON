#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#Created on Fri Jul  5 15:11:26 2019
#Updated for JH on 17th May 2021
#
#@author: earjcti
#"""
#
#   This program will obtain the d18osw and d18op data from near Niels site and
#   write to a netcdf.  It will extract the data from HadCM3 experiments
#
#
# This program has been ammended from 
#extract_HadCM3_data_near_site
#
#

#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import iris
#import xlwt
#from xlwt import Workbook
import os
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from iris.cube import CubeList
import sys




    
def extract_lon_lat(cube):
    """
    extracts the points near LON_REQ and LAT_REQ
    """ 

    
    #cubelats = cube.coord('latitude').points
    #cubelons = cube.coord('longitude').points

    # find nearest latitude and lontiude to the value
    #latix = (np.abs(cubelats - LAT_REQ)).argmin()
    #lonix = (np.abs(cubelons - LON_REQ)).argmin()

    lonmin = LON_REQ - 10.0
    lonmax = LON_REQ + 10.0
    latmin = LAT_REQ - 10.0
    latmax = LAT_REQ + 10.0
    lon_constraint = iris.Constraint(longitude = lambda cell: 
                                     lonmin < cell < lonmax)
    lat_constraint = iris.Constraint(latitude = lambda cell: 
                                     latmin < cell < latmax)

    lon_slice  =  cube.extract(lon_constraint)
    data_slice = lon_slice.extract(lat_constraint)
   
    return data_slice

def extract_data():
    """
    extracts the data within ??deg of the required lat and long
    """
   
    allcubes = CubeList([])
  
    # this file was obtained from program regrid_HCM3_50_year_avg
    filename = ('/nfs/hera1/earjcti/um/' + model + '/' + FIELD + '/means/mean_month.nc')
    cube = iris.load_cube(filename)
    cube_slice = extract_lon_lat(cube)
    name = FIELD 
    cube_slice.long_name = name
    cube_slice.data = np.ma.where(cube_slice.data.mask == 1.0, -99999, cube_slice.data)
    cube_slice.data.mask = np.where(cube_slice.data == -99999, 1.0, 0.0)
    try:
        cube_slice.remove_coord('time')
    except:
        pass
    try:
        cube_slice.remove_coord('year')
    except:
        pass
    cube_slice.cell_methods=None
                            
    
    return cube_slice  
        


#################
# MAIN PROGRAM
################

###################################
# get initial data including the lats and longs we require

period = {'xozzb':'EOI400','xozza':'PI'}
linuxwin = 'l'
LAT_REQ = 52.5 # 51.2N, 4.4E
LON_REQ = 3.0
FIELD = 'd18osw'  # field = d18osw

model = 'xozzb'

alldata = extract_data()
fileoutstart = '/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/Niels/'
iris.save(alldata,fileoutstart  + FIELD +'_'+ period.get(model) + '_'+ model + '.nc', netcdf_format = "NETCDF3_CLASSIC",fill_value = -99999)
