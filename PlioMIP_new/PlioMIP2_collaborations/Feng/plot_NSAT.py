#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#5.10.2023
#@author: earjcti
#"""
#
# This is a sample program which will plot the PlioMIP2 MMM NSAT
#

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import iris
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
   
    timeslice = {'xozzb':'KM5c','xozzc':'K1','xozzd':'G17','xozze':'KM3',
                 'xozzf':'3.053Ma'}
    orbital_allcubes = CubeList([])
  
    for i, model in enumerate(MODELNAMES):
        filename = ('/nfs/hera1/earjcti/um/' + model + '/' + FIELD + '/means/mean_month.nc')
        cube = iris.load_cube(filename)
        cube_slice = extract_lon_lat(cube)
        name = timeslice.get(model) + '_' + FIELD 
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
        orbital_allcubes.append(cube_slice)
                            
    
    return orbital_allcubes  
        


#################
# MAIN PROGRAM
################


FIELD = 'NearSurfaceTemperature'



orbital_alldata = extract_data()
fileoutstart = '/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/Niels/'
iris.save(orbital_alldata,fileoutstart + 'orbital_' + FIELD + '.nc', netcdf_format = "NETCDF3_CLASSIC",fill_value = -99999)
