#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#Created on Fri Jul  5 15:11:26 2019
#Updated for Friso September 2023
#
#@author: earjcti
#"""
#
#   This program will obtain the SST temperature and Salinity
#   from Frisos site ODP1209  long 158.506E, lat 32.65N 2387.4m
#                    ODP1208  lon 158.2E, 36.1N, depth 3345.7
#
#
# This program has been ammended from 
# Niels program
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

    lonmin = LON_REQ - 0.7
    lonmax = LON_REQ + 0.7
    latmin = LAT_REQ - 0.7
    latmax = LAT_REQ + 0.7
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
   
    eoi400_alltempcubes = CubeList([])
    e280_alltempcubes = CubeList([])
    eoi400_allsalcubes = CubeList([])
    e280_allsalcubes = CubeList([])

    eoi400_temp_filename = {
        'HadCM3': '/nfs/a103/palaeo_share/PlioMIP2/processed/Eoi400_2450-2499_Annual.nc',  
        'CESM1.2':'/nfs/b0164/Data/NCAR/b.e12.B1850.f09_g16.PMIP4-pliomip2.pop.h.TEMP.1100-1199.annavg.nc'}

    e280_temp_filename = {
        'HadCM3': '/nfs/a103/palaeo_share/PlioMIP2/processed/Preind_E280_2950-2999_Annual.nc',  
        'CESM1.2':'/nfs/b0164/Data/NCAR/b.e12.B1850.f09_g16.preind.pop.h.TEMP.0707-0806.annavg.nc'}

    eoi400_sal_filename = {
        'HadCM3': '/nfs/a103/palaeo_share/PlioMIP2/processed/Eoi400_2450-2499_Annual.nc',  
        'CESM1.2':'/nfs/b0164/Data/NCAR/b.e12.B1850.f09_g16.PMIP4-pliomip2.pop.h.SALT.1100-1199.annavg.nc'}
    
    e280_sal_filename = {
        'HadCM3': '/nfs/a103/palaeo_share/PlioMIP2/processed/Preind_E280_2950-2999_Annual.nc',  
        'CESM1.2':'/nfs/b0164/Data/NCAR/b.e12.B1850.f09_g16.preind.pop.h.SALT.0707-0806.annavg.nc'}

   
    tempname = {'HadCM3' : 'TEMPERATURE (OCEAN) DEG.C',
                'CESM1.2' : 'Potential Temperature'}
    salname = {'HadCM3' : 'SALINITY (OCEAN) (PSU)',
               'CESM1.2': 'Salinity'}
    
    for i, model in enumerate(MODELNAMES):
        print(eoi400_temp_filename.get(model),tempname.get(model))
        eoi400_temp_cube = iris.load_cube(eoi400_temp_filename.get(model),tempname.get(model))
        eoi400_sal_cube = iris.load_cube(eoi400_sal_filename.get(model),salname.get(model))
        eoi400_temp_cube_slice = extract_lon_lat(eoi400_temp_cube)
        eoi400_temp_cube_slice.long_name = eoi400_temp_cube_slice.long_name + '_EOI400_' + model

        eoi400_sal_cube_slice = extract_lon_lat(eoi400_sal_cube)
        eoi400_sal_cube_slice.long_name = eoi400_sal_cube_slice.long_name + '_EOI400_' + model
        eoi400_alltempcubes.append(eoi400_temp_cube_slice)
 
        eoi400_allsalcubes.append(eoi400_sal_cube_slice)


        e280_temp_cube = iris.load_cube(e280_temp_filename.get(model), tempname.get(model))
        e280_sal_cube = iris.load_cube(e280_sal_filename.get(model), salname.get(model))
        e280_temp_cube_slice = extract_lon_lat(e280_temp_cube)
        e280_temp_cube_slice.long_name = e280_temp_cube_slice.long_name + '_E280_' + model 

        e280_sal_cube_slice = extract_lon_lat(e280_sal_cube)
        e280_sal_cube_slice.long_name = e280_sal_cube_slice.long_name + '_E280_' + model 

        e280_alltempcubes.append(e280_temp_cube_slice)
        e280_allsalcubes.append(e280_sal_cube_slice)

                            
    
    return (eoi400_alltempcubes, e280_alltempcubes,
            eoi400_allsalcubes,e280_allsalcubes)
        


#################
# MAIN PROGRAM
################

###################################
# get initial data including the lats and longs we require

linuxwin = 'l'
#LAT_REQ = 32.65 # 51.2N, 4.4E
#LON_REQ = 158.506
LAT_REQ = 36.12 # 51.2N, 4.4E
LON_REQ = 158.2
test ='y'
#FIELD = 'NearSurfaceTemperature'
FIELDS = ['Temp','Salinity']
#FIELD = 'Salinity'

if test == 'y':
     MODELNAMES = ['HadCM3']
else:
     MODELNAMES = ['CCSM4', 'CCSM4-UoT', 'CCSM4-Utr',  'CESM1.2','CESM2',
                   'COSMOS', 'EC-Earth3.3', 'GISS2.1G', 'HadCM3','HadGEM3',
                   'IPSLCM6A', 'IPSLCM5A2', 'IPSLCM5A', 'MIROC4m', 'MRI2.3',
                   'NorESM-L', 'NorESM1-F'
                               ]
(eoi400_tempdata, e280_tempdata,
 eoi400_saldata, e280_saldata)= extract_data()

print(eoi400_tempdata, e280_tempdata, eoi400_saldata, e280_saldata)
cubelist = CubeList([eoi400_tempdata[0], e280_tempdata[0], eoi400_saldata[0], e280_saldata[0]])
print(cubelist)

iris.save(cubelist, 'tempfile.nc')
