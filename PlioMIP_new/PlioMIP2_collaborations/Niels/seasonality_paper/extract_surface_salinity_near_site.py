#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#Created on Fri Jul  5 15:11:26 2019
#Updated for JH on 17th May 2021
#
#@author: earjcti
#"""
#
#   This program will obtain the SST data from near Niels site and
#   write to a netcdf
#
#
# This program has been ammended from 
#PlioMIP2/large_scale_features/extract_data_locations.py
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
import iris.coord_categorisation
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
    if lonmin < 0.0: lonmin=0.
    lonmax = LON_REQ + 10.0
    latmin = LAT_REQ - 10.0
    latmax = LAT_REQ + 10.0
    print('j1',cube.coord('longitude').points,lonmin,lonmax)
    lon_constraint = iris.Constraint(longitude = lambda cell: 
                                     lonmin <= cell <= lonmax)
    #print('j2',lon_constraint)
    lat_constraint = iris.Constraint(latitude = lambda cell: 
                                     latmin <= cell <= latmax)
    #print('j3',lat_constraint)

    lon_slice  =  cube.extract(lon_constraint)
    #print('j4',lon_slice)
    data_slice = lon_slice.extract(lat_constraint)
    return data_slice

def extract_data():
    """
    extracts the data within ??deg of the required lat and long
    """
   
    eoi400_allcubes = CubeList([])
    e280_allcubes = CubeList([])
    anom_allcubes = CubeList([])
    fieldname = {'HadCM3':'SALINITY (OCEAN) (PSU)',
                 'COSMOS':'salinity','EC-Earth3.3':'sea_water_salinity',
                 'IPSL-CM6A-LR':'Sea Water Salinity',
                 'IPSLCM5A2':'sea_water_salinity',
                 'IPSLCM5A':'sea_water_salinity',
                 'MIROC4m':'so',}

    directory = '/nfs/b0164/Data/3D_ocean_salinity_Eoi400_E280/'

    for i, model in enumerate(MODELNAMES):

        if model in ['CCSM4','CESM1.2','CESM2','HadCM3',
                     'COSMOS','IPSLCM5A','IPSLCM5A2',
                     'IPSL-CM6A-LR','MIROC4m']:
            file_eoi400 = directory + 'Eoi400_' + model + '_so.nc'
            file_e280 = directory + 'E280_' + model + '_so.nc'
            eoi400_cube = iris.load_cube(file_eoi400,
                                         fieldname.get(model,'Salinity'))
            e280_cube = iris.load_cube(file_e280,fieldname.get(model,'Salinity'))
        elif model == 'CCSM4-UoT':
            # eoi400
            file_eoi400 = ('/nfs/b0164/Data/UofT/UofT-CCSM4/Eoi400/Omon/sos_Omon_UofT-CCSM4_midPliocene-eoi400_r1i1p1f1_gr1_160101-170012.nc')
            fullcube = iris.load_cube(file_eoi400,'Sea Surface Salinity')
            iris.coord_categorisation.add_month_number(fullcube,'time',name='month')
            eoi400_cube = fullcube.aggregated_by('month',iris.analysis.MEAN)
            iris.util.promote_aux_coord_to_dim_coord(eoi400_cube,'month')
            eoi400_cube.coord('longitude').rename('lonaux')
            points=eoi400_cube.coord('lonaux').points[0]
            eoi400_cube.add_dim_coord(iris.coords.DimCoord(points,
                 standard_name = 'longitude',  long_name = 'longitude',
                 var_name = 'longitude', units='degrees',
                 coord_system = None,  circular = True), 2)
            eoi400_cube.coord('latitude').rename('lataux')
            points=eoi400_cube.coord('lataux').points[:,0]
            eoi400_cube.add_dim_coord(iris.coords.DimCoord(points,
                 standard_name = 'latitude',  long_name = 'latitude',
                 var_name = 'latitude',
                 coord_system = None,  circular = False), 1)

            # e280
            file_e280 = ('/nfs/b0164/Data/UofT/UofT-CCSM4/E280/Omon/sos_Omon_UofT-CCSM4_piControl_r1i1p1f1_gr1_150101-160012.nc')
            fullcube = iris.load_cube(file_e280,'Sea Surface Salinity')
            iris.coord_categorisation.add_month_number(fullcube,'time',name='month')
            e280_cube = fullcube.aggregated_by('month',iris.analysis.MEAN)
            iris.util.promote_aux_coord_to_dim_coord(e280_cube,'month')
            e280_cube.coord('longitude').rename('lonaux')
            points=e280_cube.coord('lonaux').points[0]
            e280_cube.add_dim_coord(iris.coords.DimCoord(points,
                 standard_name = 'longitude',  long_name = 'longitude',
                 var_name = 'longitude', units='degrees',
                 coord_system = None,  circular = True), 2)
            e280_cube.coord('latitude').rename('lataux')
            points=e280_cube.coord('lataux').points[:,0]
            e280_cube.add_dim_coord(iris.coords.DimCoord(points,
                 standard_name = 'latitude',  long_name = 'latitude',
                 var_name = 'latitude',
                 coord_system = None,  circular = False), 1)



        elif model == 'HadGEM3':
            # eoi400
            file_eoi400 = ('/nfs/b0164/Data/HadGEM3_new/climatologies/Eoi400/ocean/clims_hadgem3_pliocene_sal_final.nc')
            eoi400_cube = iris.load_cube(file_eoi400,'sal')

            # e280
            file_e280 = ('/nfs/b0164/Data/HadGEM3_new/climatologies/E280/ocean/clims_hadgem3_pi_sal.nc')
            e280_cube = iris.load_cube(file_e280,'sal')

        elif model == 'NorESM-L':
            # eoi400
            file_eoi400 = ('/nfs/b0164/Data/'+model+'/' + model + '_Eoi400.sss.climo.nc')
            eoi400_cube = iris.load_cube(file_eoi400,'Ocean surface salinity')
            print('eoi400',model,eoi400_cube)

            # e280
            file_e280 = ('/nfs/b0164/Data/'+model+'/' + model + '_E280.sss.climo.nc')
            e280_cube = iris.load_cube(file_e280,'Ocean surface salinity')
            print('e280',model,e280_cube)

        elif model == 'NorESM1-F':
            # eoi400
            file_eoi400 = ('/nfs/b0164/Data/'+model+'/' + model + '_Eoi400_sss_climo.nc')
            eoi400_cube = iris.load_cube(file_eoi400,'Ocean surface salinity')
            print('eoi400',model,eoi400_cube)

            # e280
            file_e280 = ('/nfs/b0164/Data/'+model+'/' + model + '_E280_sss_climo.nc')
            e280_cube = iris.load_cube(file_e280,'Ocean surface salinity')
            print('e280',model,e280_cube)


        else:
            print('setup model',model)
            sys.exit(0)


        eoi400_cube_slice = extract_lon_lat(eoi400_cube)
        name = model + '_Salinity_EOI400'
        eoi400_cube_slice.long_name = name
        eoi400_cube_slice.data = np.ma.where(eoi400_cube_slice.data.mask == 1.0, -99999, eoi400_cube_slice.data)
        eoi400_cube_slice.data.mask = np.where(eoi400_cube_slice.data == -99999, 1.0, 0.0)
        eoi400_cube_slice.cell_methods=None
        eoi400_cube_slice.coord('time').bounds = None
        try:
            eoi400_cube_slice.coord('depth').bounds = None
        except:
            pass
        eoi400_allcubes.append(eoi400_cube_slice)
        

        e280_cube_slice = extract_lon_lat(e280_cube)
        name = model + '_Salinity_E280'
        e280_cube_slice.long_name = name
        e280_cube_slice.data = np.ma.where(e280_cube_slice.data.mask == 1.0, -99999, e280_cube_slice.data)
        e280_cube_slice.data.mask = np.where(e280_cube_slice.data == -99999, 1.0, 0.0)
        e280_cube_slice.cell_methods=None
        e280_cube_slice.coord('time').bounds = None
        try:
            e280_cube_slice.coord('depth').bounds = None
        except:
            pass
        e280_allcubes.append(e280_cube_slice)

         
        anom_data = eoi400_cube_slice.data - e280_cube_slice.data
        anom_data.mask = np.maximum(eoi400_cube_slice.data.mask,
                                e280_cube_slice.data.mask)
        anom_data = np.where(anom_data.mask == 1.0, -99999, anom_data)
        anom_cube = e280_cube_slice.copy(data=anom_data)
        name = model + '_' + FIELD + '_EOI400-E280'
        anom_cube.long_name = name
      
        anom_allcubes.append(anom_cube)
                            
    
    return eoi400_allcubes, e280_allcubes,anom_allcubes    
        


#################
# MAIN PROGRAM
################

###################################
# get initial data including the lats and longs we require

linuxwin = 'l'
LAT_REQ = 52.5 # 51.2N, 4.4E
LON_REQ = 3.0
test ='n'
FIELD='Salinity'
if test == 'y':
     MODELNAMES = ['HadGEM3','HadCM3','NorESM-L','NorESM1-F']
else:
     MODELNAMES = ['COSMOS', 'HadCM3',
                   'IPSL-CM6A-LR', 'IPSLCM5A2', 'IPSLCM5A', 'MIROC4m', 
                    'NorESM1-F','CCSM4-UoT','HadGEM3','NorESM-L']
     julia = ['CCSM4-Utr','GISS2.1G',
                   'EC-Earth3.3',  'HadGEM3', 'MRI2.3']
eoi400_alldata, e280_alldata, anom_alldata = extract_data()
fileoutstart = '/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/Niels/'
iris.save(eoi400_alldata,fileoutstart + 'EOI400_' + FIELD + '.nc', netcdf_format = "NETCDF3_CLASSIC",fill_value = -99999)
iris.save(e280_alldata,fileoutstart + 'E280_' + FIELD + '.nc', netcdf_format = "NETCDF4", fill_value = -99999)
iris.save(anom_alldata,fileoutstart + 'EOI400-E280_anom_' + FIELD + '.nc', netcdf_format = "NETCDF4", fill_value=-99999)
