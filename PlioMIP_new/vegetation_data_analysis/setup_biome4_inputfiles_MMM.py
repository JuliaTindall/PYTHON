#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created January 2021 by Julia

This program will create the input files for biome4 for each 'regridded' model or the multimodel mean

"""

import numpy as np
from iris.experimental.equalise_cubes import equalise_attributes
import iris
import sys


def cube_reformat(fileend, field, long_name, units):
    """
    loads in renames and reformats that cube
    """
    print(FILESTART, fileend)
    filename = (FILESTART + fileend)
    print(filename, field)
    cube = iris.load_cube(filename, field)
    #cube.remove_coord('time')
    #cube.remove_coord('year')
    cube.coord('longitude').rename('lon')
    cube.coord('latitude').rename('lat')
    #cube.coord('month').rename('time')
    cube.long_name = long_name
    cube.short_name = None
    cube.var_name = None
    cube.standard_name = None
    cube.units = units
    return cube


def cube_cloud(modelnames, fileend,long_name, units):
    """
    loads in renames and reformats that cube
    """

  
    for i, model in enumerate(modelnames):
        filename = (FILESTART + model + fileend)
        cube = iris.load_cube(filename)
        cube.remove_coord('time')
        cube.remove_coord('year')
        cube.coord('longitude').rename('lon')
        cube.coord('latitude').rename('lat')
        cube.coord('month').rename('time')

        for coord in cube.coords():
            coord.points = coord.points.astype('float32')
            coord.bounds = None
            coord.attributes = None
        cube.data = cube.data.astype('float32')

        cube.long_name = long_name
        cube.var_name = None
        cube.standard_name = None
        cube.cell_methods = None
   
        cube.units = units
        if model != 'HadCM3':
            cube.data = cube.data / 100.
        cube.data = (cube.data * (-1.0) + 1.0) * 1000.
        cube.attributes = None

        if i==0:
            cubedata = cube.data
            count = 1
        else:
            cubedata = cubedata + cube.data
            count = count + 1
        
  
    meandata = cubedata / count
   
    meancube = cube.copy(data = meandata)
    meancube.short_name = 'sun'
      

    return meancube

def get_anom(field, plio_cube, pi_cube):
    """
    gets the standard data from biome 4 and regrids it onto our grid
    """
    filename = ('/nfs/see-fs-02_users/earjcti/BIOME4/'+
                'biome4_pliomip2/inputdata.nc')
    cube = iris.load_cube(filename, field)
    cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
    
    biomecube = cube.regrid(cubegrid, iris.analysis.Linear())
    biomecube_data = biomecube.data
    
    anom_data = np.where(biomecube_data.mask, 
                         plio_cube.data,
                         plio_cube.data - pi_cube.data + biomecube.data)

    if (field == 'monthly total precipitation'
        or field ==  'mean monthly percent of possible sunshine'):
        anom_data = np.where(anom_data > 0, anom_data, 0.0)
    if (field ==  'mean monthly percent of possible sunshine'):
        anom_data = np.where(anom_data < 1000, anom_data, 1000.0)
   
    mpwp_anomcube = plio_cube.copy(data=anom_data)
   
    return mpwp_anomcube

    

def get_temp_precip():
    """
    get temperature and precipitation data in correct units
    plio_cube is from pliocene data
    anom_cube is model_plio - model_pi + observed_pi
    everything is from the multimodel mean
    """
    allcubes_plio = iris.cube.CubeList([])
    allcubes_anom = iris.cube.CubeList([])

    # temperature
    cube = cube_reformat('/NearSurfaceTemperature_multimodelmean_month.nc',
                         'NearSurfaceTemperaturemean_plio', 
                         'monthly mean temperature','degC')
    cube.data = cube.data * 10.0
    cube.var_name = 'tas'
   
    allcubes_plio.append(cube)

    pi_cube = cube_reformat('/NearSurfaceTemperature_multimodelmean_month.nc',
                         'NearSurfaceTemperaturemean_pi', 
                          'monthly mean temperature','degC')
    pi_cube.data = pi_cube.data * 10.0
    pi_cube.var_name = 'tas'

    anom_cube = get_anom('monthly mean temperature', cube, pi_cube)
    allcubes_anom.append(anom_cube)


    # get minimum temperature from the pliocene and the anomaly method
    mintemp = cube.collapsed('time', iris.analysis.MIN)
    mintemp.long_name = 'annual absolute minimum temperature'
    mintemp.var_name = 'tas_0'
    mintemp.units = 'degC'
    mintemp.remove_coord('time')
    #mintemp.remove_coord('surface')

    mintemp_anom = anom_cube.collapsed('time', iris.analysis.MIN)
    mintemp_anom.long_name = 'annual absolute minimum temperature'
    mintemp_anom.var_name = 'tas_0'
    mintemp_anom.units = 'degC'
    mintemp_anom.remove_coord('time')
    #mintemp_anom.remove_coord('surface')
   


    # precipitation
    cube = cube_reformat('TotalPrecipitation_multimodelmean_month.nc',
                         'TotalPrecipitationmean_plio',
                         'monthly total precipitation', 'mm')
    cube.data = cube.data * 30.0
    cube.var_name = 'pr'
    allcubes_plio.append(cube)

    pi_cube = cube_reformat('TotalPrecipitation_multimodelmean_month.nc',
                         'TotalPrecipitationmean_pi',
                         'monthly total precipitation', 'mm')
    pi_cube.data = pi_cube.data * 30.0
    pi_cube.var_name = 'pr'

    anom_cube = get_anom('monthly total precipitation', cube, pi_cube)
    allcubes_anom.append(anom_cube)

       
    return allcubes_plio, mintemp, allcubes_anom, mintemp_anom
    
def get_sunshine():
    """
    gets mean monthly percent of possible sunshine.  Steve P and James did
    this by:
    1. get pd total cloud (this is field30)
    2. multiplies this by -1 to get negative total cloud
    3. adds 1 to get total sun.
    4. multiplies by 1000 to get biome units.
    """
    modelnames = ['CESM2', 'CCSM4-UoT', 'CCSM4-Utr', 'HadCM3',
                  'IPSLCM6A','MRI2.3']

    # mPWP
    cube_mPWP = cube_cloud(modelnames, '/EOI400.totcloud.mean_month.nc',
                           'mean monthly percent of possible sunshine',
                           'percent')
    cube_mPWP.var_name = 'clt'

   
    # PI
    cube_PI = cube_cloud(modelnames, '/E280.totcloud.mean_month.nc',
                            'mean monthly percent of possible sunshine',
                            'percent')
    cube_PI.var_name = 'clt'
 
    # find anomaly from observations
    cube_anom = get_anom('mean monthly percent of possible sunshine',
                         cube_mPWP, cube_PI)
   

    return cube_mPWP, cube_anom
  
    

def get_soils():
    """
    gets the soils for input to biome4
    input cubegrid: the grid to put the soils on
    """

    def process_soils(filename):
        """
        processes each file for the soils
        """

        cube1 = iris.load_cube(filename)
        cube = cube1.regrid(cubegrid, iris.analysis.Linear())
        cube.attributes["missing_value"] = -9999

        # temporarily change sea points to average because we
        # will be adding a lsm later

        cube_data = cube.data
        avg_data = np.mean(cube_data)
        newcube_data = np.where(cube_data.mask, avg_data, cube_data)
        cube_full = cube.copy(data = newcube_data) 

        return cube_full

    
    # get grid
    cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')

    # get mPWP

    file_soils_mPWP = ('/nfs/hera2/scripts/BIOME4/reference_HadGEM/' + 
                       'PRISM3_soil_alternative_whc.nc')
    soil_whc_cube_mPWP = process_soils(file_soils_mPWP)

    file_perc_mPWP = ('/nfs/hera2/scripts/BIOME4/reference_HadGEM/' + 
                      'PRISM3_soil_alternative_perc.nc')
    soil_perc_cube_mPWP = process_soils(file_perc_mPWP)
  

    # get PI

    file_soils_PI = ('/nfs/hera2/scripts/BIOME4/reference_HadGEM/' + 
                       'MODERN_soil_alternative_whc.nc')
    soil_whc_cube_PI = process_soils(file_soils_mPWP)

    file_perc_PI = ('/nfs/hera2/scripts/BIOME4/reference_HadGEM/' + 
                      'MODERN_soil_alternative_perc.nc')
    soil_perc_cube_PI = process_soils(file_perc_mPWP)
  
    # get anomaly

    whc_anom_cube = get_anom('soil water holding capacity',
                         soil_whc_cube_mPWP, soil_whc_cube_PI)
    
    perc_anom_cube = get_anom('soil water percolation index',
                         soil_perc_cube_mPWP, soil_perc_cube_PI)

    return (soil_whc_cube_mPWP, soil_perc_cube_mPWP, 
            whc_anom_cube, perc_anom_cube)

def apply_lsm(cubelist):
    """
    apply a lsm to the cubes
    """

# NOTE THE MASK IS NOT ON THE SAME GRID AS THE DATA

    masked_cubelist = iris.cube.CubeList([])
    lsm = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    lsmcube_temp = iris.load_cube(lsm)
    cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
    lsmcube = lsmcube_temp.regrid(cubegrid, iris.analysis.Linear())


   
    for cube in cubelist:
        newcube_data = np.ma.where(lsmcube.data == 0, -9999, cube.data)
        newcube = cube.copy(data = newcube_data)
        masked_cubelist.append(newcube)
        print(newcube.data)
            
   
    return masked_cubelist


def main():
    """
    driver for program to get biome4 intput field
    """

    
    (allcubes_plio, abs_min_cube_plio,
     allcubes_anom, abs_min_cube_anom) = get_temp_precip()
    
    suncube_plio, suncube_anom = get_sunshine()
    allcubes_plio.append(suncube_plio)
    allcubes_anom.append(suncube_anom)
    
    (soil_whc_plio, soil_perc_plio,
     soil_whc_anom, soil_perc_anom)= get_soils()

    allcubes_plio.append(soil_whc_plio)
    allcubes_plio.append(soil_perc_plio)
    allcubes_anom.append(soil_whc_anom)
    allcubes_anom.append(soil_perc_anom)

    allcubes_plio.append(abs_min_cube_plio)
    allcubes_anom.append(abs_min_cube_anom)
                                  # following SPickering this is monthly
                                  # mean minimum not absolute minimum
                                  # best to run in anomaly mode.
    allcubes_land_plio = apply_lsm(allcubes_plio)
  
    iris.save(allcubes_land_plio,OUTFILE + 'absolute.nc', 
              netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  


    allcubes_land_anom = apply_lsm(allcubes_anom)
    iris.save(allcubes_land_anom,OUTFILE + 'anomaly.nc', 
              netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  
    
 

FILESTART = '/nfs/hera1/earjcti/regridded/'
OUTFILE = FILESTART +  'BIOME4/inputdataMMM_'

main()
