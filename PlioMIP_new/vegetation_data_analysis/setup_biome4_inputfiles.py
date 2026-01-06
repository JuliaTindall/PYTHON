#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created January 2021 by Julia

This program will create the input files for biome4 for each 'regridded' model or the multimodel mean

ammended 2025.  we also want to be able to do a preindustrial

"""

import numpy as np
import iris
from iris.cube import CubeList
import sys


def cube_reformat(fileend, long_name, units):
    """
    loads in renames and reformats that cube
    """
    filename = (FILESTART + MODELNAME + fileend)
    cube = iris.load_cube(filename)
    cube.remove_coord('time')
    cube.remove_coord('year')
    cube.coord('longitude').rename('lon')
    cube.coord('latitude').rename('lat')
    cube.coord('month').rename('time')
    cube.long_name = long_name
    cube.units = units
    return cube

def get_anom(field, plio_cube, pi_cube):
    """
    gets the standard data from biome 4 and regrids it onto our grid
    """
    filename = ('~earjcti/BIOME4/'+
                'biome4_pliomip2/inputdata.nc')
    cube = iris.load_cube(filename, field)
    cubegrid = iris.load_cube(FILESTART + MODELNAME + '/E280.NearSurfaceTemperature.allmean.nc')
    
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
    """
    allcubes_plio = CubeList([])
    allcubes_anom = CubeList([])
    allcubes_pi = CubeList([])

    # temperature
    cube = cube_reformat('/EOI400.NearSurfaceTemperature.mean_month.nc',
                         'monthly mean temperature', 'degC')
    cube.data = cube.data * 10.0
    allcubes_plio.append(cube)

    pi_cube = cube_reformat('/E280.NearSurfaceTemperature.mean_month.nc',
                         'monthly mean temperature', 'degC')
    pi_cube.data = pi_cube.data * 10.0
    allcubes_pi.append(pi_cube)


    anom_cube = get_anom('monthly mean temperature', cube, pi_cube)
    allcubes_anom.append(anom_cube)


    # get minimum temperature from the pliocene and the anomaly method
    mintemp = cube.collapsed('time', iris.analysis.MIN)
    mintemp.long_name = 'annual absolute minimum temperature'
    mintemp.short_name = 'tmin'
    mintemp.units = 'degC'
    mintemp.remove_coord('time')
    #mintemp.remove_coord('surface')

    mintemp_pi = pi_cube.collapsed('time', iris.analysis.MIN)
    mintemp_pi.long_name = 'annual absolute minimum temperature'
    mintemp_pi.short_name = 'tmin'
    mintemp_pi.units = 'degC'
    mintemp_pi.remove_coord('time')
  
    
    mintemp_anom = anom_cube.collapsed('time', iris.analysis.MIN)
    mintemp_anom.long_name = 'annual absolute minimum temperature'
    mintemp_anom.short_name = 'tmin'
    mintemp_anom.units = 'degC'
    mintemp_anom.remove_coord('time')
    #mintemp_anom.remove_coord('surface')
   


    # precipitation
    cube = cube_reformat('/EOI400.TotalPrecipitation.mean_month.nc',
                         'monthly total precipitation', 'mm')
    cube.data = cube.data * 30.0
    allcubes_plio.append(cube)

    pi_cube = cube_reformat('/E280.TotalPrecipitation.mean_month.nc',
                         'monthly total precipitation', 'degC')
    pi_cube.data = pi_cube.data * 30.0
    allcubes_pi.append(pi_cube)


    anom_cube = get_anom('monthly total precipitation', cube, pi_cube)
    allcubes_anom.append(anom_cube)

       
    return (allcubes_plio, mintemp, allcubes_anom, mintemp_anom,
            allcubes_pi,mintemp_pi)
    
def get_sunshine():
    """
    gets mean monthly percent of possible sunshine.  Steve P and James did
    this by:
    1. get pd total cloud (this is field30)
    2. multiplies this by -1 to get negative total cloud
    3. adds 1 to get total sun.
    4. multiplies by 1000 to get biome units.
    """

    # mPWP
    cube_mPWP = cube_reformat('/EOI400.totcloud.mean_month.nc',
                              'mean monthly percent of possible sunshine',
                              'percent')
    cube_mPWP.short_name = 'sun'
    if MODELNAME != 'HadCM3':
        cube_mPWP.data = cube_mPWP.data / 100.
    cube_mPWP.data = (cube_mPWP.data * (-1.0) + 1.0) * 1000.

    # PI
    cube_PI = cube_reformat('/E280.totcloud.mean_month.nc',
                              'mean monthly percent of possible sunshine',
                              'percent')
    if MODELNAME != 'HadCM3':
        cube_PI.data = cube_PI.data / 100.
  
    cube_PI.short_name = 'sun'
    cube_PI.data = (cube_PI.data * (-1.0) + 1.0) * 1000.

    # find anomaly from observations
    cube_anom = get_anom('mean monthly percent of possible sunshine',
                         cube_mPWP, cube_PI)
   

    return cube_mPWP, cube_anom, cube_PI
  
    

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
    cubegrid = iris.load_cube(FILESTART + MODELNAME + '/E280.NearSurfaceTemperature.allmean.nc')
   
    # get mPWP

    file_soils_mPWP = ('/uolstore/Research/a/hera2/scripts/BIOME4/reference_HadGEM/' + 
                       'PRISM3_soil_alternative_whc.nc')
    soil_whc_cube_mPWP = process_soils(file_soils_mPWP)

    file_perc_mPWP = ('/uolstore/Research/a/hera2/scripts/BIOME4/reference_HadGEM/' + 
                      'PRISM3_soil_alternative_perc.nc')
    soil_perc_cube_mPWP = process_soils(file_perc_mPWP)
  

    # get PI

    file_soils_PI = ('/uolstore/Research/ahera2/scripts/BIOME4/reference_HadGEM/' + 
                       'MODERN_soil_alternative_whc.nc')
    soil_whc_cube_PI = process_soils(file_soils_mPWP)

    file_perc_PI = ('/uolstore/Research/ahera2/scripts/BIOME4/reference_HadGEM/' + 
                      'MODERN_soil_alternative_perc.nc')
    soil_perc_cube_PI = process_soils(file_perc_mPWP)
  
    # get anomaly

    whc_anom_cube = get_anom('soil water holding capacity',
                         soil_whc_cube_mPWP, soil_whc_cube_PI)
    
    perc_anom_cube = get_anom('soil water percolation index',
                         soil_perc_cube_mPWP, soil_perc_cube_PI)

    return (soil_whc_cube_mPWP, soil_perc_cube_mPWP, 
            whc_anom_cube, perc_anom_cube,
            soil_whc_cube_PI, soil_perc_cube_PI)

def apply_lsm(cubelist,PLIO_PI_IND):
    """
    apply a lsm to the cubes and regrid
    """

# NOTE THE MASK IS NOT ON THE SAME GRID AS THE DATA

    masked_cubelist = CubeList([])
    if PLIO_PI_IND == 'Plio':
        lsm = '/uolstore/Research/a/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    if PLIO_PI_IND == 'PI':
        lsm = '/uolstore/Research/a/hera1/earjcti/PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc'
    
    lsmcube_temp = iris.load_cube(lsm)
    

    #cubegrid = iris.load_cube('/uolstore/Research/a/hera1/earjcti/regridded/' + MODELNAME + '/E280.NearSurfaceTemperature.allmean.nc')
    
    cubegrid=lsmcube_temp
    lsmcube = lsmcube_temp
    #lsmcube = lsmcube_temp.regrid(cubegrid, iris.analysis.Linear())

    print(cubegrid)

   
    for cube in cubelist:
        print(cube.name)
        try:
            cube.coord('lat').rename('latitude')
            cube.coord('lon').rename('longitude')
            rename='y'
        except:
            rename = 'n'
    
        if cube.ndim == 2:
            cube = cube.regrid(cubegrid,iris.analysis.Linear())
        if cube.ndim == 3:
            allsmall_cubelist = CubeList([])
            for i in range(0,len(cube.data[:,0,0])):
                smallcube = cube[i,:,:]
                smallcube_r = smallcube.regrid(cubegrid,iris.analysis.Linear())
                allsmall_cubelist.append(smallcube_r)
            cube = allsmall_cubelist.merge_cube()
            if rename =='y':
                cube.coord('latitude').rename('lat')
                cube.coord('longitude').rename('lon')
        newcube_data = np.ma.where(lsmcube.data == 0, -9999, cube.data)
        newcube = cube.copy(data = newcube_data)
        masked_cubelist.append(newcube)
            
   
    return masked_cubelist


def main():
    """
    driver for program to get biome4 intput field
    """
    
    (allcubes_plio, abs_min_cube_plio,
     allcubes_anom, abs_min_cube_anom,
     allcubes_pi,   abs_min_cube_pi) = get_temp_precip()
    
    suncube_plio, suncube_anom, suncube_pi = get_sunshine()
    allcubes_plio.append(suncube_plio)
    allcubes_anom.append(suncube_anom)
    allcubes_pi.append(suncube_pi)
    
    (soil_whc_plio, soil_perc_plio,
     soil_whc_anom, soil_perc_anom,
     soil_whc_pi, soil_perc_pi)= get_soils()

    allcubes_plio.append(soil_whc_plio)
    allcubes_plio.append(soil_perc_plio)
    allcubes_anom.append(soil_whc_anom)
    allcubes_anom.append(soil_perc_anom)
    allcubes_pi.append(soil_whc_pi)
    allcubes_pi.append(soil_perc_pi)

    allcubes_plio.append(abs_min_cube_plio)
    allcubes_anom.append(abs_min_cube_anom)
    allcubes_pi.append(abs_min_cube_pi)
                                  # following SPickering this is monthly
                                  # mean minimum not absolute minimum
                                  # best to run in anomaly mode.
    allcubes_land_plio = apply_lsm(allcubes_plio,'Plio')
  
    iris.save(allcubes_land_plio,OUTFILE + 'absolute.nc', 
              netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  


    allcubes_land_anom = apply_lsm(allcubes_anom,'Plio')
    iris.save(allcubes_land_anom,OUTFILE + 'anomaly.nc', 
              netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  

    allcubes_land_anom = apply_lsm(allcubes_pi,'PI')
    iris.save(allcubes_land_anom,OUTFILE + 'pi.nc', 
              netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  

  
FILESTART = '/uolstore/Research/a/hera1/earjcti/regridded/'
MODELNAME = 'IPSLCM6A'
OUTFILE = FILESTART + MODELNAME + '/biome4/inputdata_'

main()
