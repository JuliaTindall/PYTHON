#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created January 2021 by Julia

This program will create the input files for biome4 for each 'regridded' model or the multimodel mean

This was altered in November 2021 to do BIOME4 for a HadCM3 sensitivity study

"""

import numpy as np
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import sys
import warnings
warnings.filterwarnings("ignore")



def get_cube_avg(expt, long_name, units, long_name_req):
    """
    loads in renames and reformats that cube
    """
    monthnames = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
#    monthnames = ['ja','fb']

    allmonths_cubes = iris.cube.CubeList([])
    for i, month in enumerate(monthnames):
        cubes = iris.cube.CubeList([])

        for year in range(STARTYEAR,ENDYEAR):
            filename = (FILESTART + expt + '/netcdf/' + expt + 
                        'a@pdk' + np.str(year) + month + '.nc')
            indcube = iris.load_cube(filename, long_name)
            indcube.coord('t').points = year
           
            cubes.append(indcube)
        iris.util.equalise_attributes(cubes)
        iris.util.unify_time_units(cubes)
        cube = cubes.concatenate_cube()
        monthcube = cube.collapsed('t', iris.analysis.MEAN)
        
        if monthcube.coord('t') != None:
            monthcube.remove_coord('t')
        nmoncube = iris.util.new_axis(monthcube)
        nmoncube.add_dim_coord(iris.coords.DimCoord(i+1, 
                standard_name='time', long_name='t', 
                var_name='t', 
                units=None,
                bounds=None,
                coord_system=None, circular=False),0)
       
        nmoncube.cell_methods=None
       
        allmonths_cubes.append(nmoncube)

    iris.util.equalise_attributes(allmonths_cubes)
    iris.util.unify_time_units(allmonths_cubes)
  
    monthscube = allmonths_cubes.concatenate_cube()
    monthscube.coord('longitude').rename('lon')
    monthscube.coord('latitude').rename('lat')
    #'monthscube.coord('t').rename('time')
    monthscube.long_name = long_name_req
    monthscube.units = units
  
    return monthscube

def get_anom(field, plio_cube, pi_cube):
    """
    gets the standard data from biome 4 and regrids it onto our grid
    """
    filename = ('/nfs/see-fs-02_users/earjcti/BIOME4/'+
                'biome4_pliomip2/inputdata.nc')
    cube = iris.load_cube(filename, field)
    cubegrid = iris.load_cube('/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrparm.mask.nc', 'LAND MASK (LOGICAL: LAND=TRUE)')
    
    biomecube = cube.regrid(cubegrid, iris.analysis.Nearest())
    iris.util.squeeze(biomecube)
    sq_plio = iris.util.squeeze(plio_cube)
    sq_pi = iris.util.squeeze(pi_cube)
    biomecube_data = biomecube.data
    
    anomdata = sq_plio.data - sq_pi.data + biomecube_data  
    if field[0:4] != 'soil':
        anom_data = np.ma.masked_array(anomdata, biomecube.data.mask) 
    else:
        print('here')
        cubes = iris.load(filename)
        print(cubes)
        cubetemp = iris.load_cube(filename, 
                                  'annual absolute mimimum temperature')
        biomecube2 = cubetemp.regrid(cubegrid, iris.analysis.Nearest())
        biomearr3 = np.zeros((2,) +  np.shape(biomecube2.data))
        biomearr3[0, :, :] = biomecube2.data.mask
        biomearr3[1, :, :] = biomecube2.data.mask
        anom_datat = np.where(biomecube_data >= 0, anomdata, sq_plio.data)
        anom_data = np.ma.masked_array(anom_datat, biomearr3)          
    

    
    if (field == 'monthly total precipitation'
        or field ==  'mean monthly percent of possible sunshine'):
        anom_data = np.ma.where(anom_data > 0, anom_data, 0.0)
    if (field ==  'mean monthly percent of possible sunshine'):
        anom_data = np.ma.where(anom_data < 1000, anom_data, 1000.0)
   
    mpwp_anomcube = sq_plio.copy(data=anom_data)
    
    return mpwp_anomcube

    

def get_temp_precip():
    """
    get temperature and precipitation data in correct units
    plio_cube is from pliocene data
    anom_cube is model_plio - model_pi + observed_pi
    """
    allcubes_plio = iris.cube.CubeList([])
    allcubes_anom = iris.cube.CubeList([])

    # temperature
    cube = get_cube_avg(EXPTNAME, 'TEMPERATURE AT 1.5M', 'degC', 
                        'monthly mean temperature')
       
    cube.data = cube.data * 10.0
    allcubes_plio.append(cube)

    pi_cube = get_cube_avg(PI_EXPT, 'TEMPERATURE AT 1.5M', 'degC', 
                           'monthly mean temperature')
    pi_cube.data = pi_cube.data * 10.0

    anom_cube = get_anom('monthly mean temperature', cube, pi_cube)
    allcubes_anom.append(anom_cube)
  
    # get minimum temperature from the pliocene and the anomaly method
    print(cube)
    mintemp = cube.collapsed('time', iris.analysis.MIN)
    mintemp.long_name = 'annual absolute minimum temperature'
    mintemp.short_name = 'tmin'
    mintemp.units = 'degC'
    mintemp.remove_coord('time')
    #mintemp.remove_coord('surface')

    mintemp_anom = anom_cube.collapsed('time', iris.analysis.MIN)
    mintemp_anom.long_name = 'annual absolute minimum temperature'
    mintemp_anom.short_name = 'tmin'
    mintemp_anom.units = 'degC'
    mintemp_anom.remove_coord('time')
    #mintemp_anom.remove_coord('surface')
   


    # precipitation
    cube = get_cube_avg(EXPTNAME,'TOTAL PRECIPITATION RATE     KG/M2/S', 'mm',
                         'monthly total precipitation')
    cube.data = cube.data * 30.0
    allcubes_plio.append(cube)

    pi_cube = get_cube_avg(PI_EXPT,'TOTAL PRECIPITATION RATE     KG/M2/S', 'mm',
                         'monthly total precipitation')
    pi_cube.data = pi_cube.data * 30.0

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

    # mPWP
    cube_mPWP = get_cube_avg(EXPTNAME, 'TOTAL CLOUD AMOUNT - RANDOM OVERLAP',
                             'percent', 
                             'mean monthly percent of possible sunshine')
    cube_mPWP.short_name = 'sun'
    cube_mPWP.data = cube_mPWP.data / 100.
    cube_mPWP.data = (cube_mPWP.data * (-1.0) + 1.0) * 1000.

    # PI
    cube_PI = get_cube_avg(PI_EXPT, 'TOTAL CLOUD AMOUNT - RANDOM OVERLAP',
                           'percent', 
                           'mean monthly percent of possible sunshine')
    cube_PI.data = cube_PI.data / 100.
  
    cube_PI.short_name = 'sun'
    cube_PI.data = (cube_PI.data * (-1.0) + 1.0) * 1000.

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

        cube = iris.load_cube(filename)
        cube.attributes["missing_value"] = -9999

        # temporarily change sea points to average because we
        # will be adding a lsm later

        cube_data = cube.data
        avg_data = np.mean(cube_data)
        newcube_data = np.where(cube_data.mask, avg_data, cube_data)
        cube_full = cube.copy(data = newcube_data) 

        return cube_full

    
    # get grid
    cubegrid = iris.load_cube('/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrparm.mask.nc', 'LAND MASK (LOGICAL: LAND=TRUE)')
    

    # get mPWP

    file_soils_mPWP = ('/nfs/hera2/scripts/BIOME4/reference/' + 
                       'PRISM3_soil_alternative_whc.nc')
    soil_whc_cube_mPWP = process_soils(file_soils_mPWP)

    file_perc_mPWP = ('/nfs/hera2/scripts/BIOME4/reference/' + 
                      'PRISM3_soil_alternative_perc.nc')
    soil_perc_cube_mPWP = process_soils(file_perc_mPWP)
  


    # get PI

    file_soils_PI = ('/nfs/hera2/scripts/BIOME4/reference/' + 
                       'MODERN_soil_alternative_whc.nc')
    soil_whc_cube_PI = process_soils(file_soils_mPWP)

    file_perc_PI = ('/nfs/hera2/scripts/BIOME4/reference/' + 
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
    cubegrid = iris.load_cube('/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrparm.mask.nc', 'LAND MASK (LOGICAL: LAND=TRUE)')
   

    lsmcube = lsmcube_temp.regrid(cubegrid, iris.analysis.Linear())

    
    for cube in cubelist:
        print(np.shape(cube.data),np.shape(lsmcube.data))
        newcube_data = np.ma.where(lsmcube.data == 0, -9999, cube.data)
        newcube = cube.copy(data = newcube_data)
        masked_cubelist.append(newcube)
        print(newcube.data)
            
   
    return masked_cubelist


def main():
    """
    driver for program to get biome4 intput field
    """

    print('GETTING TEMP AND PRECIP')
    (allcubes_plio, abs_min_cube_plio,
     allcubes_anom, abs_min_cube_anom) = get_temp_precip()

    print('GETTING SUNSHINE')
    suncube_plio, suncube_anom = get_sunshine()
    allcubes_plio.append(suncube_plio)
    allcubes_anom.append(suncube_anom)

    print('GETTING SOIL')
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
#    allcubes_land_plio = apply_lsm(allcubes_plio)
  
    print('saving',OUTFILE)
    iris.save(allcubes_plio,OUTFILE + 'absolute.nc', 
              netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  


    allcubes_land_anom = apply_lsm(allcubes_anom)
    iris.save(allcubes_anom,OUTFILE + 'anomaly.nc', 
              netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  
    #iris.save(allcubes_anom,OUTFILE + 'allanomaly.nc', 
    #          netcdf_format="NETCDF3_CLASSIC", fill_value = -9999)  
    
 

FILESTART = '/nfs/hera1/earjcti/um/'
EXPTNAME = 'xozzh'
PI_EXPT = 'xozza'
STARTYEAR = 40
ENDYEAR=60
OUTFILE = FILESTART + EXPTNAME + '/biome4/inputdata_'

main()
