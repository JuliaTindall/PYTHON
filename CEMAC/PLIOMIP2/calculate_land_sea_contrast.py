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

    def get_ipsl_lsm(file, fieldnames):
        # get's the ipsl lsm which is sum of terrestrial and land ice
        cubes = iris.load(file, fieldnames)
        cube = cubes[0] + cubes[1]
        lsm_cube = cube.collapsed('time_counter', iris.analysis.MEAN)
        return lsm_cube

    def change_to_2d(cube):
        # if cube is 3d then extract the first time dimension only
        if cube.ndim == 2:
            cube_2d = cube
        else:
            cube_2d = cube[0, :, :]
       
        return cube_2d

    ############################################
    if MODELNAME == 'IPSLCM5A' or MODELNAME == 'IPSLCM5A2':
        plio_lsm_cube = get_ipsl_lsm(LSM_PLIO, FIELDLSM)
        pi_lsm_cube = get_ipsl_lsm(LSM_PI, FIELDLSM)
    elif MODELNAME == 'HadGEM3':
        test = iris.fileformats.netcdf.load_cubes(LSM_PLIO, callback=None)
        print(test)
        for data in test:
            print(data)
        sys.exit(0)
        f = netCDF4.Dataset(LSM_PLIO, "r")
        print(f)
        sys.exit(0)
    else:
        plio_lsm_cube = iris.util.squeeze(iris.load_cube(LSM_PI, FIELDLSM))
        pi_lsm_cube = iris.util.squeeze(iris.load_cube(LSM_PLIO, FIELDLSM))
   
    plio_lsm_cube2 = change_to_2d(plio_lsm_cube)
    pi_lsm_cube2 = change_to_2d(pi_lsm_cube)
   
    plio_lsm_data = plio_lsm_cube2.data
    pi_lsm_data = pi_lsm_cube2.data

    if MODELNAME == 'IPSLCM6A':
        plio_lsm_data = plio_lsm_data / 100.0
        pi_lsm_data = pi_lsm_data / 100.0

    if MODELNAME == 'EC-Earth3.3':
        plio_lsm_data = np.where(plio_lsm_data > 0.5, 1.0, 0.0)
        pi_lsm_data = np.where(pi_lsm_data > 0.5, 1.0, 0.0)
  
  
    land_mask = np.zeros(np.shape(plio_lsm_data))
    sea_mask = np.zeros(np.shape(plio_lsm_data))

    for ix, plio_mask in np.ndenumerate(plio_lsm_data):
        if plio_mask == 1.0 and pi_lsm_data[ix] == 1.0:
            land_mask[ix] = 1.0
        #if pi_lsm_data[ix] > 0:
        #    land_mask[ix] = 1.0
        if plio_mask == 0.0 and pi_lsm_data[ix] == 0.0:
            sea_mask[ix] = 1.0

    land_cube = plio_lsm_cube2.copy(data=land_mask)
    land_cube.var_name = 'land_mask'
    land_cube.long_name = 'land_mask'
    print('getting sea mask')
    sea_cube = pi_lsm_cube2.copy(data=sea_mask)
    sea_cube.var_name = 'sea_mask'
    sea_cube.long_name = 'sea_mask'
    print('got sea mask')

    return land_cube, sea_cube

def get_hadcm3_data(filestart):
    """
    gets the nsat data from HadCM3 and MRI
    called by get_nsat_data
    """
    allcubes = iris.cube.CubeList([])
    startyear = 0
    endyear = 100
    if MODELNAME  == 'MRI-CGCM2.3':
        startyear = startyear+1
        endyear = endyear+1

    for i in range(startyear, endyear):
        yearuse = str(i).zfill(3)
        filenameuse = (filestart + yearuse + '.nc')
        cubetemp = iris.load_cube(filenameuse)

        u  =  unit.Unit('days since 0800-01-01 00:00:00',
                                  calendar = unit.CALENDAR_360_DAY)
        if MODELNAME  == 'HadCM3':
            cubetemp.coord('t').rename('time')
        cubetemp.coord('time').points = (np.arange(0, 12)+((i-startyear)*12))*30.

        cubetemp.coord('time').units = u

        allcubes.append(cubetemp)


    equalise_attributes(allcubes)
    cube_temp = allcubes.concatenate_cube()

    if MODELNAME  == 'MRI-CGCM2.3':
        cube_temp.coord('pressure level').rename('surface')
  
    if MODELNAME  == 'HadCM3':
        cube_temp.coord('ht').rename('surface')

    cube_temp.coord('surface').points = 0.
    cube  =  cube_temp.extract(iris.Constraint(surface = 0.))

    return cube

def get_ipslcm6a_data(file):
    """
    for ipslcm6a we have 200 years in the file.  but we only need 100 years
    """
    cube = iris.load_cube(file, FIELDNAME)
    # reduce number of years
    cubelist = iris.cube.CubeList([])
    for i,  t_slice in enumerate(cube.slices(['latitude', 'longitude'])):
        if i >= 1200:
            t_slice.coord('time').bounds = None
            t_slice2 = iris.util.new_axis(t_slice, 'time')
            cubelist.append(t_slice2)

    cube100yr = cubelist.concatenate_cube()
  
    return cube100yr

def get_ipsl5_data(filename, exptname):
    """
    gets nsat data from ipsl
    there is a bit of an error in the file calendar so we will
    """
# copy the data to a new file but without the error
    with Dataset(filename) as src,  Dataset("temporary.nc",  "w", format = 'NETCDF3_CLASSIC') as dst:
        # copy attributes
        for name in src.ncattrs():
            dst.setncattr(name,  src.getncattr(name))
        # copy dimensions
        #for name,  dimension in src.dimensions.iteritems():
        for name,  dimension in src.dimensions.items():

            if name !=  'tbnds':   # don't copy across time counter bounds
                dst.createDimension(name,  (len(dimension)))

        # copy all file data
        for name,  variable in src.variables.items():
            if name != 'time_counter_bnds' and name!= 'time_centered':
                x  =  dst.createVariable(name,  variable.datatype,
                                       variable.dimensions)
                if name  == 'time_counter':
                    # convert from seconds to days and start at middle of month
                    dst.variables[name][:]  =  (src.variables[name][:] / (60.*60.*24))-(src.variables[name][0] / (60.*60.*24))+15.
                else:
                    dst.variables[name][:]  =  src.variables[name][:]
                # copy attributes for this variable
                for ncattr in src.variables[name].ncattrs():
                    attribute = src.variables[name].getncattr(ncattr)
                    if ncattr  == 'calendar' and exptname  == 'Eoi400':
                        dst.variables[name].setncattr(ncattr, '360_day')
                    else:
                        if (ncattr  == 'units' and name  == 'time_counter'):
                    # change units from seconds to days
                            dst.variables[name].setncattr(ncattr, attribute.replace('seconds', 'days'))
                        elif ncattr !='_FillValue':
                            dst.variables[name].setncattr(ncattr, attribute)

        fieldreq = 'Temperature 2m'

        cube = iris.load_cube('temporary.nc', fieldreq)

        cube.convert_units('Celsius')

        if exptname  == 'Eoi400':
            u  =  unit.Unit('days since 0800-01-01 00:00:00',
                              calendar = unit.CALENDAR_360_DAY)
        else:
            u  =  unit.Unit('days since 0800-01-01 00:00:00',
                          calendar = unit.CALENDAR_365_DAY)
        cube.coord('time').units = u

        return(cube)

def get_ipslcm5a2_data(filename):
    """
    gets the data for ipslcm5a2
    and removes all auxillary coordinates
    """
    cubelist = iris.cubeList
    cube = iris.load_cube(filename, FIELDNAME)
    for coord in cube.aux_coords:
        coord.rename('toremove')
        cube.remove_coord('toremove')
    return cube

def get_giss_data(filenames):
    """
    gets giss data: this is in two files
    """ 
    allcubes = iris.cube.CubeList([])
    for file in filenames:
        cubetemp = iris.load_cube(file, FIELDNAME)
        allcubes.append(cubetemp)

    equalise_attributes(allcubes)

    cube = allcubes.concatenate_cube()
  
    return(cube)
    
def get_nsat_data():
    """
    get the average temperature from the pliocene and the preindustrial
    """
        

    if MODELNAME == 'HadCM3' or MODELNAME == 'MRI-CGCM2.3':
        cube_plio = get_hadcm3_data(FILENAME_PLIO)
        cube_pi = get_hadcm3_data(FILENAME_PI)
    elif MODELNAME == 'IPSLCM6A':
        cube_plio = get_ipslcm6a_data(FILENAME_PLIO)
        cube_pi = get_ipslcm6a_data(FILENAME_PI)
    elif MODELNAME == 'IPSLCM5A':
        cube_pi = iris.load(FILENAME_PI)[0]
        cube_plio = get_ipsl5_data(FILENAME_PLIO,'Eoi400')
    elif MODELNAME == 'IPSLCM5A2':
        cube_plio = get_ipslcm5a2_data(FILENAME_PLIO)
        cube_pi = get_ipslcm5a2_data(FILENAME_PI)
    elif MODELNAME == 'GISS2.1G':
        cube_plio = get_giss_data(FILENAME_PLIO)
        cube_pi = get_giss_data(FILENAME_PI)
    else:
        cube_plio = iris.load_cube(FILENAME_PLIO, FIELDNAME)
        cube_pi = iris.load_cube(FILENAME_PI, FIELDNAME)
    
   
    cube_plio_avg = cube_plio.collapsed('time', iris.analysis.MEAN)
    cube_pi_avg = cube_pi.collapsed('time', iris.analysis.MEAN)
  
    return cube_plio_avg, cube_pi_avg


def get_global_avg(land_cube, sea_cube, cube_nsat):
    """
    get's global average temperature, and also global avg for
    the land and the ocean
    """


    if cube_nsat.coord('latitude').has_bounds():
        cube_nsat.coord('latitude').bounds
    else:
        cube_nsat.coord('latitude').guess_bounds()

    if cube_nsat.coord('longitude').has_bounds():
        cube_nsat.coord('longitude').bounds
    else:
        cube_nsat.coord('longitude').guess_bounds()

  
    grid_areas = iris.analysis.cartography.area_weights(cube_nsat)
    grid_areas_land = grid_areas * land_cube.data
    grid_areas_sea = grid_areas * sea_cube.data

    avg_cube = cube_nsat.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas)

    avg_cube_land = cube_nsat.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_land)

    avg_cube_sea = cube_nsat.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_sea)

    avg_temp_data = avg_cube.data
    avg_temp_land = avg_cube_land.data
    avg_temp_sea = avg_cube_sea.data

    if avg_temp_data > 200.:
        avg_temp_data = avg_temp_data - 273.15
        avg_temp_land = avg_temp_land -273.15
        avg_temp_sea = avg_temp_sea -273.15
   
    return avg_temp_data, avg_temp_land, avg_temp_sea, grid_areas


def get_regional_landsea(rmax, rmin, land_cube, sea_cube, cube_nsat,
                         grid_areas_region):
    """
    gets the mean temperature for latitude bands for average and for
    land and sea  
    """

    lats = cube_nsat.coord('latitude').points
    grid_areas_use = grid_areas_region * 1.0
    for j, lat in enumerate(lats):
        if rmin > lat or rmax < lat:
            grid_areas_use[j, :] = 0.0
    grid_areas_land = grid_areas_use * land_cube.data
    grid_areas_sea = grid_areas_use * sea_cube.data
    

    avg_cube = cube_nsat.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_use)
   
    avg_cube_land = cube_nsat.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_land)
   
    avg_cube_sea = cube_nsat.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_sea)

    avg_temp_data = avg_cube.data
    avg_temp_land = avg_cube_land.data
    avg_temp_sea = avg_cube_sea.data

    if avg_temp_data > 200.:
        avg_temp_data = avg_temp_data - 273.15
        avg_temp_land = avg_temp_land -273.15
        avg_temp_sea = avg_temp_sea -273.15

    return [avg_temp_data, avg_temp_land, avg_temp_sea]
           
                                                       
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
    print('lsm is', LSM_PLIO)

    # get land and sea mask
    land_mask_cube, sea_mask_cube = get_land_sea_mask()

    # get temporally averaged nsat data
    print('getting temporally averaged nsat data')
    cube_nsat_plio, cube_nsat_pi = get_nsat_data()
    

    # get land and sea temperatures
    print('getting land sea temperatures')
    (avg_T_plio, 
     avg_land_T_plio, 
     avg_sea_T_plio,
     gridareas) = get_global_avg(land_mask_cube, sea_mask_cube,
                                      cube_nsat_plio)

 

    (avg_T_pi, 
     avg_land_T_pi,
     avg_sea_T_pi,
     grid_areas) = get_global_avg(land_mask_cube, sea_mask_cube,
                                     cube_nsat_pi)

    avg_T_anom = avg_T_plio - avg_T_pi
    avg_T_landanom = avg_land_T_plio - avg_land_T_pi
    avg_T_seaanom = avg_sea_T_plio - avg_sea_T_pi



    regionmax = [90.0, 60.0, 30.0, 0.0, -30.0, -60.0]
    regionmin = [60.0, 30.0, 0.0, -30.0, -60.0, -90.0]
    land_sea_region_pi = np.zeros((3, len(regionmax)))
    land_sea_region_plio = np.zeros((3, len(regionmax)))

    #plt.subplot(2,1,1)
    #V = [0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2]
    #qplt.contourf(cube_nsat_plio - cube_nsat_pi, levels=V)
    #plt.subplot(2,1,2)
    #qplt.contourf((sea_mask_cube + land_mask_cube),  levels=V)
    #plt.show()
    #sys.exit(0)
   


    for i, rmax in enumerate(regionmax):
        land_sea_region_pi[:, i] = get_regional_landsea(rmax, regionmin[i],
                                                        land_mask_cube,
                                                        sea_mask_cube,
                                                        cube_nsat_pi,
                                                        grid_areas)
        land_sea_region_plio[:, i] = get_regional_landsea(rmax, regionmin[i],
                                                        land_mask_cube,
                                                        sea_mask_cube,
                                                        cube_nsat_plio,
                                                        grid_areas)

       
    
 # write to a spreadsheet
    land_sea_anom = land_sea_region_plio - land_sea_region_pi

    write_to_spreadsheet(avg_T_anom, avg_T_landanom, avg_T_seaanom,
                         land_sea_anom, regionmax, regionmin)
    


    

   

#############################################################################
def getnames():

# this program will get the names of the files and the field for each
# of the model

  
    # get names for each model

    if MODELNAME == 'CESM2':
        file_e280 = (FILESTART + 'NCAR/b.e21.B1850.' + 
                         'f09_g17.CMIP6-piControl.' + 
                         '001.cam.h0.TREFHT.110001-120012.nc')
        file_eoi400 = (FILESTART + 'NCAR/b.e21.B1850.' + 
                           'f09_g17.PMIP4-midPliocene-eoi400.' + 
                           '001.cam.h0.TREFHT.1101.1200.nc')
        fielduse = 'Reference height temperature'
        lsm_e280 = (FILESTART + 'NCAR/b.e12.B1850.' + 
                    'f09_g16.preind.cam.h0.LANDFRAC.0701.0800.nc')
        lsm_eoi400 = (FILESTART + 'NCAR/b.e21.B1850.' + 
                      'f09_g17.PMIP4-midPliocene-eoi400.001.' + 
                      'cam.h0.LANDFRAC.1101.1200.nc')
        fieldlsm = 'Fraction of sfc area covered by land'

    if MODELNAME == 'COSMOS':
        file_e280 = (FILESTART + 'AWI/COSMOS/E280/E280.tas'
                      '_2650-2749_monthly_mean_time_series.nc')
        file_eoi400 = (FILESTART + 'AWI/COSMOS/Eoi400/Eoi400.tas'
                      '_2650-2749_monthly_mean_time_series.nc')
        fielduse =  "2m temperature"
        lsm_e280 = ("/nfs/hera1/pliomip2/data/AWI/COSMOS/land_sea_masks/" + 
                    "E280_et_al/E280.slf.atm.nc")
        lsm_eoi400 = ("/nfs/hera1/pliomip2/data/AWI/COSMOS/land_sea_masks/" + 
                      "Eoi400_et_al/Eoi400.slf.atm.nc")
        fieldlsm = "SLF"

    if MODELNAME == 'EC-Earth3.3':
        file_e280 = FILESTART + 'EC-Earth3.3/EC-Earth3.3_PI_surface.nc'
        file_eoi400 = FILESTART + 'EC-Earth3.3/EC-Earth3.3_mPlio_surface.nc'
        fielduse = 'Air temperature at 2m'
        lsm_e280 =  FILESTART + 'EC-Earth3.3/EC-Earth3.3_PI_LSM.nc'
        lsm_eoi400 =  FILESTART + 'EC-Earth3.3/EC-Earth3.3_mPlio_LSM.nc'
        fieldlsm = 'Land/sea mask'

    if MODELNAME == 'CESM1.2':
        file_e280 = (FILESTART + 'NCAR/b.e12.B1850.' + 
                         'f09_g16.preind.cam.h0.TREFHT.0701.0800.nc')
        file_eoi400 = (FILESTART + 'NCAR/b.e12.B1850.' + 
                       'f09_g16.PMIP4-pliomip2.cam.h0.TREFHT.1101.1200.nc')
        fielduse = 'Reference height temperature'
        lsm_e280 = (FILESTART + 'NCAR/b.e12.B1850.' + 
                    'f09_g16.preind.cam.h0.LANDFRAC.0701.0800.nc')
        lsm_eoi400 = (FILESTART + 'NCAR/b.e12.B1850.' + 
                      'f09_g16.PMIP4-pliomip2.cam.h0.LANDFRAC.1101.1200.nc')
        fieldlsm = 'Fraction of sfc area covered by land'       
  
    if MODELNAME   ==  'MIROC4m':
        file_e280 = (FILESTART + 'MIROC4m/tas/MIROC4m_E280_Amon_tas.nc')
        file_eoi400 = (FILESTART + 'MIROC4m/tas/MIROC4m_Eoi400_Amon_tas.nc')
        fielduse = "tas"
        lsm_e280 = (FILESTART + 'MIROC4m/sftlf/MIROC4m_Exxx_fx_sftlf.nc')
        lsm_eoi400 = (FILESTART + 'MIROC4m/sftlf/MIROC4m_Eoixxx_fx_sftlf.nc')
        fieldlsm = "sftlf"

    if MODELNAME  == 'HadCM3':
        file_e280 = (FILESTART+'LEEDS/HadCM3/e280/NearSurfaceTemperature/' + 
                     'e280.NearSurfaceTemperature.')
        file_eoi400 = (FILESTART+'LEEDS/HadCM3/eoi400/NearSurfaceTemperature/' + 
                     'eoi400.NearSurfaceTemperature.')
        fielduse = "TEMPERATURE AT 1.5M"
        lsm_e280 = (FILESTART+'LEEDS/HadCM3/e280/qrparm.mask.nc')
        lsm_eoi400 = (FILESTART+'LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc')
        fieldlsm = 'LAND MASK (LOGICAL: LAND=TRUE)'

    if MODELNAME == 'CCSM4':
        file_e280 = (FILESTART + 'NCAR/b40.B1850.' + 
                         'f09_g16.preind.cam.h0.TREFHT.0081.0180.nc')
        file_eoi400 = (FILESTART + 'NCAR/b40.B1850.' + 
                       'f09_g16.PMIP4-pliomip2.TREFHT.1001.1100.nc')
        fielduse = 'Reference height temperature'
        lsm_e280 = (FILESTART + 'NCAR/b40.B1850.' + 
                    'f09_g16.preind.cam.h0.LANDFRAC.0081.0180.nc')
        lsm_eoi400 = (FILESTART + 'NCAR/b40.B1850.' + 
                      'f09_g16.PMIP4-pliomip2.LANDFRAC.1001.1100.nc')
        fieldlsm = 'Fraction of sfc area covered by land'

    if MODELNAME == 'CCSM4_Utr':
        file_e280 = (FILESTART + 'Utrecht/CESM1.0.5/E280/' + 
                     'tas_Amon_CESM1.0.5_E280_r1i1p1f1_gn_275001-285012.nc')  
        file_eoi400 = (FILESTART + 'Utrecht/CESM1.0.5/Eoi400/' +
                       'tas_Amon_CESM1.0.5_Eoi400_r1i1p1f1_gn_190001-200012.nc')
        fielduse = 'Reference height temperature'
        lsm_e280 = (FILESTART + 'Utrecht/CESM1.0.5/E280/' + 
                    'land_sea_mask_Amon_CESM1.0.5_b.PI_1pic_f19g16_NESSC' + 
                    '_control_r1i1p1f1_gn.nc')
        lsm_eoi400 = (FILESTART + 'Utrecht/CESM1.0.5/Eoi400/' +
                      'land_sea_mask_Amon_CESM1.0.5_b.PLIO_5Ma_Eoi400_' + 
                      'f19g16_NESSC_control_r1i1p1f1_gn.nc')
        fieldlsm = 'LANDMASK[D=1]'
  
    if MODELNAME == 'CCSM4_UoT':
        start = FILESTART + 'UofT/UofT-CCSM4/'
        file_e280 = (start + '/E280/Amon/native_grid/tas_Amon_' + 
                     'UofT-CCSM4_piControl_r1i1p1f1_gn_150101-160012.nc')  
        file_eoi400 = (start + '/Eoi400/Amon/native_grid/tas_Amon_' + 
                     'UofT-CCSM4_midPliocene-eoi400_r1i1p1f1_gn_' + 
                       '160101-170012.nc') 
        fielduse = 'air_temperature'
        lsm_e280 = start + 'for_julia/E_mask.nc'
        lsm_eoi400 = start + 'for_julia/Eoi_mask.nc'
        fieldlsm = 'gridbox land fraction'
      
    if MODELNAME == 'NorESM-L':
       file_e280 = (FILESTART + 'NorESM-L/NorESM-L_E280_TREFHT.nc')
       file_eoi400 = (FILESTART + 'NorESM-L/NorESM-L_Eoi400_TREFHT.nc')
       fielduse = 'Reference height temperature'
       lsm_e280 = (FILESTART + 'NorESM-L/NorESM-L_E280_land_sea_mask.nc')
       lsm_eoi400 = (FILESTART + 'NorESM-L/NorESM-L_Eoi400_land_sea_mask.nc')
       fieldlsm = 'Fraction of sfc area covered by land'


    if MODELNAME  == 'MRI-CGCM2.3':
        file_e280 = (FILESTART + 'MRI-CGCM2.3/tas/e280.tas.')
        file_eoi400 = (FILESTART + 'MRI-CGCM2.3/tas/eoi400.tas.')
        fielduse = 'near surface air temperature [degC]'
        lsm_e280 = (FILESTART + 'MRI-CGCM2.3/sftlf.nc')
        lsm_eoi400 = lsm_e280
        fieldlsm = 'landsea mask [0 - 1]'


    if MODELNAME  == 'GISS2.1G':
        start = '/nfs/hera1/earjcti/PLIOMIP2/GISS2.1G/'
        mid = 'e280/tas_Amon_GISS-E2-1-G_piControl_r1i1p1f1_gn_'
        file_e280 = ([start + mid + '490101-495012.nc', 
                      start + mid + '495101-500012.nc'])
        mid = 'eoi400/tas_Amon_GISS-E2-1-G_midPliocene-eoi400_r1i1p1f1_gn_'
        file_eoi400 = ([start + mid + '305101-310012.nc',
                        start + mid + '310101-315012.nc'])
        fielduse = 'air_temperature'
        lsm_e280 = start + 'e280/NASA-GISS_PIctrl_all_fland.nc'
        lsm_eoi400 = start + 'eoi400/NASA-GISS_PlioMIP2_all_fland.nc'
        fieldlsm = 'fland'

    if MODELNAME == 'NorESM1-F':
        file_e280 = FILESTART + 'NorESM1-F/NorESM1-F_E280_TREFHT.nc'
        file_eoi400 = FILESTART + 'NorESM1-F/NorESM1-F_Eoi400_TREFHT.nc'
        lsm_e280 = FILESTART + 'NorESM1-F/NorESM1-F_E280_land_sea_mask.nc'
        lsm_eoi400 = FILESTART + 'NorESM1-F/NorESM1-F_Eoi400_land_sea_mask.nc'
        fielduse = 'Reference height temperature'
        fieldlsm =  'Fraction of sfc area covered by land'

        
    if MODELNAME == 'IPSLCM6A':
        start = '/nfs/hera1/earjcti/PLIOMIP2/IPSLCM6A/'
        file_e280 = start + 'tas_Amon_IPSL-CM6A-LR_piControl_r1i1p1f1_gr_285001-304912.nc'
        file_eoi400 = (start + 'tas_Amon_IPSL-CM6A-LR_midPliocene-eoi400_' + 
                       'r1i1p1f1_gr_185001-204912.nc')
        lsm_e280 = start + 'sftlf_fx_IPSL-CM6A-LR_piControl_r1i1p1f1_gr.nc'
        lsm_eoi400 = start + 'sftlf_fx_IPSL-CM6A-LR_midPliocene-eoi400_r1i1p1f1_gr.nc'
        fielduse = 'air_temperature'
        fieldlsm = 'land_area_fraction'

    if MODELNAME == 'IPSLCM5A':
        start = '/nfs/hera1/earjcti/PLIOMIP2/IPSLCM5A/'
        file_e280 = FILESTART + 'IPSLCM5A/PI.NearSurfaceTemp_tas_3600_3699_monthly_TS.nc'
        file_eoi400 = (FILESTART + 'IPSLCM5A/Eoi400.NearSurfaceTemp_tas_3581_3680_monthly_TS.nc')
        lsm_e280 = start + 'E280_LSM_IPSLCM5A.nc'
        lsm_eoi400 = start + 'Eoi400_LSM_IPSLCM5A.nc'
        fielduse = 'Tas'
        fieldlsm = ['Fraction ter', 'Fraction lic']
   

    if MODELNAME == 'IPSLCM5A2':
        start = '/nfs/hera1/earjcti/PLIOMIP2/IPSLCM5A/'
        file_e280 = FILESTART + 'IPSLCM5A2/PI.NearSurfaceTemp_tas_6110_6209_monthly_TS.nc'
        file_eoi400 = (FILESTART + 'IPSLCM5A2/Eoi400.NearSurfaceTemp_tas_3381_3480_monthly_TS.nc')
        lsm_e280 = start + 'E280_LSM_IPSLCM5A.nc'
        lsm_eoi400 = start + 'Eoi400_LSM_IPSLCM5A.nc'
        fielduse = 'Temperature 2m'
        fieldlsm = ['Fraction ter', 'Fraction lic']

    if MODELNAME == 'HadGEM3':
        start = '/nfs/hera1/pliomip2/data/HadGEM3_new/'
        file_e280 = start + 'climatologies/E280/clims_hadgem3_pi_airtemp_final.nc'
        file_eoi400 = start + 'climatologies/Eoi400/clims_hadgem3_pliocene_airtemp_final.nc'
        fielduse = 'temp'
        lsm_e280 = start + 'hadgem3.mask.nc'
        lsm_eoi400 = lsm_e280
        fieldlsm = 'land_binary_mask'
            
            
      
    retdata = [fielduse, file_e280, file_eoi400,
               fieldlsm, lsm_e280, lsm_eoi400]
    return(retdata)


##########################################################
# main program

FILENAME  =  ' '
LINUX_WIN  =  'l'
MODELNAME  = "HadGEM3" # MIROC4m  COSMOS CCSM4_UoT EC-Earth3.1
                   # HadCM3 MRI-CGCM2.3
                   # IPSLCM5A,  IPSLCM5A2
                   # NorESM1-F NorESM-L
                   # IPSLCM6A GISS2.1G
                   # CCSM4-2deg, CESM1.2
                   # CCSM4
                   # EC-Earth3.3 CESM2 (b.e21)
                  
FIELDNAMEIN = ['tas']

if LINUX_WIN  == 'l':
    FILESTART = '/nfs/hera1/pliomip2/data/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'


# call program to get model dependent names
# fielduse,  and  filename
retdata = getnames()

FIELDNAME = retdata[0]
FILENAME_PI = retdata[1]
FILENAME_PLIO = retdata[2]
FIELDLSM = retdata[3]
LSM_PI = retdata[4]
LSM_PLIO = retdata[5]

print('fieldname is',FIELDNAME)

get_land_sea_contrast()

#sys.exit(0)
