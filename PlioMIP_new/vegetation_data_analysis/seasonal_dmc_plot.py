#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created January 2021 by Julia

This program will produce a lat /lon dmc plot from Ulrichs spreadsheet

"""

import numpy as np
import pandas as pd
import iris
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.quickplot as qplt
import iris.plot as iplt
import cartopy.crs as ccrs
import netCDF4

import sys



def get_MMM_data(latreq, lonreq):
    """
    read in MMM data from the pliocene and the preindustrial 
    return the temperature at the list of sites
    """

    plio_cube = iris.load_cube(NSAT_MMM_FILE,
                              'NearSurfaceTemperaturemean_mPWP')
    pi_cube = iris.load_cube(NSAT_MMM_FILE,
                              'NearSurfaceTemperaturemean_pi')
   
    nsites = len(latreq)
    plio_mmm_array = np.zeros(nsites)
    pi_mmm_array = np.zeros(nsites)

    for i in range(0,nsites):
        # modellon is whole numbers from 0-360
        # lat is half numbers from -89.5 to 89.5

        modlon = np.around(lonreq[i])
        if modlon < 0: modlon = modlon + 360.


        lat_ix = ((np.abs(plio_cube.coord('latitude').points 
                         - latreq[i])).argmin())
        lon_ix = ((np.abs(plio_cube.coord('longitude').points 
                         - modlon)).argmin())
    
        plio_mmm_array[i] = plio_cube.data[lat_ix, lon_ix]
        pi_mmm_array[i] = pi_cube.data[lat_ix, lon_ix]
   
    return plio_mmm_array, pi_mmm_array

def get_lsm_names(model, period):
    """
    gets the names for each of the land sea masks
    period 0 =e280, period 1 = eoi400
    """
    if model == 'CESM2':
        lsm = [DATABASE + 'NCAR/b.e12.B1850.f09_g16.preind.cam.h0.LANDFRAC.0701.0800.nc', DATABASE + 'NCAR/b.e21.B1850.f09_g17.PMIP4-midPliocene-eoi400.001.cam.h0.LANDFRAC.1101.1200.nc']
        fieldlsm = 'Fraction of sfc area covered by land'

    if model == 'COSMOS':
        lsm = ["/nfs/hera1/pliomip2/data/AWI/COSMOS/land_sea_masks/E280_et_al/E280.slf.atm.nc", "/nfs/hera1/pliomip2/data/AWI/COSMOS/land_sea_masks/Eoi400_et_al/Eoi400.slf.atm.nc"]
        fieldlsm = "SLF"

    if model == 'EC-Earth3.3':
        lsm =  [DATABASE + 'EC-Earth3.3/EC-Earth3.3_PI_LSM.nc',
                DATABASE + 'EC-Earth3.3/EC-Earth3.3_mPlio_LSM.nc']
        fieldlsm = 'Land/sea mask'

    if model == 'CESM1.2':
        lsm = [DATABASE + 'NCAR/b.e12.B1850.f09_g16.preind.cam.h0.LANDFRAC.0701.0800.nc', DATABASE + 'NCAR/b.e12.B1850.f09_g16.PMIP4-pliomip2.cam.h0.LANDFRAC.1101.1200.nc']
        fieldlsm = 'Fraction of sfc area covered by land'       
  
    if model   ==  'MIROC4m':
        lsm = [DATABASE + 'MIROC4m/sftlf/MIROC4m_Exxx_fx_sftlf.nc', 
               DATABASE + 'MIROC4m/sftlf/MIROC4m_Eoixxx_fx_sftlf.nc']
        fieldlsm = "sftlf"

    if model  == 'HadCM3':
        lsm = [DATABASE+'LEEDS/HadCM3/e280/qrparm.mask.nc',
               DATABASE+'LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc']
        fieldlsm = 'LAND MASK (LOGICAL: LAND=TRUE)'

    if model == 'CCSM4':
        lsm = [DATABASE + 'NCAR/b40.B1850.f09_g16.preind.cam.h0.LANDFRAC.0081.0180.nc', DATABASE + 'NCAR/b40.B1850.f09_g16.PMIP4-pliomip2.LANDFRAC.1001.1100.nc']
        fieldlsm = 'Fraction of sfc area covered by land'

    if model == 'CCSM4-Utr':
        lsm = [DATABASE + 'Utrecht/CESM1.0.5/E280/land_sea_mask_Amon_CESM1.0.5_b.PI_1pic_f19g16_NESSC_control_r1i1p1f1_gn.nc', DATABASE + 'Utrecht/CESM1.0.5/Eoi400/land_sea_mask_Amon_CESM1.0.5_b.PLIO_5Ma_Eoi400_f19g16_NESSC_control_r1i1p1f1_gn.nc']
        fieldlsm = 'LANDMASK[D=1]'
  
    if model == 'CCSM4-UoT':
        start = DATABASE + 'UofT/UofT-CCSM4/'
        lsm = [start + 'for_julia/E_mask.nc', start + 'for_julia/Eoi_mask.nc']
        fieldlsm = 'gridbox land fraction'
      
    if model == 'NorESM-L':
       lsm = [DATABASE + 'NorESM-L/NorESM-L_E280_land_sea_mask.nc',
              DATABASE + 'NorESM-L/NorESM-L_Eoi400_land_sea_mask.nc']
       fieldlsm = 'Fraction of sfc area covered by land'


    if model  == 'MRI2.3':
        lsm = [DATABASE + 'MRI-CGCM2.3/sftlf.nc', 
               DATABASE + 'MRI-CGCM2.3/sftlf.nc']
        fieldlsm = 'landsea mask [0 - 1]'


    if model  == 'GISS2.1G':
        start = '/nfs/hera1/earjcti/PLIOMIP2/GISS2.1G/'
        lsm = [start + 'e280/NASA-GISS_PIctrl_all_fland.nc',
               start + 'eoi400/NASA-GISS_PlioMIP2_all_fland.nc']
        fieldlsm = 'fland'

    if model == 'NorESM1-F':
        lsm = [DATABASE + 'NorESM1-F/NorESM1-F_E280_land_sea_mask.nc',
               DATABASE + 'NorESM1-F/NorESM1-F_Eoi400_land_sea_mask.nc']
        fieldlsm =  'Fraction of sfc area covered by land'

        
    if model == 'IPSLCM6A':
        start = '/nfs/hera1/earjcti/PLIOMIP2/IPSLCM6A/'
        lsm = [start + 'sftlf_fx_IPSL-CM6A-LR_piControl_r1i1p1f1_gr.nc',
              start + 'sftlf_fx_IPSL-CM6A-LR_midPliocene-eoi400_r1i1p1f1_gr.nc']
        fieldlsm = 'land_area_fraction'

    if model == 'IPSLCM5A':
        start = '/nfs/hera1/earjcti/PLIOMIP2/IPSLCM5A/'
        lsm = [start + 'E280_LSM_IPSLCM5A.nc',start + 'Eoi400_LSM_IPSLCM5A.nc']
        fieldlsm = ['Fraction ter', 'Fraction lic']
   

    if model == 'IPSLCM5A2':
        start = '/nfs/hera1/earjcti/PLIOMIP2/IPSLCM5A/'
        lsm = [start + 'E280_LSM_IPSLCM5A.nc',
               start + 'Eoi400_LSM_IPSLCM5A.nc']
        fieldlsm = ['Fraction ter', 'Fraction lic']

    if model == 'HadGEM3':
        start = '/nfs/hera1/pliomip2/data/HadGEM3_new/'
        lsm = [start + 'hadgem3.mask.nc', start + 'hadgem3.mask.nc']
        fieldlsm = 'land_binary_mask'
            
            
    return lsm[period], fieldlsm


def get_cru_temp(lats, lons):
    """
    get's the cru temperature at the given latitude and longitude
    """
    
    crufile = ('/nfs/hera1/earjcti/regridded/CRUTEMP/' + 
               'E280.NearSurfaceTemperature.mean_month.nc')
    cube = iris.load_cube(crufile)
    
    
    cru_min_temp = np.zeros(len(lats))
    cru_max_temp = np.zeros(len(lats))
    for i, lat in enumerate(lats):
        if lons[i] > 180.:
            lon = lons[i]-360.
        else:
            lon = lons[i]
        lat_ix = (np.abs(cube.coord('latitude').points - lat)).argmin()
        lon_ix = (np.abs(cube.coord('longitude').points - lon)).argmin()
        
   
        cru_temp = cube.data[:, lat_ix, lon_ix]
        if np.isfinite(cru_temp[0]):
            pass
        else:
            # get an average of surrounding ones
            surround = [cube.data[:, lat_ix + 1, lon_ix],
                        cube.data[:, lat_ix - 1, lon_ix],
                        cube.data[:, lat_ix, lon_ix + 1],
                        cube.data[:, lat_ix, lon_ix -1],
                        ]
            cru_temp = np.ma.mean(surround, axis=0)
            cru_max_temp[i] = np.max(cru_temp)
        cru_min_temp[i] = np.min(cru_temp)
      
     
    return cru_min_temp, cru_max_temp

###########################
def get_land_sea_mask(model, period):
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


    lsm, fieldlsm = get_lsm_names(model,period)

    ############################################
    if model == 'IPSLCM5A' or model == 'IPSLCM5A2':
        lsm_cube = get_ipsl_lsm(lsm, fieldlsm)
    elif model == 'HadGEM3':
        f = netCDF4.Dataset(lsm, "r")
        print(f.variables['longitude'])
        longitude = iris.coords.DimCoord(f.variables['longitude'], 
                             standard_name = 'longitude', units='degrees')
        latitude = iris.coords.DimCoord(f.variables['latitude'], 
                             standard_name = 'latitude', units='degrees')
        lsm_cube = iris.cube.Cube(np.squeeze(f.variables['lsm'][:]),
                             long_name='lsm', var_name='lsm', units=None, 
                             attributes=None, cell_methods=None, 
                             dim_coords_and_dims=[(latitude,0), (longitude,1)])
     
    else:
        lsm_cube = iris.util.squeeze(iris.load_cube(lsm, fieldlsm))
     
    lsm_cube2 = change_to_2d(lsm_cube)
   
   
    if model == 'IPSLCM6A':
        lsm_cube2.data = lsm_cube2.data / 100.0
       

    # regrid
    cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
    lsm_cube3 = lsm_cube2.regrid(cubegrid, iris.analysis.Linear())
   

    lsm_cube3.var_name = 'land_mask'
    lsm_cube3.long_name = 'land_mask'
    
    return lsm_cube3

def check_lsm(lsm_lons, lsm_lats, lsm_data, latrq, lonrq):
    """
    if our model is a sea point then set index to nan
    """

    lat_ix = (np.abs(lsm_lats - latrq)).argmin()
    lon_ix = (np.abs(lsm_lons - lonrq)).argmin()
       
    if lsm_data[lat_ix, lon_ix] <  0.5:
        # check to south, north, east, west
#        print(lsm_data[lat_ix - 1, lon_ix],lsm_data[lat_ix + 1, lon_ix],lsm_data[lat_ix, lon_ix-1],lsm_data[lat_ix, lon_ix+1], lsm_data[lat_ix - 1, lon_ix - 1],lsm_data[lat_ix + 1, lon_ix+1 ],lsm_data[lat_ix + 1, lon_ix-1],lsm_data[lat_ix-1, lon_ix+1])
#        if lsm_data[lat_ix - 1, lon_ix] >  0.5:
#            lat_ix = lat_ix -1
#        elif lsm_data[lat_ix + 1, lon_ix] >  0.5:
#            lat_ix = lat_ix + 1
#        elif lsm_data[lat_ix, lon_ix - 1] > 0.5:
#            lon_ix = lon_ix - 1
#        elif lsm_data[lat_ix, lon_ix + 1] > 0.5:
#            lon_ix = lon_ix + 1
#        elif lsm_data[lat_ix - 1, lon_ix + 1] > 0.5:
#            lat_ix = lat_ix -1
#            lon_ix = lon_ix + 1
#        elif lsm_data[lat_ix - 1, lon_ix - 1] > 0.5:
#            lat_ix = lat_ix -1
#            lon_ix = lon_ix - 1
#        elif lsm_data[lat_ix + 1, lon_ix + 1] > 0.5:
#            lat_ix = lat_ix +1
#            lon_ix = lon_ix + 1
#        elif lsm_data[lat_ix + 1, lon_ix - 1] > 0.5:
#            lat_ix = lat_ix +1
#            lon_ix = lon_ix - 1
       # elif lsm_data[lat_ix, lon_ix - 2] > 0.5:
       #     lon_ix = lon_ix - 2
       # elif lsm_data[lat_ix, lon_ix + 2] > 0.5:
       #     lon_ix = lon_ix + 2
       # elif lsm_data[lat_ix, lon_ix - 3] > 0.5:
       #     lon_ix = lon_ix - 3
       # elif lsm_data[lat_ix, lon_ix + 3] > 0.5:
       #     lon_ix = lon_ix + 3
       # elif lsm_data[lat_ix - 2, lon_ix] > 0.5:
       #     lat_ix = lat_ix -2
       # elif lsm_data[lat_ix + 2, lon_ix] > 0.5:
       #     lat_ix = lat_ix + 2
      
       
  #      else:
            lat_ix = np.nan
            lon_ix = np.nan
            
       # print('new',lsm_data[lat_ix, lon_ix], lat_ix, lon_ix, latrq, lonrq, lsm_lons[lon_ix], lsm_lats[lat_ix])
       # sys.exit(0)

            
        

    return lat_ix, lon_ix

def get_single_model(model, latreq, lonreq, period):
    """
    read in the pliocene data from 'model'  return the temperatures
    at the list of sites
    """
    # get lsm
    if period == 'E280':
        lsm_cube  = get_land_sea_mask(model, 0)
    if period == 'EOI400':
        lsm_cube  = get_land_sea_mask(model, 1)

    filename = ('/nfs/hera1/earjcti/regridded100/' + model +
                '/' + period + '.NearSurfaceTemperature.mean_month.nc')
  
    print(filename)
    plio_cube = iris.load_cube(filename)
   
    nsites = len(latreq)
    plio_minval_array = np.zeros(nsites)
    plio_maxval_array = np.zeros(nsites)

    plio_cube_lats = plio_cube.coord('latitude').points
    plio_cube_lons = plio_cube.coord('longitude').points

    lsm_cube_lats = lsm_cube.coord('latitude').points
    lsm_cube_lons = lsm_cube.coord('longitude').points
    lsm_cube_data = lsm_cube.data
   
  
    if plio_cube_lats.any() != lsm_cube_lats.any():
        print('data cube does not match lsm lat')
    if plio_cube_lats.any() != lsm_cube_lons.any():
        print('data cube does not match lsm lat')

    for i in range(0,nsites):
        # modellon is whole numbers from 0-360
        # lat is half numbers from -89.5 to 89.5

        modlon = np.around(lonreq[i])
        if modlon < 0: modlon = modlon + 360.

        (lat_ix, 
         lon_ix) = check_lsm(lsm_cube_lons, lsm_cube_lats, 
                                    lsm_cube_data, latreq[i], modlon)

        if np.isfinite(lat_ix):
            plio_array = plio_cube.data[:, lat_ix, lon_ix]
        else:
            plio_array = np.zeros(12)
            plio_array[:] = np.nan
        plio_minval_array[i] = np.min(plio_array)
        plio_maxval_array[i] = np.max(plio_array)
   
    return plio_minval_array, plio_maxval_array

 
def get_land_obs():
    """
    these have been obtained from various sources so I am just typing them in
    """
    sitedata = []
    # site data is
    # sitename, sitelat, sitelon, min WMMT veg, max WMMT veg,
    # min WMMT beetle, max WMMT beetle, min CMMT veg, max CMMT veg
    # min CMMT beetle, max CMMT beetle, modern obs WMMT, modern obs CMMT
    #sitedata.append(['Lake El\'gygytgyn', 67, 172, 15.0, 16.0,
    #                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
    #                 8.0, np.nan])
    sitedata.append(['Lake Baikal', 56, 108, 15.28, 17.52,
                     np.nan, np.nan, -1.67, 1.07, np.nan, np.nan, 
                     15.3, -17.4])
    sitedata.append(['Lake El\'gygytgyn', 67, 172, -36.8, -30.4,
                    np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, 
                     8.0, np.nan])
    sitedata.append(['Near Meighen Island', 77.5, 261, 19.6, 20.5,
                     11.5, 13.5, -11.6, -11.4,
                     -33.0, -18.5, 4.1, -42.5])
    sitedata.append(['Beaver Pond', 79, 278, 18.4, 20.9,
                     np.nan, np.nan, -12.2, -11.5,
                     np.nan, np.nan,  7.1, -39.7])
    sitedata.append(['Flyes Leaf Bed', 79, 278, 19.7, 21.1,
                     np.nan, np.nan, -12.8, -9.1,
                     np.nan, np.nan,  7.1, -39.7])
    sitedata.append(['Lost Chicken Mine', 64, 218, 12.0, 12.0, 
                     13.5, 16.0, -2.0, -2.0, -27.75, -19.25, 15.3, -25.1])
   

    sites_relaxed_coex = []
    # sitename WMMTmin WMMT max, CMMT min, CMMTmax, BMA WMMTmin, BMA_WMMTmax,
    # BMA_CMMTmin, BMA_CMMTmax
    sites_relaxed_coex.append(['Lake Baikal', np.nan, np.nan, np.nan, np.nan,
                               np.nan, np.nan, np.nan, np.nan])
    sites_relaxed_coex.append(['Lake El\'gygytgyn', np.nan, np.nan, np.nan, np.nan, 15.0, 16.0, -36.8, -30.4 ])
    sites_relaxed_coex.append(['Near Meighen Island', 18.1, 22.8, -21.7, -7.9,
                               np.nan, np.nan, np.nan, np.nan])
    sites_relaxed_coex.append(['Beaver Pond', 18.1, 22.4, -21.7, -8.1,
                               np.nan, np.nan, np.nan, np.nan])
    sites_relaxed_coex.append(['Flyes Leaf Bed', 18.1, 22.7, -16.9, -6.4,
                               np.nan, np.nan, np.nan, np.nan])
    sites_relaxed_coex.append(['Lost Chicken Mine', np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
             

    # lcm observations from elias and matthews 2002
    # lake baikal is from what ulrich sent me.  MI and BP and FLBfrom Fletcher
    # Lake E from Brigette-greeme mean july temp of +8 and average winter lows of 35degC
    # Meighen Island is 'upper'
    
    sites = []
    lats = []
    lons = []
    WMMT_veg_min = []
    WMMT_veg_max = []
    WMMT_coex_min = []
    WMMT_coex_max = []
    WMMT_BMA_min = []
    WMMT_BMA_max = []
    WMMT_beetle_min = []
    WMMT_beetle_max = []
    WMMT_modern_obs = []
    CMMT_veg_min = []
    CMMT_veg_max = []
    CMMT_coex_min = []
    CMMT_coex_max = []
    CMMT_BMA_min = []
    CMMT_BMA_max = []
    CMMT_beetle_min = []
    CMMT_beetle_max = []
    CMMT_modern_obs = []
    
    for i, info in enumerate(sitedata):
        info2 = sites_relaxed_coex[i]                          
        sites.append(info[0])
        lats.append(info[1])
        lons.append(info[2])
        WMMT_veg_min.append(info[3])
        WMMT_veg_max.append(info[4])
        WMMT_beetle_min.append(info[5])
        WMMT_beetle_max.append(info[6])
        CMMT_veg_min.append(info[7])
        CMMT_veg_max.append(info[8])
        CMMT_beetle_min.append(info[9])
        CMMT_beetle_max.append(info[10])
        WMMT_modern_obs.append(info[11])
        CMMT_modern_obs.append(info[12])
        if info[0] == info2[0]:
           WMMT_coex_min.append(info2[1])
           WMMT_coex_max.append(info2[2])
           CMMT_coex_min.append(info2[3])
           CMMT_coex_max.append(info2[4])
           WMMT_BMA_min.append(info2[5])
           WMMT_BMA_max.append(info2[6])
           CMMT_BMA_min.append(info2[7])
           CMMT_BMA_max.append(info2[8])
        else:
           print,'check names match',info[0],info2[0]
           sys.exit(0)      
   
    labels = []
    deg= u'\N{DEGREE SIGN}'
    for i, site in enumerate(sites):
  #     label = ''.join([c for c in site if c.isupper()])
        latstr = np.str(lats[i]) + deg + 'N'
        if lons[i] >180:
           lonstr = np.str((lons[i] - 360.) * -1.0) + deg +  'W'
        else:
           lonstr = np.str(lons[i]) + deg + 'E'
         
        label = site + '\n (' +  latstr + ',' +  lonstr + ')'
        labels.append(label)
   
    return  (labels, lats, lons, np.asarray(WMMT_veg_min),
             np.asarray(WMMT_veg_max), np.asarray(WMMT_beetle_min),
             np.asarray(WMMT_beetle_max), np.asarray(CMMT_veg_min),
             np.asarray(CMMT_veg_max), np.asarray( CMMT_beetle_min),
             np.asarray(CMMT_beetle_max), np.asarray(WMMT_modern_obs),
             np.asarray(CMMT_modern_obs), np.asarray(WMMT_coex_min),
             np.asarray(WMMT_coex_max),np.asarray(CMMT_coex_min),
             np.asarray(CMMT_coex_max),np.asarray(WMMT_BMA_min),
             np.asarray(WMMT_BMA_max),np.asarray(CMMT_BMA_min),
             np.asarray(CMMT_BMA_max))


def plot_figure(WMMT_veg_min, WMMT_veg_max, WMMT_beetle_min, WMMT_beetle_max,
                CMMT_veg_min, CMMT_veg_max, CMMT_beetle_min, CMMT_beetle_max,
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                labels, pi_ind, cru_min, cru_max,WMMT_coex_min,
                WMMT_coex_max,CMMT_coex_min,CMMT_coex_max,WMMT_BMA_min,
                WMMT_BMA_max,CMMT_BMA_min,CMMT_BMA_max):
    """
    this subroutine tries to plot the figure for the paper which shows a nice
    DMC 
    """


    fig1 = plt.figure(figsize=[10.0, 8.0])
    ax1 = plt.axes(frameon=False)
    ax1.get_xaxis().tick_top()
    ax1.axes.get_yaxis().set_visible(False)

    nsites=len(WMMT_veg_min)
    yarray = np.arange(1, nsites+1, 1)
 
    # plot warm and cold month temperature
    ax1.hlines(y=yarray, xmin=WMMT_veg_min, xmax= WMMT_veg_max,color='red')
    ax1.hlines(y=yarray, xmin=CMMT_veg_min, xmax= CMMT_veg_max,color='blue')
    if pi_ind == 'y': # plot cru temps
        plt.scatter((WMMT_veg_max + WMMT_veg_min) / 2.0, 
                yarray, color='red', marker='^',
                label='modern WMMT at site')
        plt.scatter((CMMT_veg_max + CMMT_veg_min) / 2.0, 
                 yarray, color='blue', marker='^',
                label =  'modern CMMT at site')
  
        plt.scatter(cru_max, yarray, color='red', marker = 'v', label = 'CRU WMMT')
        plt.scatter(cru_min, yarray, color='blue', marker = 'v', label = 'CRU CMMT')   
    else:
        plt.scatter((WMMT_veg_max + WMMT_veg_min) / 2.0, 
                yarray, color='red', marker='^',
                label='paleovegetation WMMT')
        plt.scatter((CMMT_veg_max + CMMT_veg_min) / 2.0, 
                 yarray, color='blue', marker='^',
                label =  'paleovegetation CMMT')
  
        # put an arrow on LCM
        plt.arrow(CMMT_veg_max[5], yarray[5], -10, 0, linestyle='dotted', 
                  color='blue', head_width=0.1, head_length=0.5)
   
    # try plotting axis
  
    for j in range(0, nsites):
        if pi_ind == 'y':
           plt.text(20.0,
                 yarray[j], labels[j], ha='left')
        else:
            plt.text(-35.0,
                 yarray[j], labels[j], ha='right')

     

    
    plt.scatter(mmm_WMT, yarray + 0.2, color='black', s=50)
    plt.scatter(mmm_WMT, yarray + 0.2, color='red', 
                s=25, label='MMM WMMT')
    plt.scatter(mmm_CMT, yarray + 0.2, color='black', s=50)
    plt.scatter(mmm_CMT, yarray + 0.2, color='blue', 
                s=25, label='MMM CMMT')
  
    plt.ylim(8, -0.5)

    # plot individual models for pliocene
    for i in range(0, len(MODELNAMES)):
        if i == 0:
            plt.scatter(all_models_plio_WMT[:, i], yarray+0.2, color='red', 
                        marker = 'x', s=10, label='models WMMT')
            plt.scatter(all_models_plio_CMT[:, i], yarray+0.2, color='blue', 
                        marker = 'x', s=10, label='models CMMT')
        else:
            plt.scatter(all_models_plio_WMT[:, i], yarray+0.2, color='red', 
                        marker = 'x',s=10)
            plt.scatter(all_models_plio_CMT[:, i], yarray+0.2, color='blue', 
                        marker = 'x', s=10)
        
    if pi_ind == 'y':
        plt.xlim(-50, 25)
        plt.hlines(y=-0.45, xmin=-50., xmax=25., linewidth=0.5)
        plt.text(-45.0, 0.0, 'Cold Month Temperature (degC)', ha='left', color='blue')
        plt.text(20.0, 0.0, 'Warm Month Temperature (deg C)', ha='right', color='red') 
    else:
        plt.xlim(-40, 30)
        plt.hlines(y=-0.45, xmin=-40., xmax=30., linewidth=0.5)
        plt.text(-35.0, 0.0, 'Cold Month Temperature (degC)', ha='left', color='blue')
        plt.text(25.0, 0.0, 'Warm Month Temperature (deg C)', ha='right', color='red')
        # plot all beetle assemblage results
        plt.hlines(yarray, xmin=WMMT_beetle_min, 
               xmax=WMMT_beetle_max, color='green',linestyle='dashed',
               label='beetle assemblage data')
        plt.hlines(yarray, xmin=CMMT_beetle_min, 
               xmax=CMMT_beetle_max, color='green', linestyle='dashed')
        # plot coexistence results
        plt.hlines(yarray+0.1, xmin=WMMT_coex_min, 
               xmax=WMMT_coex_max, color='red',linestyle='dotted',
               label='relaxed coexistence approach data')
        plt.hlines(yarray+0.1, xmin=CMMT_coex_min, 
               xmax=CMMT_coex_max, color='blue', linestyle='dotted')
        # plot BMA results
        plt.scatter((WMMT_BMA_max + WMMT_BMA_min) / 2.0, 
                yarray, color='red', marker='v',
                label='WMMT BMA')
        plt.scatter((CMMT_BMA_max + CMMT_BMA_min) / 2.0, 
                yarray, color='blue', marker='v',
                label='WMMT BMA')
        plt.hlines(yarray, xmin=WMMT_BMA_min, 
               xmax=WMMT_BMA_max, color='red')
        plt.hlines(yarray, xmin=CMMT_BMA_min, 
               xmax=CMMT_BMA_max, color='blue')
    
    
    plt.legend(loc='lower center', ncol=3)

    if pi_ind == 'y':
        fig1.suptitle('b) modern/preindustrial DMC', x=0.1, ha='left', fontsize=16)
        fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
                   'vegetation/seasonal_dmc_plot_pi.eps')
        plt.savefig(fileout)
        fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
                   'vegetation/seasonal_dmc_plot_pi.png')
        plt.savefig(fileout)
    else:
        fig1.suptitle('a) Pliocene DMC', x=0.1, ha='left', fontsize=16)
        fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/seasonal_dmc_plot.eps')
        plt.savefig(fileout)
        fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
                   'vegetation/seasonal_dmc_plot.png')
        plt.savefig(fileout)
    plt.close()



####################################################
def plot_chg_seas_cyc(WMMT_veg,CMMT_veg, all_models_plio_WMT, 
                      all_models_plio_CMT, all_models_pi_WMT,
                      all_models_pi_CMT,
                      sites, mmm_WMT, mmm_CMT, mmm_WMT_pi, mmm_CMT_pi,cru):
    """
    plots the change in the seasonal cycle between pi and mpwp
    """
    colors = {'CESM2' : 'green',
              'HadGEM3': 'green',
              'IPSLCM6A' : 'green', 
              'COSMOS': 'red', 
              'MIROC4m': 'red', 
              'HadCM3': 'green',
              'GISS2.1G' : 'black', 
              'CCSM4' : 'black', 
              'CCSM4-Utr' : 'black', 
              'CCSM4-UoT' : 'purple', 
              'NorESM-L' : 'purple', 
              'NorESM1-F' : 'purple',
              'MRI2.3' : 'red'}

    marker = {'CESM2' : '^',
              'HadGEM3': 'x',
              'IPSLCM6A' : 'o', 
              'COSMOS': '^', 
              'MIROC4m': 'x', 
              'HadCM3': 'o',
              'GISS2.1G' : '^', 
              'CCSM4' : 'x', 
              'CCSM4-Utr' : 'o', 
              'CCSM4-UoT' : '^', 
              'NorESM-L' : 'x', 
              'NorESM1-F' : 'o',
              'MRI2.3' : 'v',
              'EC-Earth3.3' : '^',
              'CESM1.2' : 'o'}

    print('julia')
    fig1 = plt.figure(figsize=[10.0, 8.0])
    ax1 = plt.subplot(211)
   # ax1.get_xaxis().tick_top()
    ax1.axes.get_xaxis().set_visible(False)

    xarray = np.arange(0,len(sites))
    plt.scatter(xarray, np.asarray(WMMT_veg) - np.asarray(CMMT_veg), 
                color='red', marker = 'x', s=20, label='veg data')
    plt.scatter(xarray, np.asarray(cru),color='black',marker = 's',
                label='cru')
    

    for i, model in enumerate(MODELNAMES):
        plt.scatter(xarray+(0.02 * i), np.asarray(all_models_plio_WMT[:,i])
                    -np.asarray(all_models_plio_CMT[:,i]), 
                    color=colors.get(model,'blue'), 
                    marker = marker.get(model,'x'),label=model)
   
    plt.xlim(-1, 8)
    plt.hlines(y=0, xmin=-1., xmax=6., linewidth=0.5)
    plt.vlines(x=-1, ymin=-25, ymax=7, linewidth=0.5)

    plt.legend(loc='right')

    for i in range(0, len(sites)):
        plt.text(xarray[i], 0, sites[i], ha='right',rotation=90)


    # plot change in seasonal cycle
    seas_chg = (np.asarray(all_models_plio_WMT) 
                - np.asarray(all_models_plio_CMT) 
                - np.asarray(all_models_pi_WMT)
                + np.asarray(all_models_pi_CMT))
    ax2 = plt.subplot(212)
    for i, model in enumerate(MODELNAMES):
        plt.scatter(xarray+(0.02 * i),seas_chg[:,i],
                    color=colors.get(model,'blue'), 
                    marker = marker.get(model,'x'),label=model)
        print(model, all_models_plio_WMT[3,i],all_models_plio_CMT[3,i],
              all_models_pi_WMT[3,i],all_models_pi_CMT[3,i])
    sys.exit(0)
    plt.hlines(y=0, xmin=-1., xmax=6., linewidth=0.5)
   

    plt.show()
            
   

  
def main():
    """
    calling structure
    a) get's model data
    b) get's proxy data
    c) plots model data with proxy data on top
    d) plots change in seasonal cycle at MI, BP and LB
    """

   
    # get land observations and cru temperature at land points 
    (sites, land_lats, land_lons, WMMT_veg_min, WMMT_veg_max, WMMT_beetle_min,
     WMMT_beetle_max, CMMT_veg_min, CMMT_veg_max, CMMT_beetle_min,
     CMMT_beetle_max, WMMT_modern_obs, CMMT_modern_obs,WMMT_coex_min,
     WMMT_coex_max,CMMT_coex_min,CMMT_coex_max,WMMT_BMA_min,
     WMMT_BMA_max,CMMT_BMA_min,CMMT_BMA_max) =  get_land_obs()
    
    cru_min_temp, cru_max_temp = get_cru_temp(land_lats, land_lons)

    
    
    # get ind models data
    print(land_lons)
    all_models_plio_WMT = np.zeros((len(land_lons), len(MODELNAMES)))
    all_models_plio_CMT = np.zeros((len(land_lons), len(MODELNAMES)))
    all_models_pi_WMT = np.zeros((len(land_lons), len(MODELNAMES)))
    all_models_pi_CMT = np.zeros((len(land_lons), len(MODELNAMES)))
    all_models_seas_cyc_anom = np.zeros((len(land_lons),len(MODELNAMES)))
    all_models_CMMT_anom = np.zeros((len(land_lons),len(MODELNAMES)))
    all_models_WMMT_anom = np.zeros((len(land_lons),len(MODELNAMES)))
    for i, model in enumerate(MODELNAMES):
        (ind_CMT, ind_WMT) = get_single_model(model, land_lats, 
                                              land_lons, 'EOI400')
        all_models_plio_WMT[:, i] = ind_WMT
        all_models_plio_CMT[:, i] = ind_CMT

        (ind_CMT, ind_WMT) = get_single_model(model, land_lats, 
                                              land_lons, 'E280')
        all_models_pi_WMT[:, i] = ind_WMT
        all_models_pi_CMT[:, i] = ind_CMT

        all_models_seas_cyc_anom[:, i] = (
            (all_models_plio_WMT[:, i] - all_models_plio_CMT[:, i]) -
            (all_models_pi_WMT[:, i] - all_models_pi_CMT[:, i]))
        
        all_models_CMMT_anom[:, i] = (all_models_plio_CMT[:, i] -
                                      all_models_pi_CMT[:, i])

        all_models_WMMT_anom[:, i] = (all_models_plio_WMT[:, i] -
                                      all_models_pi_WMT[:, i])



    mmm_WMT = np.nanmean(all_models_plio_WMT, axis=1)
    mmm_CMT = np.nanmean(all_models_plio_CMT, axis=1)
    mmm_WMT_pi = np.nanmean(all_models_pi_WMT, axis=1)
    mmm_CMT_pi = np.nanmean(all_models_pi_CMT, axis=1)


    print('MMM WMMT avg all sites',np.nanmean(mmm_WMT), 'all',mmm_WMT)
    print('veg data WMMT avg all sites',
          np.nanmean(WMMT_veg_min + WMMT_veg_max) / 2.0, 'all',
          (WMMT_veg_min + WMMT_veg_max) / 2.0)
    print('MMM CMMT avg all sites',np.nanmean(mmm_CMT), 'all',mmm_CMT)
    print('veg data CMMT avg all sites',
          np.nanmean(CMMT_veg_min + CMMT_veg_max) / 2.0, 'all',
          (CMMT_veg_min + CMMT_veg_max) / 2.0)

    # plot data
    dummy = np.zeros(len(WMMT_modern_obs))
    dummy[:] = np.nan
 
    plot_figure(WMMT_veg_min, WMMT_veg_max, WMMT_beetle_min, WMMT_beetle_max,
                CMMT_veg_min, CMMT_veg_max, CMMT_beetle_min, CMMT_beetle_max,
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                sites,'n', dummy, dummy,WMMT_coex_min,
                WMMT_coex_max,CMMT_coex_min,CMMT_coex_max,WMMT_BMA_min,
                WMMT_BMA_max,CMMT_BMA_min,CMMT_BMA_max)
    # plot a dmc for the pi
    plot_figure (WMMT_modern_obs, WMMT_modern_obs, dummy, dummy,
                 CMMT_modern_obs, CMMT_modern_obs, dummy, dummy,
                 all_models_pi_WMT, all_models_pi_CMT, mmm_WMT_pi, mmm_CMT_pi,
                 sites,'y', cru_min_temp, cru_max_temp,dummy, dummy,
                 dummy, dummy, dummy, dummy, dummy, dummy)

    # plot change in seasonal cycle
    plot_chg_seas_cyc((WMMT_veg_min + WMMT_veg_max)/2.0, 
                (CMMT_veg_min +  CMMT_veg_max) / 2.0, 
                all_models_plio_WMT, all_models_plio_CMT,
                all_models_pi_WMT, all_models_pi_CMT,
                      sites, mmm_WMT, mmm_CMT, mmm_WMT_pi, mmm_CMT_pi,
                (cru_max_temp - cru_min_temp) / 2.0   )

##########################################################
# main program

LINUX_WIN = 'l'
FILESTART = '/nfs/hera1/earjcti/'
DATABASE = '/nfs/hera1/pliomip2/data/'

MODELNAMES = [
               'HadGEM3', 'CESM2',
              'IPSLCM6A', 
              'COSMOS', 
              'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
              'MIROC4m', 'IPSLCM5A2', 'HadCM3',
              'GISS2.1G', 'CCSM4', 
              'CCSM4-Utr', 'CCSM4-UoT', 
              'NorESM-L',  'NorESM1-F'
           #  ,  'MRI2.3'
              ]

#MODELNAMES = ['CESM2']

NSAT_MMM_FILE = (FILESTART + 
                 'regridded100/NearSurfaceTemperature_multimodelmean.nc')

LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')

main()

#sys.exit(0)
