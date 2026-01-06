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
    # sitename, sitelat, sitelon, min WMMT CA, max WMMT CA,
    # min CMMT CA, max CMMT CA

    # LP is late pliocene, EP is early pliocene
    # these are late pliocene
    # these are from popova et al 2012 using coexistence approach.
    sitedata.append(['Mirny (LP)', 55, 82, 18.8, 24.6, -0.3, 0.7])
    sitedata.append(['Merkutlinskiy (LP)', 56, 72, 17.3, 23.8, -3.8, 6.2])
    sitedata.append(['Kabinet (LP)', 55, 80, 21.6, 24.4, -4.4, 4.6])
    sitedata.append(['Delyankir (LP)', 63, 133, 18.9, 24.9, -6.9, 1.3])
    sitedata.append(['Chernoluche (LP)', 55, 73, 19.6, 20.3, -5.9, 0.7])
    sitedata.append(['Blizkiy (LP)', 64, 162, 15.6, 23.3, -12.8, 5.2])
    sitedata.append(['42km (LP)', 55, 80, 21.6, 23.3, -4.4, 0.7])
    sitedata.append(['Tnekveem (EP)', 66, 177, 18.9, 25.6, -11.8, 5.8])
    sitedata.append(['Hydzhak (EP)', 63, 147, 18.8, 24.9, -8.7, 1.3])

    # lake baikal is from demske 2002.
    sitedata.append(['Lake Baikal (3.57-3.15Ma)', 56, 108, 13.0, 24.0, -15, 5])
    
   

   
    sites = []
    lats = []
    lons = []
    WMMT_veg_min = []
    WMMT_veg_max = []
    CMMT_veg_min = []
    CMMT_veg_max = []
    
    for site in sitedata:
        sites.append(site[0])
        lats.append(site[1])
        lons.append(site[2])
        WMMT_veg_min.append(site[3])
        WMMT_veg_max.append(site[4])
        CMMT_veg_min.append(site[5])
        CMMT_veg_max.append(site[6])
   
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
             np.asarray(WMMT_veg_max), np.asarray(CMMT_veg_min),
             np.asarray(CMMT_veg_max))


def plot_figure(WMMT_veg_min, WMMT_veg_max, CMMT_veg_min, CMMT_veg_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                labels):

    """
    this subroutine tries to plot the figure for the paper which shows a nice
    DMC 
    """


    fig1 = plt.figure(figsize=[10.0, 10])
    ax1 = plt.axes(frameon=False)
    ax1.get_xaxis().tick_top()
    ax1.axes.get_yaxis().set_visible(False)

    nsites=len(WMMT_veg_min)
    yarray = np.arange(1, nsites+1, 1)
 
    # plot warm and cold month temperature
    ax1.hlines(y=yarray, xmin=WMMT_veg_min, xmax= WMMT_veg_max,color='red')
    ax1.hlines(y=yarray, xmin=CMMT_veg_min, xmax= CMMT_veg_max,color='blue')
    plt.scatter((WMMT_veg_max + WMMT_veg_min) / 2.0, 
                yarray, color='red', marker='^',
                label='CA WMMT')
    plt.scatter((CMMT_veg_max + CMMT_veg_min) / 2.0, 
                 yarray, color='blue', marker='^',
                label =  'CA CMMT')
  
       
    # try plotting axis
  
    for j in range(0, nsites):
        plt.text(-35.0, yarray[j], labels[j], ha='right')

     

    
    plt.scatter(mmm_WMT, yarray + 0.2, color='black', s=50)
    plt.scatter(mmm_WMT, yarray + 0.2, color='red', 
                s=25, label='MMM WMMT')
    plt.scatter(mmm_CMT, yarray + 0.2, color='black', s=50)
    plt.scatter(mmm_CMT, yarray + 0.2, color='blue', 
                s=25, label='MMM CMMT')
  
    plt.ylim(12, -0.5)

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
        
    
    
    plt.legend(loc='lower center', ncol=3)

    fig1.suptitle('a) Pliocene DMC', x=0.1, ha='left', fontsize=16)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/seasonal_dmc_alternative_sites.eps')
    plt.savefig(fileout)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/seasonal_dmc_alternative_sites.png')
    plt.savefig(fileout)
    plt.close()



   

  
def main():
    """
    calling structure
    a) get's model data
    b) get's proxy data
    c) plots model data with proxy data on top
    d) plots change in seasonal cycle at MI, BP and LB
    """

   
    # get land observations and cru temperature at land points 
    (sites, land_lats, land_lons, WMMT_veg_min, WMMT_veg_max, 
     CMMT_veg_min, CMMT_veg_max) =  get_land_obs()
    
    
    
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



    # plot data
  
    plot_figure(WMMT_veg_min, WMMT_veg_max, 
                CMMT_veg_min, CMMT_veg_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                sites)
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
             ,  'MRI2.3'
              ]

#MODELNAMES = ['CESM2']

NSAT_MMM_FILE = (FILESTART + 
                 'regridded100/NearSurfaceTemperature_multimodelmean.nc')

LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')

main()

#sys.exit(0)
