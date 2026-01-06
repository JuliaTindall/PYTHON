
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Updated from ~earjcti/PYTHON/PROGRAMS/PlioMIP_new/vegetation_data_analysis/seasonal_dmc_plot_v2.py

Changed April 2023 to just show HadCM3 and HadCM3 with the new parameters



FROM ~earjcti/PYTHON/PROGRAMS/PlioMIP_new/vegetation_data_analysis/seasonal_dmc_plot_v2.py
Created January 2021 by Julia

This program will produce a lat /lon dmc plot from Ulrichs spreadsheet
The difference between this and version 1 is that we will group the data
by dating / proxy

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

       

    
    return lsm_cube3

def check_lsm(lsm_lons, lsm_lats, lsm_data, latrq, lonrq):
    """
    if our model is a sea point then set index to nan
    """

    lat_ix = (np.abs(lsm_lats - latrq)).argmin()
    lon_ix = (np.abs(lsm_lons - lonrq)).argmin()

    if lsm_data[lat_ix, lon_ix] <  0.5:
            lat_ix = np.nan
            lon_ix = np.nan

    return lat_ix, lon_ix

def get_single_model(model, latreq, lonreq, period):
    """
    read in the pliocene data from 'model'  return the temperatures
    at the list of sites
    """
    # get lsm
    if period == 'E280':
        lsmfile = DATABASE+'LEEDS/HadCM3/e280/qrparm.mask.nc'

    if period == 'EOI400':
        lsmfile = DATABASE+'LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc'


    lsm_cube2  = iris.load_cube(lsmfile,'LAND MASK (LOGICAL: LAND=TRUE)')
    # regrid
    cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
    lsm_cube = lsm_cube2.regrid(cubegrid, iris.analysis.Linear())
    lsm_cube = iris.util.squeeze(lsm_cube)
   
    lsm_cube.var_name = 'land_mask'
    lsm_cube.long_name = 'land_mask'


    # get temperature
    filename = ('/nfs/hera1/earjcti/regridded/' + model +
                '/' + period + '.NearSurfaceTemperature.mean_month.nc')
  
    print(filename)
    plio_cube = iris.load_cube(filename)
    if model == 'HadCM3_new':
        plio_cube = plio_cube - 273.15
   
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

def land_reformat(sitedata):
    """
    reformats the land data into different arrays
    """
    
    sites = []
    lats = []
    lons = []
    WMMT_data_min = []
    WMMT_data_max = []
    WMMT_modern_obs = []
    CMMT_data_min = []
    CMMT_data_max = []
    CMMT_modern_obs = []
    refs = []
    
    for info in sitedata:
        sites.append(info[0])
        lats.append(info[1])
        lons.append(info[2])
        WMMT_data_min.append(info[3])
        WMMT_data_max.append(info[4])
        CMMT_data_min.append(info[5])
        CMMT_data_max.append(info[6])
        WMMT_modern_obs.append(info[7])
        CMMT_modern_obs.append(info[8])
        refs.append(info[9])
   
    labels = []
    deg= u'\N{DEGREE SIGN}'
    for i, site in enumerate(sites):
  #     label = ''.join([c for c in site if c.isupper()])
        latstr = np.str(lats[i]) + deg + 'N'
        if lons[i] >180:
           lonstr = np.str((lons[i] - 360.) * -1.0) + deg +  'W'
        else:
           lonstr = np.str(lons[i]) + deg + 'E'
         
        label = site + ' (' +  latstr + ',' +  lonstr + ')'
        labels.append(site)
   
    return  (labels, lats, lons, np.asarray(WMMT_data_min),
             np.asarray(WMMT_data_max),  np.asarray(CMMT_data_min),
             np.asarray(CMMT_data_max),  np.asarray(WMMT_modern_obs),
             np.asarray(CMMT_modern_obs),refs)

 
def get_land_km5c():
    """
    these have been obtained from various sources so I am just typing them in
    """
    sitedata = []
    # site data is
    # sitename, sitelat, sitelon, min WMMT veg, max WMMT veg,
    # min WMMT beetle, max WMMT beetle, min CMMT veg, max CMMT veg
    # min CMMT beetle, max CMMT beetle, modern obs WMMT, modern obs CMMT
    # reference and date

    # lake baikal is from what ulrich sent me.
    # Lake E from Brigette-greeme mean july temp of +8 and average winter lows of 35degC

    sitedata.append(['Lake El\'gygytgyn', 67, 172, 15.0, 16.0,
                    -36.8, -30.4, 
                     8.0, np.nan,'CMMT 3.199Ma - 3.209Ma; Pavel Tarasov (pers. comm) \n WMMT Brigham-Grette et al. 2013'])
    sitedata.append(['Lake Baikal', 56, 108, 15.28, 17.52,
                     -1.67, 1.07, 
                     15.3, -17.4,'Km5c - unpublished \n (Method of Klage et al 2020)'])
    
   
    (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs) = land_reformat(sitedata)

    return (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs)

def get_land_LP():
    """
    these have been obtained from various sources so I am just typing them in
    LP means late pliocene
    """
    sitedata = []
    # lake baikal is from demske 2002.
    sitedata.append(['Lake Baikal', 56, 108, 13.0, 24.0, -15, 5,
                     np.nan, np.nan,'Prior to 3.5Ma (Demske et al 2002)'])
    # these are from popova et al 2012 using coexistence approach.
  
    sitedata.append(['Mirny', 55, 82, 18.8, 24.6, -0.3, 0.7,
                     np.nan, np.nan,'Popova et al 2012'])
    sitedata.append(['Merkutlinskiy', 56, 72, 17.3, 23.8, -3.8, 6.2,
                     np.nan, np.nan,' --"--'])
    sitedata.append(['Kabinet', 55, 80, 21.6, 24.4, -4.4, 4.6,
                     np.nan, np.nan,' --"--'])
    sitedata.append(['Delyankir', 63, 133, 18.9, 24.9, -6.9, 1.3,
                     np.nan, np.nan,' --"--'])
    sitedata.append(['Chernoluche', 55, 73, 19.6, 20.3, -5.9, 0.7,
                     np.nan, np.nan,' --"--'])
    sitedata.append(['Blizkiy', 64, 162, 15.6, 23.3, -12.8, 5.2,
                     np.nan, np.nan,' --"--'])
    sitedata.append(['42km', 55, 80, 21.6, 23.3, -4.4, 0.7,
                     np.nan, np.nan,' --"--'])
  
    sitedata.append(['Lost Chicken Mine', 64, 218, 12.0, 12.0, 
                     -2.0, -2.0, 15.3, -25.1, '2.9 +/- 0.4Ma: Ager et al. 1994'])
  
   
    (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs) = land_reformat(sitedata)

    return (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs)

def get_land_EP():
    """
    these have been obtained from various sources so I am just typing them in
    EP means early pliocene
    """
    sitedata = []
    # these are from popova et al 2012 using coexistence approach.
  
    sitedata.append(['Tnekveem', 66, 177, 18.9, 25.6, -11.8, 5.8, 
                     np.nan, np.nan,'Popova et al 2012'])
    sitedata.append(['Hydzhak', 63, 147, 18.8, 24.9, -8.7, 1.3,
                     np.nan, np.nan,' --"--'])

   
    (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs) = land_reformat(sitedata)

    return (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs)

def get_land_fletcher():
    """
    these have been obtained from various sources so I am just typing them in
    """
    sitedata = []
    # site data is
    # sitename, sitelat, sitelon, min WMMT data, max WMMT data,
    # min CMMT data, max CMMT data
    # modern obs WMMT, modern obs CMMT
    # we think FLB is 3.8Ma, BP originally was 3.4 but has been redated
   
    # for TF data we are reporting here temperature of warmest month
    # and temperature of warmest quater
    
    sitedata.append(['Near Meighen Island', 77.5, 261, [19.6, 12.8], 
                     [20.5, 13.3],
                    [-11.6, -6.8], [-11.4, -6.2],
                     4.1, -42.5, 'Fletcher et al. 2017'])
    sitedata.append(['Beaver Pond', 79, 278, [18.4, 12.4], [20.9, 13.1],
                     [-12.2, -7.3], [-11.5, -6.8],
                      7.1, -39.7, '3.9 +1.5 / -0.5Ma: Fletcher et al. 2017'])
    sitedata.append(['Fyles Leaf Beds', 79, 277, [19.7, 12.6], [21.1, 13.4],
                     [-12.8, -7.2], [-9.1, -5.5],
                      np.nan, np.nan, '3.8 +1/-0.7Ma: Fletcher et al. 2017'])
     
    (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs) = land_reformat(sitedata)


    
   
    return  (labels, lats, lons, WMMT_data_min,
             WMMT_data_max,  CMMT_data_min,
             CMMT_data_max,  WMMT_modern_obs,
             CMMT_modern_obs,refs)

def get_land_fletcher_ca():
    """
    these have been obtained from various sources so I am just typing them in
    """
    sitedata = []
    # site data is
    # sitename, sitelat, sitelon, min WMMT data, max WMMT data,
    # min CMMT data, max CMMT data
    # modern obs WMMT, modern obs CMMT
    # we think FLB is 3.8Ma, BP originally was 3.4 but has been redated
   
    # for tf data we are reporting hwere temp of warmest month and
    # warmest quater
    
    sitedata.append(['Near Meighen Island', 77.5, 261, [18.1, 10.6], 
                     [22.8, 16.2],
                    [-21.7,-16.3,], [-7.9, -2.7],
                     4.1, -42.5,'Fletcher et al. 2017'])
    sitedata.append(['Beaver Pond', 79, 278, [18.1, 11.3],[22.4, 16.3],
                     [-21.7, -15.0], [-8.1, -3.5],
                      7.1, -39.7,'3.9 +1.5 / -0.5Ma: Fletcher et al. 2017'])
    sitedata.append(['Fyles Leaf Beds', 79, 277, [18.1, 10.9], [22.7, 15.0],
                     [-16.9, -12.4], [-6.4, -2.3],
                      np.nan, np.nan,'3.8 +1/-0.7Ma: Fletcher et al. 2017'])
     
    (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs) = land_reformat(sitedata)


    
   
    return  (labels, lats, lons, WMMT_data_min,
             WMMT_data_max,  CMMT_data_min,
             CMMT_data_max,  WMMT_modern_obs,
             CMMT_modern_obs,refs)

def get_land_beetle():
    """
    these have been obtained from various sources so I am just typing them in
    """
    sitedata = []
    # site data is
    # sitename, sitelat, sitelon, min WMMT veg, max WMMT veg,
    # min WMMT beetle, max WMMT beetle, min CMMT veg, max CMMT veg
    # min CMMT beetle, max CMMT beetle, modern obs WMMT, modern obs CMMT

    # all were reported in Elias and Matthews
    # note Pliocene Tmin - modern Tmin varies between 2.3 degC and 20degC

    sitedata.append(['Ballast Brook', 74, 237, 
                     14.0, 14.5, -21.0, -19.5, 2.4, -41.4,'~3-5Ma: Flyes et al 1994'])
    sitedata.append(['Strathcona Beaver Peat', 79, 278, 
                     11.7, 12.2, -28.7, -27.2, 2.4, -41.4,'>3.3Ma: Matthews and Flyes 2000'])
    sitedata.append(['Near Meighen Island', 77.5, 261, 
                     11.5, 13.5, -33.0, -18.5, 4.1, -42.5,'~3Ma: Elias and Matthews 2002'])
    sitedata.append(['Lost Chicken Mine', 64, 218,  
                     13.5, 16.0, -27.75, -19.25, 15.3, -25.1,'~3Ma: Matthews and Telka 1997'])
    sitedata.append(['Bluefish', 67, 221,  
                     12.7, 15.0, -30.0, -20.5, 16.0, -29.0,'~LP: Matthews and Telka 1997'])
   
    (labels, lats, lons, WMMT_data_min, WMMT_data_max,  
     CMMT_data_min, CMMT_data_max, WMMT_modern_obs,
     CMMT_modern_obs,refs) = land_reformat(sitedata)


    
   
    return  (labels, lats, lons, WMMT_data_min,
             WMMT_data_max,  CMMT_data_min,
             CMMT_data_max,  WMMT_modern_obs,
             CMMT_modern_obs,refs)


def plot_figure(ax1,WMMT_data_min, WMMT_data_max, 
                CMMT_data_min, CMMT_data_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                labels,titleinfo,ystart,refs):
    """
    this subroutine tries to plot the figure for the paper which shows a nice
    DMC 
    """

    titlename = {'KM5c' : 'Near KM5c',
                 'LPCA' : 'Late Pliocene',
                 'EPCA' : 'Early Pliocene',
                 'TF-CRACLE' : 'Pliocene (Coexistence Likelihood Estimation)',
                 'TF-CA' : 'Pliocene (Coexistence Approach)',
                 'BA' : 'Beetle Assemblage data'
        }
    labelname = {0 : 'palaeodata WMMT',
                 1 : 'palaeodata CMMT'}

    WMMT_data_mean = (WMMT_data_min + WMMT_data_max) / 2.0
 
    CMMT_data_mean = (CMMT_data_min + CMMT_data_max) / 2.0
    
    print('ystart in the plot prog',ystart)
    if ystart == 0:
        print('here')
        ax1.get_xaxis().tick_top()
        ax1.axes.get_yaxis().set_visible(False)
        #ax.set_xlim([xmin, xmax])
        ax1.set_ylim([0, 32])
        plt.gca().invert_yaxis()
       

   
    nsites=len(WMMT_data_min)
    yarray = np.arange(ystart+1.5, nsites+ystart+1, 1)
    print('yarray',yarray)
    yend=ystart+nsites
 
    # plot warm and cold month temperature anomalies

    if titleinfo[0:2] != 'TF': # tamara fletchers data was not wmmt
        ax1.hlines(y=yarray, xmin=WMMT_data_min, xmax= WMMT_data_max,color='tab:red')
        ax1.hlines(y=yarray, xmin=CMMT_data_min, xmax= CMMT_data_max,color='tab:blue')
   
        plt.scatter(WMMT_data_mean, 
                    yarray, color='tab:red', marker='^',
                    label=labelname.get(ystart,None),s=60)
        plt.scatter(CMMT_data_mean, 
                    yarray, color='tab:blue', marker='^',
                    label = labelname.get(ystart+1,None) ,s=60)
    else:
        print(WMMT_data_mean)
        ax1.hlines(y=yarray, xmin=WMMT_data_min[:,0], xmax= WMMT_data_max[:,0],color='tab:red')
        ax1.hlines(y=yarray, xmin=CMMT_data_min[:,0], xmax= CMMT_data_max[:,0],color='tab:blue')
        ax1.hlines(y=yarray, xmin=WMMT_data_min[:,1], xmax= WMMT_data_max[:,1],color='tab:pink',linestyle='dashed')
        ax1.hlines(y=yarray, xmin=CMMT_data_min[:,1], xmax= CMMT_data_max[:,1],color='tab:cyan',linestyle='dashed')
        if titleinfo == 'TF-CRACLE':
            plt.scatter(WMMT_data_mean[:,0], 
                        yarray, color='tab:red', marker='v',
                        label='max T warmest month',s=60)
            plt.scatter(CMMT_data_mean[:,0], 
                        yarray, color='tab:blue', marker='v',
                        label =  'min T coldest month ',s=60)
            plt.scatter(WMMT_data_mean[:,1], 
                        yarray, color='tab:pink', marker='^',
                        label='wamest quarter T',s=60)
            plt.scatter(CMMT_data_mean[:,1], 
                        yarray, color='tab:cyan', marker='^',
                        label =  'coldest quarter T ',s=60)
        else:
            plt.scatter(WMMT_data_mean[:,0], 
                        yarray, color='tab:red', marker='v',s=60)
            plt.scatter(CMMT_data_mean[:,0], 
                        yarray, color='tab:blue', marker='v',s=60)
            plt.scatter(WMMT_data_mean[:,1], 
                        yarray, color='tab:pink', marker='^',s=60)
            plt.scatter(CMMT_data_mean[:,1], 
                        yarray, color='tab:cyan', marker='^',s=60)

    # put an arrow on LCM
    for i, label in enumerate(labels):
        if label[0:4] == 'Lost' and i > 4:
            plt.arrow(CMMT_data_max[i], yarray[i], -10, 0, linestyle='dotted', 
                  color='tab:blue', head_width=0.2, head_length=1.0)
   
    # try plotting axis
  
    for j in range(0, nsites):
        plt.text(-40.0, yarray[j], labels[j], ha='right')
#        if ystart == 0:
#            plt.text(27.0, yarray[j], refs[j], ha='left')
#        else:
#            plt.text(35.0, yarray[j], refs[j], ha='left')

  
  
    # plot individual models for pliocene
    for i,model in enumerate (MODELNAMES):
        if ystart==0:
            if i ==0: 
                plt.scatter(all_models_plio_WMT[:, i], yarray+0.2, 
                            color='tab:red', 
                            marker = 'x', s=25, label='HCM3 orig WMMT')
                plt.scatter(all_models_plio_CMT[:, i], yarray+0.2, 
                            color='tab:blue', 
                            marker = 'x', s=25, label='HCM3 orig CMMT')
            if i ==1: 
                plt.scatter(all_models_plio_WMT[:, i], yarray+0.2, 
                            color='tab:green', 
                            marker = 'o', s=25)
                plt.scatter(all_models_plio_CMT[:, i], yarray+0.2, 
                            color='tab:green', 
                            marker = 'o', s=25, label='HCM3 new CMMT')
        else:
            if i ==1:  #HadCM3 new
                plt.scatter(all_models_plio_WMT[:, i], yarray+0.2, 
                            color='tab:green', 
                            marker = 'o',s=25)
                plt.scatter(all_models_plio_CMT[:, i], yarray+0.2, 
                            color='tab:green', 
                            marker = 'o', s=25)

            else:
                plt.scatter(all_models_plio_WMT[:, i], yarray+0.2, 
                            color='tab:red', 
                            marker = 'x',s=25)
                plt.scatter(all_models_plio_CMT[:, i], yarray+0.2, 
                            color='tab:blue', 
                            marker = 'x', s=25)
        
    plt.hlines(y=ystart, xmin=-50., xmax=30., linewidth=0.5)
    plt.text(-30.0, -1.75, 'Cold Month Temperature (degC)', ha='left', color='tab:blue',fontsize=12)
    plt.text(5.0, -1.75, 'Warm Month Temperature (deg C)', ha='left', color='tab:red',fontsize=12)
    plt.text(-35.0,ystart+0.5,titlename.get(titleinfo),fontsize=12)
       
    
    #plt.legend(bbox_to_anchor=(0.5, -0.04), loc='lower center', ncol=3,
    #               prop = {'size':12})
    if ystart > 23:
       # handles, labels = plt.gca().get_legend_handles_labels()
        handles, labels = ax1.get_legend_handles_labels()
        order = [0,1,2,3,4,5,6,7,8] # 0: data WMMT
        print('handles',handles,'labels',labels)
        plt.legend([handles[idx] for idx in order],
                   [labels[idx] for idx in order],
               bbox_to_anchor=(0.5, -0.04), loc='lower center', ncol=3,
               prop = {'size':12})
   
   
 
    return ax1,yend+2

  
def get_model_data(land_lats, land_lons):
    """
    get the model data for these latitude and longitudes
    """
    
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


    return  (all_models_plio_WMT, all_models_plio_CMT,
     all_models_pi_WMT, all_models_pi_CMT, 
     all_models_CMMT_anom, 
     all_models_WMMT_anom, mmm_WMT, mmm_CMT, mmm_WMT_pi, mmm_CMT_pi)

def main():
    """
    calling structure
    a) get's model data
    b) get's proxy data
    c) plots model data with proxy data on top
    d) plots change in seasonal cycle at MI, BP and LB
    """

   
    # get land observations and cru temperature at km5c points 
    (sites, land_lats, land_lons, WMMT_data_min, WMMT_data_max,
     CMMT_data_min, CMMT_data_max, 
     WMMT_modern_obs, CMMT_modern_obs,refs) =  get_land_km5c()

    cru_min_temp, cru_max_temp = get_cru_temp(land_lats, land_lons)
    (all_models_plio_WMT, all_models_plio_CMT,
     all_models_pi_WMT, all_models_pi_CMT, 
     all_models_CMMT_anom,  all_models_WMMT_anom,  mmm_WMT, mmm_CMT, 
     mmm_WMT_pi, mmm_CMT_pi) = get_model_data(land_lats, land_lons)


   
    fig1 = plt.figure(figsize=[11.7, 11.7])
    ax1 = plt.axes(frameon=False)
   
    ax1,ystart = plot_figure(ax1,WMMT_data_min, WMMT_data_max, 
                CMMT_data_min, CMMT_data_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                sites,'KM5c',0,refs)

   
    #####################################################################
    # get land observations and cru temperature at lp points 
    (sites, land_lats, land_lons, WMMT_data_min, WMMT_data_max,
     CMMT_data_min, CMMT_data_max, 
     WMMT_modern_obs, CMMT_modern_obs,refs) =  get_land_LP()

    cru_min_temp, cru_max_temp = get_cru_temp(land_lats, land_lons)
    (all_models_plio_WMT, all_models_plio_CMT,
     all_models_pi_WMT, all_models_pi_CMT, 
     all_models_CMMT_anom,  all_models_WMMT_anom,  mmm_WMT, mmm_CMT, 
     mmm_WMT_pi, mmm_CMT_pi) = get_model_data(land_lats, land_lons)

    # plot data
    ax1,ystart = plot_figure(ax1,WMMT_data_min, WMMT_data_max, 
                CMMT_data_min, CMMT_data_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                sites,'LPCA',ystart,refs)
   
    #####################################################################
    # get land observations and cru temperature at ep points 
    (sites, land_lats, land_lons, WMMT_data_min, WMMT_data_max,
     CMMT_data_min, CMMT_data_max, 
     WMMT_modern_obs, CMMT_modern_obs,refs) =  get_land_EP()

    cru_min_temp, cru_max_temp = get_cru_temp(land_lats, land_lons)
    (all_models_plio_WMT, all_models_plio_CMT,
     all_models_pi_WMT, all_models_pi_CMT, 
     all_models_CMMT_anom,  all_models_WMMT_anom,  mmm_WMT, mmm_CMT, 
     mmm_WMT_pi, mmm_CMT_pi) = get_model_data(land_lats, land_lons)


   

    # plot data
    ax1,ystart = plot_figure(ax1,WMMT_data_min, WMMT_data_max, 
                CMMT_data_min, CMMT_data_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                sites,'EPCA',ystart,refs)
   

    ####################################################################
    # get land observations and cru temperature at early pliocene points 
    (sites, land_lats, land_lons, WMMT_data_min, WMMT_data_max,
     CMMT_data_min, CMMT_data_max, 
     WMMT_modern_obs, CMMT_modern_obs,refs) =  get_land_fletcher()

    cru_min_temp, cru_max_temp = get_cru_temp(land_lats, land_lons)
    (all_models_plio_WMT, all_models_plio_CMT,
     all_models_pi_WMT, all_models_pi_CMT, 
     all_models_CMMT_anom,  all_models_WMMT_anom,  mmm_WMT, mmm_CMT, 
     mmm_WMT_pi, mmm_CMT_pi) = get_model_data(land_lats, land_lons)

    # plot data

    ax1,ystart = plot_figure(ax1,WMMT_data_min, WMMT_data_max, 
                CMMT_data_min, CMMT_data_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                sites,'TF-CRACLE',ystart,refs)

     ####################################################################
    # get land observations and cru temperature at early pliocene points 
    (sites, land_lats, land_lons, WMMT_data_min, WMMT_data_max,
     CMMT_data_min, CMMT_data_max, 
     WMMT_modern_obs, CMMT_modern_obs,refs) =  get_land_fletcher_ca()

    cru_min_temp, cru_max_temp = get_cru_temp(land_lats, land_lons)
    (all_models_plio_WMT, all_models_plio_CMT,
     all_models_pi_WMT, all_models_pi_CMT, 
     all_models_CMMT_anom,  all_models_WMMT_anom,  mmm_WMT, mmm_CMT, 
     mmm_WMT_pi, mmm_CMT_pi) = get_model_data(land_lats, land_lons)

    # plot data

    for i,site in enumerate(sites):
        print('warm',site, WMMT_data_min[i], WMMT_data_max[i],mmm_WMT[i])
        print('cold',site, CMMT_data_min[i], CMMT_data_max[i],mmm_CMT[i])
        print(' ')
 #       for j, model in enumerate(MODELNAMES):
 #           print(model, all_models_plio_CMT[i,j], all_models_plio_CMT[i,j] - all_models_plio_CMT[1,j])
  #  sys.exit(0)

    ax1,ystart = plot_figure(ax1,WMMT_data_min, WMMT_data_max, 
                CMMT_data_min, CMMT_data_max, 
                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
                sites,'TF-CA',ystart,refs)



     ####################################################################
    # get land observations and cru temperature at early pliocene points 
#    (sites, land_lats, land_lons, WMMT_data_min, WMMT_data_max,
#     CMMT_data_min, CMMT_data_max, 
#     WMMT_modern_obs, CMMT_modern_obs,refs) =  get_land_beetle()

#    cru_min_temp, cru_max_temp = get_cru_temp(land_lats, land_lons)
#    (all_models_plio_WMT, all_models_plio_CMT,
#     all_models_pi_WMT, all_models_pi_CMT, 
#     all_models_CMMT_anom,  all_models_WMMT_anom,  mmm_WMT, mmm_CMT, 
#     mmm_WMT_pi, mmm_CMT_pi) = get_model_data(land_lats, land_lons)

    # plot data

#    ax1,ystart = plot_figure(WMMT_data_min, WMMT_data_max, 
#                CMMT_data_min, CMMT_data_max, 
#                all_models_plio_WMT, all_models_plio_CMT, mmm_WMT, mmm_CMT,
#                sites,'BA',ystart,refs)
#    ax1.set_ylim(None,0)    

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'newparams/seasonal_dmc_plot.eps')
    plt.savefig(fileout)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
                   'newparams/seasonal_dmc_plot.png')
    plt.savefig(fileout)
    plt.close()

 
##########################################################
# main program


LINUX_WIN = 'l'
FILESTART = '/nfs/hera1/earjcti/'
DATABASE = '/nfs/hera1/pliomip2/data/'

MODELNAMES = ['HadCM3','HadCM3_new']


LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')

main()

#sys.exit(0)
