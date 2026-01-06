#!/usr/bin/env python3
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
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.quickplot as qplt
import iris.plot as iplt
import cartopy.crs as ccrs

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


def get_single_model(model, latreq, lonreq, exptid):
    """
    read in the pliocene data from 'model'  return the temperatures
    at the list of sites
    """

    filename = ('/nfs/hera1/earjcti/regridded100/' + model +
                '/' + exptid + '.NearSurfaceTemperature.mean_month.nc')
  
    print(filename)
    plio_cube = iris.load_cube(filename)
   
    nsites = len(latreq)
    plio_mean_array = np.zeros(nsites)
    plio_min_array = np.zeros(nsites)
    plio_max_array = np.zeros(nsites)
   
    for i in range(0,nsites):
        # modellon is whole numbers from 0-360
        # lat is half numbers from -89.5 to 89.5

        modlon = np.around(lonreq[i])
        if modlon < 0: modlon = modlon + 360.

        lat_ix = ((np.abs(plio_cube.coord('latitude').points 
                         - latreq[i])).argmin())
        lon_ix = ((np.abs(plio_cube.coord('longitude').points 
                         - modlon)).argmin())
    
        plio_array = plio_cube.data[:, lat_ix, lon_ix]
        plio_min_array[i] = np.min(plio_array)
        plio_max_array[i] = np.max(plio_array)
        plio_mean_array[i] = np.mean(plio_array)

   
    return plio_mean_array, plio_max_array, plio_min_array

 
def get_land_obs():
    """
    reads in the spredsheet from ulrich and returns temperatures
    """

    dfs = pd.read_excel(LAND_DATAFILE)
    sites = []
    lats = []
    lons = []
    temps = []
    temp_modern = []
    temp_uncert = []

    row_locs = [2, 3, 4, 5, 6, 7, 8, 9, 11, 12]
    for rl in row_locs:
        # if temp ne nan then move to array
        temp = dfs.iloc[rl, 9]
        
        print(temp,'julia')
        if np.isfinite(temp):
            sites.append(dfs.iloc[rl, 0])
            lats.append(dfs.iloc[rl, 2])
            lons.append(dfs.iloc[rl, 3])
            temp_modern.append(dfs.iloc[rl, 4])
            temp_uncert.append(dfs.iloc[rl,10])
            temps.append(temp)

    print(temp_uncert)
    for i, temp in enumerate(temp_uncert):
        if i > 0:
            temp2 = temp[2:]
        else:
            temp2=0.0
        print(temp, temp2)
        temp_uncert[i]=np.float(temp2)
     
    labels = []
    deg= u'\N{DEGREE SIGN}'
    for i, site in enumerate(sites):
        label = ''.join([c for c in site if c.isupper()])
        if lats[i] < 0:
            latstr = np.str(np.int(np.round(lats[i] * -1.0, 0))) + deg +  'S'
        else:
            latstr = np.str(np.int(np.around(lats[i], 0))) + deg + 'N'
        if lons[i] < 0:
            lonstr = np.str(np.int(np.round(lons[i] * -1.0, 0))) + deg +  'W'
        else:
            lonstr = np.str(np.int(np.around(lons[i], 0))) + deg + 'E'
        
        label = site + '\n (' +  latstr + ',' +  lonstr + ')'
        labels.append(label)
   
    return lats, lons, temps, temp_modern, temp_uncert, labels

 

def plot_figure(plio_temp_obs, plio_model_400, plio_model_450, labels, ax, fig):
    """
    this subroutine tries to plot the figure for the paper which shows a nice
    DMC 
    """


    #ax1 = ax.axes(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_yaxis().set_ticks([])
    #ax.axes.get_yaxis().set_visible(False)
    #ax.get_xaxis().tick_bottom()
   

    yarray = np.arange(1, len(plio_temp_obs) + 1, 1)

 
   # try plotting model data anomaly
    plt.vlines(x=0, ymin=-0, ymax=9, linewidth=0.5)
   
    # plot individual models for pliocene
    model_400_anom = np.zeros(np.shape(plio_model_400))
    model_450_anom = np.zeros(np.shape(plio_model_450))
    colors = ['black','green','orange']
    for i, model in enumerate(MODELNAMES):
        model_400_anom[:, i] = plio_model_400[:, i] - plio_temp_obs  
        model_450_anom[:, i] = plio_model_450[:, i] - plio_temp_obs  
        plt.scatter(model_400_anom[:, i], yarray, marker = 'o', 
                        color = colors[i], s=30)
        plt.scatter(model_450_anom[:, i], yarray-0.2, marker = '^',
                       color = colors[i], s=30)
       
  
    # add site labels
    plt.text(-5.0, yarray[7], labels[7], ha='right')
    for j in range(0, 7):
        plt.text(np.min(model_400_anom[j, :]) - 1.0,
                 yarray[j], labels[j], ha='right')

    plt.xlabel('Temperature difference from observations (deg C)')
    #fig.legend(loc = 'center left')
    plt.title('a) Annual Mean Temperature', loc='left')

  
    plt.xlim(-25, 7.5)
    plt.ylim(9, 0)
   
def get_land_warm_cold():
    """
    these have been obtained from various sources so I am just typing them in
    """
    sitedata = []
    # site data is
    # sitename, sitelat, sitelon, min WMMT veg, max WMMT veg,
    # min WMMT beetle, max WMMT beetle, min CMMT veg, max CMMT veg
    # min CMMT beetle, max CMMT beetle
    sitedata.append(['Meighen Island', 80, 261, 19.6, 20.5,
                     11.5, 13.5, -11.6, -11.4,
                     -33.0, -18.5])
    sitedata.append(['Beaver Pond', 79, 278, 18.4, 20.9,
                     np.nan, np.nan, -12.2, -11.5, np.nan, np.nan])
    sitedata.append(['Lake El\'gygytgyn', 67, 172, 15.0, 16.0,
                    np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    sitedata.append(['Lost Chicken Mine', 64, 142, 12.0, 12.0, 
                     13.5, 16.0, -2.0, -2.0, -27.75, -19.25])
    sitedata.append(['Lake Baikal', 56, 108, 15.28, 17.52,
                     np.nan, np.nan, -1.67, 1.07, np.nan, np.nan])
                     

    sites = []
    lats = []
    lons = []
    WMMT_veg_min = []
    WMMT_veg_max = []
    WMMT_beetle_min = []
    WMMT_beetle_max = []
    CMMT_veg_min = []
    CMMT_veg_max = []
    CMMT_beetle_min = []
    CMMT_beetle_max = []
    
    for info in sitedata:
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
             np.asarray(CMMT_beetle_max))


def plot_seas_fig(veg_temp, beetle_temp, plio_model_400,
                     plio_model_450, labels, ax, fig, wc_ind):
    """
    this subroutine tries to plot the figure for the paper which shows how 
    different values of CO2 affect the seasonal anomaly
    """

    nmods = len(veg_temp)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_yaxis().set_ticks([])
   
    yarray = np.arange(1, nmods + 1, 1)

 
   # try plotting model data anomaly
    plt.vlines(x=0, ymin=-0, ymax=9, linewidth=0.5)
    plt.hlines(y=8.95, xmin=-25, xmax=7.5, linewidth=0.5)
   
    # plot individual models for pliocene
    model_400_anom = np.zeros(np.shape(plio_model_400))
    model_450_anom = np.zeros(np.shape(plio_model_450))
    model_450_beetle_anom = np.zeros(np.shape(plio_model_450))
    colors = ['black','green','orange']
    for i, model in enumerate(MODELNAMES):
       
        model_400_anom[:, i] = plio_model_400[:, i] - veg_temp 
        model_450_anom[:, i] = plio_model_450[:, i] - veg_temp  
        model_450_beetle_anom[:, i] = plio_model_450[:, i] - beetle_temp
        ax.scatter(model_400_anom[:, i], yarray - 0.2, marker = 'o', 
                        color = colors[i], s=30, label=model + ' 400ppmv - vegdata')
        ax.scatter(model_450_anom[:, i], yarray, marker = '^',
                       color = colors[i], s=30, label = model + ' 450ppmv - vegdata')
        ax.scatter(model_450_beetle_anom[:, i], yarray + 0.2, marker = 's',
                    color = colors[i], s=30, label= model + ' 450ppmv - beetledata')
   
    # add site labels
    if wc_ind == 'c':
        for j in [0, 1, 3, 4]:
            plt.text(0.5,
                     yarray[j], labels[j], ha='left')
        plt.title('c) Cold Month Temperature', loc='left')
        print(labels[j], model_400_anom[j, :], model_450_anom[j, :], model_450_anom[j, :] -  model_400_anom[j, :], MODELNAMES)

    else:
        for j in range(0, 5):
            plt.text(-20, yarray[j], labels[j], ha='right')
        plt.title('b) Warm Month Temperature', loc='left')



        
    if wc_ind == 'c':
        plt.xlabel('Temperature difference from observations (deg C)')
        handles, labs = fig.gca().get_legend_handles_labels()
        order = [0, 3, 6, 1, 4, 7, 2, 5, 8, ]
        fig.legend([handles[i] for i in order], 
                   [labs[i] for i in order],
                   loc = 'center left')
       # fig.legend(loc='center left')
  
    plt.xlim(-35, 15)
    plt.ylim(nmods + 1, 0)
   
  

  
def main():
    """
    calling structure
    a) get's model data
    b) get's proxy data
    c) plots model data with proxy data on top
    """

   
    # get land observations and cru temperature at land points
    
    (land_lats, land_lons, land_temp, 
     modern_temp, plio_unc, land_labels)= get_land_obs()
   
    
    # get ind models data
    all_models_plio_400 = np.zeros((len(land_lons), len(MODELNAMES)))
    all_models_plio_450 = np.zeros((len(land_lons), len(MODELNAMES)))
    for i, model in enumerate(MODELNAMES):
        (ind400, ind400WMT, 
         ind400CMT) = get_single_model(model, land_lats, land_lons, 'EOI400')
        all_models_plio_400[:, i] = ind400
        
        (ind450, ind450WMT,
         ind450CMT) = get_single_model(model, land_lats, land_lons, 'EOI450')
        all_models_plio_450[:, i] = ind450
        

  
    # get warm month and cold month temperatures from data
    (sites, land_lats, land_lons, WMMT_veg_min, WMMT_veg_max, WMMT_beetle_min,
     WMMT_beetle_max, CMMT_veg_min, CMMT_veg_max, CMMT_beetle_min,
     CMMT_beetle_max) =  get_land_warm_cold()

    # get warm month and cold month temperatures from model
    allmod_plio_400_wmt = np.zeros((len(land_lons), len(MODELNAMES)))
    allmod_plio_450_wmt = np.zeros((len(land_lons), len(MODELNAMES)))
    allmod_plio_400_cmt = np.zeros((len(land_lons), len(MODELNAMES)))
    allmod_plio_450_cmt = np.zeros((len(land_lons), len(MODELNAMES)))
   
    for i, model in enumerate(MODELNAMES):
        (ind400, ind400WMT, 
         ind400CMT) = get_single_model(model, land_lats, land_lons, 'EOI400')
        allmod_plio_400_wmt[:, i] = ind400WMT
        allmod_plio_400_cmt[:, i] = ind400CMT
        
        (ind450, ind450WMT,
         ind450CMT) = get_single_model(model, land_lats, land_lons, 'EOI450')
        allmod_plio_450_wmt[:, i] = ind450WMT
        allmod_plio_450_cmt[:, i] = ind450CMT
  

    # plot figure annual mean
    fig1 = plt.figure(figsize=[12.0, 12.0], constrained_layout=True)
    gs = gridspec.GridSpec(nrows=2, ncols=2)

    ax1 = fig1.add_subplot(gs[:,0])
    plot_figure(land_temp,  all_models_plio_400, all_models_plio_450, 
                land_labels, ax1, fig1) 

  
    ax2 = fig1.add_subplot(gs[0,1])
    plot_seas_fig((WMMT_veg_min + WMMT_veg_max) / 2.0, 
                  (WMMT_beetle_min + WMMT_beetle_max) / 2.0,
                  allmod_plio_400_wmt, allmod_plio_450_wmt, sites, 
                  ax2, fig1, 'w')

    ax3 = fig1.add_subplot(gs[1,1])
    plot_seas_fig((CMMT_veg_min + CMMT_veg_max) / 2.0, 
                  (CMMT_beetle_min + CMMT_beetle_max) / 2.0,
                  allmod_plio_400_cmt, allmod_plio_450_cmt, sites, ax3, 
                  fig1, 'c')

  
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/CO2_uncertainty.eps')
    plt.savefig(fileout)


##########################################################
# main program


LINUX_WIN = 'l'
FILESTART = '/nfs/hera1/earjcti/'

MODELNAMES = ['COSMOS', 'HadCM3', 'MIROC4m']

LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')

main()

#sys.exit(0)
