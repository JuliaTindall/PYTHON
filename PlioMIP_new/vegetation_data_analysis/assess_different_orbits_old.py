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




def get_single_model(expt, latreq, lonreq):
    """
    read in the pliocene data from 'model'  return the temperatures
    at the list of sites
    """
    
    nsites = len(latreq)
    months = ['January','February','March','April','May','June','July','August','September','October','November','December']
    monthalt = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    seas_field = np.zeros((len(months), nsites))

    plio_mean_array = np.zeros(nsites)
    plio_min_array = np.zeros(nsites)
    plio_max_array = np.zeros(nsites)
   
  
 
    # can change this to calendar corrected files
    for j, month in enumerate(months):
        if expt in ('xozzb', 'xozzc','xozzd','xozzf','xibol'):
            filename = ('/nfs/hera1/earjcti/um/' + expt +
                        '/database_averages/' + expt + 
                        '_Monthly_Average_' + month 
                        + '_a@pd_Temperature.nc')
            cubetemp = iris.load_cube(filename)
      
        else:
            filename = ('/nfs/hera1/earjcti/um/'+ expt + 
                        '/cal_cor/' + expt + 'a@pa_avg'+ monthalt[j] + '.nc')
            cube2 = iris.load_cube(filename, 
                                      'SURFACE TEMPERATURE AFTER TIMESTEP')
            cubetemp = cube2.collapsed(['t'], iris.analysis.MEAN)
        cube = iris.util.squeeze(cubetemp)
       
        for i in range(0,nsites):
            # modellon is whole numbers from 0-360
            # lat is half numbers from -89.5 to 89.5

            modlon = np.around(lonreq[i])
            if modlon < 0: modlon = modlon + 360.

            lat_ix = ((np.abs(cube.coord('latitude').points 
                              - latreq[i])).argmin())
            lon_ix = ((np.abs(cube.coord('longitude').points 
                         - modlon)).argmin())

            seas_field[j, i] = cube.data[lat_ix, lon_ix]
    
    for i in range(0, nsites):
        plio_min_array[i] = np.min(seas_field[:, i]) - 273.15
        plio_max_array[i] = np.max(seas_field[:, i]) - 273.15
        plio_mean_array[i] = np.mean(seas_field[:, i]) - 273.15

   
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

 

def plot_figure(plio_temp_obs, plio_model, labels, ax, fig):
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
    model_anom = np.zeros(np.shape(plio_model))
    colors = ['blue','orange','green', 'red', 'black','purple']
  
    for i, model in enumerate(EXPTNAMES):
        model_anom[:, i] = plio_model[:, i] - plio_temp_obs  
        plt.scatter(model_anom[:, i], yarray - (0.01 * i), marker = 'o', 
                    color=colors[i],s=30)
       
  
    # add site labels
    plt.text(-5.0, yarray[7], labels[7], ha='right')
    for j in range(0, 7):
        plt.text(np.min(model_anom[j, :]) - 1.0,
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


def plot_seas_fig(veg_temp, beetle_temp, plio_model,
                     labels, ax, fig, wc_ind):
    """
    this subroutine tries to plot the figure for the paper which shows how 
    different values of orbit affect the seasonal anomaly
    """
    print(plio_model[4,:],labels[4],wc_ind)

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
    model_anom = np.zeros(np.shape(plio_model))
    model_beetle_anom = np.zeros(np.shape(plio_model))
    colors = ['blue','green','orange', 'red', 'black','purple']
    for i, model in enumerate(EXPTNAMES):
       
        model_anom[:, i] = plio_model[:, i] - veg_temp 
        model_beetle_anom[:, i] = plio_model[:, i] - beetle_temp
        if wc_ind == 'c':
            ax.scatter(model_anom[:, i], yarray - 0.2 - (0.01 * i) , 
                       marker = 'o', 
                       color = colors[i],
                       s=30, label=PERIOD.get(model))
        else:
            ax.scatter(model_anom[:, i], yarray - 0.2- (0.01 * i),
                       marker = 'o', 
                       color = colors[i], s=30)
        ax.scatter(model_beetle_anom[:, i], yarray + 0.2- (0.01 * i),
                   marker = 's',
                   color = colors[i], 
                   s=30)
   
    # add site labels
    if wc_ind == 'c':
        for j in [0, 1, 3, 4]:
            plt.text(0.5,
                     yarray[j], labels[j], ha='left')
        plt.title('c) Cold Month Temperature', loc='left')
        print(labels[j], model_anom[j, :],  EXPTNAMES)


    else:
        for j in range(0, 5):
            plt.text(-20, yarray[j], labels[j], ha='right')
        plt.title('b) Warm Month Temperature', loc='left')

        
    


        
    if wc_ind == 'c':
        plt.xlabel('Temperature difference from observations (deg C)')
        #handles, labs = fig.gca().get_legend_handles_labels()
        #order = [0, 3, 6, 1, 4, 7, 2, 5, 8, ]
        #fig.legend([handles[i] for i in order], 
        #           [labs[i] for i in order],
        #           loc = 'center left')
        fig.legend(loc='center left')
  
    #plt.xlim(-35, 15)
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
    all_models_plio = np.zeros((len(land_lons), len(EXPTNAMES)))
    for i, model in enumerate(EXPTNAMES):
        (ind, indWMT, 
         indCMT) = get_single_model(model, land_lats, land_lons)
        all_models_plio[:, i] = ind
        
        

  
    # get warm month and cold month temperatures from data
    (sites, land_lats, land_lons, WMMT_veg_min, WMMT_veg_max, WMMT_beetle_min,
     WMMT_beetle_max, CMMT_veg_min, CMMT_veg_max, CMMT_beetle_min,
     CMMT_beetle_max) =  get_land_warm_cold()

    # get warm month and cold month temperatures from model
    allmod_plio_wmt = np.zeros((len(land_lons), len(EXPTNAMES)))
    allmod_plio_cmt = np.zeros((len(land_lons), len(EXPTNAMES)))
   
    for i, model in enumerate(EXPTNAMES):
        (ind, indWMT, 
         indCMT) = get_single_model(model, land_lats, land_lons)
        allmod_plio_wmt[:, i] = indWMT
        allmod_plio_cmt[:, i] = indCMT
        
     

    # plot figure annual mean
    fig1 = plt.figure(figsize=[12.0, 12.0], constrained_layout=True)
    gs = gridspec.GridSpec(nrows=2, ncols=2)

    ax1 = fig1.add_subplot(gs[:,0])
    plot_figure(land_temp,  all_models_plio, 
                land_labels, ax1, fig1) 

   
    ax2 = fig1.add_subplot(gs[0,1])
    plot_seas_fig((WMMT_veg_min + WMMT_veg_max) / 2.0, 
                  (WMMT_beetle_min + WMMT_beetle_max) / 2.0,
                  allmod_plio_wmt, sites, 
                  ax2, fig1, 'w')

    ax3 = fig1.add_subplot(gs[1,1])
    plot_seas_fig((CMMT_veg_min + CMMT_veg_max) / 2.0, 
                  (CMMT_beetle_min + CMMT_beetle_max) / 2.0,
                  allmod_plio_cmt, sites, ax3, 
                  fig1, 'c')

  
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/orbital_uncertainty.eps')
    plt.savefig(fileout)


##########################################################
# main program


LINUX_WIN = 'l'
FILESTART = '/nfs/hera1/earjcti/'

EXPTNAMES = ['xibol', 'xozzf', 'xozzd', 'xozze', 'xozzc']

#EXPTNAMES = ['xogzc', 'xogzb']
PERIOD = {'xiboi' : 'pi (old)',
          'xibol' : 'Km5c(old)',
          'xozzc' : 'K1 (3.0560)',
          'xozzd' : 'G17 (2.950)',
          'xozze' : 'KM3 (3.155)',
          'xozzf' : '    (3.053)'}


LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')

main()

#sys.exit(0)
