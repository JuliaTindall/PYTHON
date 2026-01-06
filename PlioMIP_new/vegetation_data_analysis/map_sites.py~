#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on December 2020
# 
#  This program will plot a map showing the amplitude of the seasonal cycle 
#  over Land for the pliocene and the preindustrial
#
#
#import os
import numpy as np
import pandas as pd
#import scipy as sp
#import cf
import iris
#import iris.util
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import netCDF4
from mpl_toolkits.basemap import Basemap, shiftgrid
#from netCDF4 import Dataset, MFDataset
#import iris.analysis.cartography
#import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
#import cf_units as unit
#from iris.experimental.equalise_cubes import equalise_attributes
import cartopy
import cartopy.crs as ccrs
#import matplotlib.ticker as mticker
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from mpl_toolkits.basemap import Basemap

import sys

 




def mask_ocean_for_cube(cube, lsmfile):
    """
    get's the land sea mask and sets the cube value to missing where the  
    lsm is 0

    """
    lsm_cube = iris.load_cube(lsmfile)
    cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
    
    regrid_lsm_cube = lsm_cube.regrid(cubegrid, iris.analysis.Linear())
    
    for i, lsm in enumerate(regrid_lsm_cube.coord('longitude').points):
        if lsm != cube.coord('longitude').points[i]:
            print('lsm lon error')
            print(lsm, cube.coord('longitude').points[i])
            sys.exit(0)
       
    for j, lsm in enumerate(regrid_lsm_cube.coord('latitude').points):
         if lsm != cube.coord('latitude').points[j]:
            print('lsm lat error')
            print(lsm, cube.coord('latitude').points[i])
            sys.exit(0)
    
    masked_cube = iris.util.mask_cube(cube, np.where(regrid_lsm_cube.data == 0))
    
    print(masked_cube.data)
        
    return masked_cube
        
  

def get_max_seascyc(fieldname):
    """
    gets the warm month temperature minus the cold month temperature
    """

    temp_cube = iris.load_cube(FILENAME, fieldname)
    warm_cube = temp_cube.collapsed(['time'], iris.analysis.MAX)
    cold_cube = temp_cube.collapsed(['time'], iris.analysis.MIN)

    ann_mean_cube = temp_cube.collapsed(['time'], iris.analysis.MEAN)
    seas_cyc_cube = warm_cube - cold_cube


    return seas_cyc_cube, ann_mean_cube, warm_cube, cold_cube

def plotcube(cube, period, outname, V, ax):
    """
    plots the cube to a file
    """

    ax.coastlines()
    ax.set_extent([-180, 180, 30, 80], ccrs.PlateCarree())
   
    cs = iplt.contourf(cube, V, extend = 'both')
    plt.title(period + ': ' + outname)
    if period == 'mPWP':
        cbar = plt.colorbar(cs, orientation = 'horizontal')
        cbar.set_label('degC')
        
   

def plot_locations2(match_cube, alt_seasoncube, precip_cube,
                    matchno, subscript):
    """
    we have a cube (match_cube) where the temperature has to match the matchno.
    This could be the warm_cube matching 20degC, or the cold_cube matching -10degC
    We find all the locations in match_cube which are close to match number
    We then plot all the locations (and their temperatures in the alt_seasoncube).  The alternative-season_cube could be either the warmcube or the cold cube
    """  

    lons = {'Lake Baikal' : 108., 'Lake Baikal M' : 108.,
            'Lost Chicken Mine' : -142., 'Lost Chicken Mine M' : -142.,
            'Lost Chicken Mine B' : -142.}
    lats = {'Lake Baikal' : 56., 'Lake Baikal M' : 56.,
            'Lost Chicken Mine' : 64., 'Lost Chicken Mine M' : 64.,
            'Lost Chicken Mine B' : 64.}
  
    
    
    lon = lons.get(SITE)
    lat = lats.get(SITE)
    lonalt = lon
    if lonalt > 180: lonalt = lon - 360.
   

    treqmin = matchno - 1.0
    treqmax = matchno + 1.0
    
  
    cube2 = match_cube.copy()
    cube3 = alt_seasoncube.copy()
    for index, point in np.ndenumerate(cube2.data):
        if treqmax < point  or treqmin > point:
            cube2.data[index] = 0
    for i, loncube in enumerate(cube2.coord('longitude').points):
        if lon  < loncube  or loncube < lon - 50:
            cube2.data[:, i] = 0
    for j, latcube in enumerate(cube2.coord('latitude').points):
        if 85 < latcube  or latcube < lat - 30:
            cube2.data[j, :] = 0
       
    masked_cube = iris.util.mask_cube(cube2, 
                                      np.where(cube2.data == 0))   
    masked_alt = iris.util.mask_cube(cube3, 
                                      np.ma.where(masked_cube.data.mask == 1))
            
    
    # plot locations and winter temperatures.
    
    ax1 = plt.subplot(221, projection = ccrs.PlateCarree())
    ax1.coastlines()
    ax1.set_extent([np.max([-180., lonalt - 50.]), lonalt + 50., lat-30., 85.], ccrs.PlateCarree())
    cs=iplt.contourf(masked_alt, extend='both')
    cbar = plt.colorbar(cs, orientation = 'horizontal')
    plt.text(lonalt, lat, 'x', fontsize=20)
    plt.text(lon+5.0, lat, SITE)        
    plt.title(subscript)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/' + SITE.replace(' ','_') + '_' + subscript + '.eps')
    
    # plot an histogram of the temperatures
    allvals = []
    precipvals = []
    for index, point in np.ndenumerate(masked_alt.data):
        val = masked_alt.data[index]
        if np.isfinite(val):
            allvals.append(val)
            precipvals.append(precip_cube.data[index] * 360.)
    ax2 = plt.subplot(222)
    num_bins = 10
    n, bins, patches = plt.hist(allvals, num_bins, alpha=0.5)
   
    # plot a scatter plot of the cold month temperatures vs the precip
    ax3 = plt.subplot(223)
    plt.scatter(allvals,precipvals)
    plt.xlabel('cold month temperature')
    plt.ylabel('precipitation')
   
    plt.tight_layout()
    plt.savefig(fileout)
    plt.close()



def plot_locations_bothmatch(warm_cube, cold_cube):
    """
    Is there a location where the warm temperature matches the warm cube
    and the cold temperature matches the cold cube
    """  

    lons = {'Lake Baikal' : 108., 'Lake Baikal M' : 108.,
            'Lost Chicken Mine' : -142., 'Lost Chicken Mine M' : -142.,
            'Lost Chicken Mine B' : -142.}
    lats = {'Lake Baikal' : 56., 'Lake Baikal M' : 56.,
            'Lost Chicken Mine' : 64., 'Lost Chicken Mine M' : 64.,
            'Lost Chicken Mine B' : 64.}
  
    
    lon = lons.get(SITE)
    lat = lats.get(SITE)
    lonalt = lon
    if lonalt > 180: lonalt = lon - 360.
   
    precision = 3.0
    locs_lon = []
    locs_lat = []
    warm_temp = []
    cold_temp = []
    
    cube2 = warm_cube.copy()
    data2 = cube2.data
    cube3 = cold_cube.copy()
    data3 = cube3.data
    for index, point_warm in np.ndenumerate(data2):
        point_cold = data3[index]
        
        if ((WARM_TEMP_REQ + precision  > point_warm 
        and WARM_TEMP_REQ - precision  < point_warm
        and COLD_TEMP_REQ + precision > point_cold
        and COLD_TEMP_REQ - precision < point_cold) and not data2.mask[index]):
            locs_lon.append(cube2.coord('longitude').points[index[1]])
            locs_lat.append(cube2.coord('latitude').points[index[0]])
            warm_temp.append(point_warm)
            cold_temp.append(point_cold)
        else:
            cube2.data[index] = 0
    
    for i, loncube in enumerate(cube2.coord('longitude').points):
        if lon  < loncube  or loncube < lon - 50:
            cube2.data[:, i] = 0
    for j, latcube in enumerate(cube2.coord('latitude').points):
        if 85 < latcube  or latcube < lat - 30:
            cube2.data[j, :] = 0
       
    for i, lon in enumerate(locs_lon):
        print(lon, locs_lat[i], warm_temp[i], cold_temp[i])
    print('MEAN', np.mean(np.asarray(warm_temp)), np.mean(np.asarray(cold_temp)), np.mean(np.asarray(warm_temp) - np.asarray(cold_temp)))
  
    # plot locations and winter temperatures.
    
    ax1 = plt.subplot(111, projection = ccrs.PlateCarree())
    ax1.coastlines()
    #ax1.set_extent([np.max([-180., lonalt - 50.]), lonalt + 50., lat-30., 85.], ccrs.PlateCarree())
    ax1.set_extent([-180, 180, -90, 90], ccrs.PlateCarree())
    #cs=iplt.contourf(masked_alt, extend='both')
    plt.scatter(locs_lon, locs_lat)
    #cbar = plt.colorbar(cs, orientation = 'horizontal')
    plt.text(lonalt, lat, 'x', fontsize=20)
    plt.text(lonalt-20.0, lat-5, SITE)        
   # plt.title(subscript)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/' + SITE.replace(' ','_') + '_bothmatch.eps')
    
    plt.tight_layout()
    plt.savefig(fileout)
    plt.close()

  
                                      

def plot_locations(cube, warm_cube, cold_cube):
    """
    we have a value for the required annual cycle.  We will plot all locations
    which have this annual cycle
    """  

    lons = {'Lake Baikal' : 108., 'Lake Baikal M' : 108.,
            'Lost Chicken Mine' : -142., 'Lost Chicken Mine M' : -142.,
            'Lost Chicken Mine B' : -142.}
    lats = {'Lake Baikal' : 56., 'Lake Baikal M' : 56.,
            'Lost Chicken Mine' : 64., 'Lost Chicken Mine M' : 64.,
            'Lost Chicken Mine B' : 64.}
    treqmin = ANN_CYC_REQ - 1.0
    treqmax = ANN_CYC_REQ + 1.0
    
  
    cube2 = cube.copy()
    for index, point in np.ndenumerate(cube2.data):
        if treqmax < point  or treqmin > point:
            cube2.data[index] = 0
           
    masked_cube = iris.util.mask_cube(cube2, 
                                      np.where(cube2.data == 0))
    masked_warm = iris.util.mask_cube(warm_cube, 
                                      np.ma.where(cube2.data.mask == 1))
    masked_cold = iris.util.mask_cube(cold_cube, 
                                      np.ma.where(cube2.data.mask == 1))
    
    # plot
    lon = lons.get(SITE)
    lat = lats.get(SITE)

    ax1 = plt.subplot(221, projection = ccrs.PlateCarree())
    ax1.coastlines()
    ax1.set_extent([np.max([lon - 50., -180.]), lon + 50., lat-50., 85.], ccrs.PlateCarree())
    #V = np.arange(13, 28, 3)
    cs=iplt.contourf(masked_warm)
    cbar = plt.colorbar(cs, orientation = 'horizontal')
    plt.text(lon, lat, 'x', fontsize=20)
    plt.text(lon+5.0, lat, SITE)
   
    ax2 = plt.subplot(222, projection = ccrs.PlateCarree())
    ax2.coastlines()
    ax2.set_extent([lon - 50, lon + 50, lat-50, 85], ccrs.PlateCarree())
    #V = np.arange(-15, 10, 5)
    cs=iplt.contourf(masked_cold)
    cbar = plt.colorbar(cs, orientation = 'horizontal')
    plt.text(lon, lat, 'x', fontsize=20)
    plt.text(lon+5.0, lat, SITE)
           
    plt.title('Annual cycle amplitude ~ ' + np.str(np.int(ANN_CYC_REQ)) + 'degC')
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/' + SITE.replace(' ','_') + '_same_seas_cyc.eps')

    plt.show()
    plt.close()
                                      

  
def main():
    """
    calling structure
    1. read in data
    2. get amplitude(warm month T - cold month T) and annual mean temperature 
    3. filter out land sea mask
    4. plot amplitude and annmeantemp
    """

    # get pliocene stuff
    (plio_amplitude_cube,
     plio_mean_cube,
     plio_warm_cube,
     plio_cold_cube)= get_max_seascyc('NearSurfaceTemperaturemean_plio')
    lsmfile = BCSTART + '/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    plio_landamp_cube  = mask_ocean_for_cube(plio_amplitude_cube, lsmfile)
    plio_landmean_cube = mask_ocean_for_cube(plio_mean_cube, lsmfile)
   
    # get pi stuff
    (pi_amplitude_cube,
     pi_mean_cube,
     pi_warm_cube,
     pi_cold_cube)= get_max_seascyc('NearSurfaceTemperaturemean_pi')
    pi_precip_cube = iris.load_cube(PRECIPFILE, 'TotalPrecipitationmean_pi')
    lsmfile = BCSTART + '/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc'
    pi_landamp_cube= mask_ocean_for_cube(pi_amplitude_cube, lsmfile)
    pi_landmean_cube = mask_ocean_for_cube(pi_mean_cube, lsmfile)
    pi_warm_cube = mask_ocean_for_cube(pi_warm_cube, lsmfile)
    pi_cold_cube = mask_ocean_for_cube(pi_cold_cube, lsmfile)
   

    #fig1 = plt.figure(figsize=[8.0, 4.0], constrained_layout=True)
   
    fig1 = plt.figure(figsize=[8.0, 4.0])
    gs = gridspec.GridSpec(nrows=2, ncols=2)

    ax1 = fig1.add_subplot(gs[0,0], projection=ccrs.PlateCarree())
    plotcube(pi_landmean_cube, 'PI', 'mean', np.arange(-20, 20, 5), ax1)
    ax2 = fig1.add_subplot(gs[0,1], projection=ccrs.PlateCarree())
    plotcube(pi_landamp_cube, 'PI','amplitude', np.arange(20, 50, 5), ax2)
    ax3 = fig1.add_subplot(gs[1,0], projection=ccrs.PlateCarree())
    plotcube(plio_landmean_cube, 'mPWP', 'mean', np.arange(-20, 20, 5), ax3)
    ax4 = fig1.add_subplot(gs[1,1], projection=ccrs.PlateCarree()) 
    plotcube(plio_landamp_cube, 'mPWP', 'amplitude', np.arange(20, 50, 5), ax4)
  
   
    
    plt.tight_layout()
    #plt.subplots_adjust(bottom=0.1, top=0.9, hspace=None)
    
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/map_seas_cyc.eps')
   
    plt.savefig(fileout)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/map_seas_cyc.png')
   
    plt.savefig(fileout)
    plt.close()
  


    # anomalies
    #plotcube(plio_landamp_cube - pi_landamp_cube, 'mPWP-PI','amplitude',
    #         np.arange(-5, 5, 1))
    #plotcube(plio_landmean_cube - pi_landmean_cube, 'mPWP-PI','mean',
    #         np.arange(-5, 5, 1))

    # plot all pi locations which have an annual cycle within 1 deg of 
    # ann_cyc_req

    # note that I think this is a bit misleading.
    #  all it shows are that these temperatures can exist in the preindustrial
    plot_locations_bothmatch(pi_warm_cube, pi_cold_cube)

    plot_locations2(pi_warm_cube, pi_cold_cube, pi_precip_cube,
                    WARM_TEMP_REQ, 'warm_match')
    plot_locations2(pi_cold_cube, pi_warm_cube, pi_precip_cube,
                    COLD_TEMP_REQ, 'cold_match')
    plot_locations(pi_landamp_cube, pi_warm_cube, pi_cold_cube)
      

##########################################################
# main program

LINUX_WIN = 'l'
EXPTNAME = 'EOI400'
PI_EXPT = 'E280'

FILENAME= ('/nfs/hera1/earjcti/regridded100/' + 
                    'NearSurfaceTemperature_multimodelmean_month.nc')
PRECIPFILE  = ('/nfs/hera1/earjcti/regridded100/' + 
                    'TotalPrecipitation_multimodelmean.nc')

BCSTART = '/nfs/hera1/earjcti/regridded/PlioMIP2_Boundary_conds/'

ANN_CYC = {'Lake Baikal' : 17.0, 'Lake Baikal M': 46.0,
           'Lost Chicken Mine' : 14.0, 'Lost Chicken Mine M' : 35.0,
           'Lost Chicken Mine B' : 39.0}
#WARM_TEMP = {'Lake Baikal' : 22.4, 'Lost Chicken Mine' : 15.4, 
#             'James Bay Lowland' : 22.6}
WARM_TEMP = {'Lake Baikal' : 16.4, 'Lake Baikal M' : 20.0, 
             'Lost Chicken Mine' : 12.0, 'Lost Chicken Mine M': 15.4,
             'Lost Chicken Mine B' : 15.0}
COLD_TEMP = {'Lake Baikal' : -0.3, 'Lake Baikal M' : -15.0,
             'Lost Chicken Mine' : -2.0, 'Lost Chicken Mine M': -20,
             'Lost Chicken Mine B' : -24.0}

SITE = 'Lake Baikal M'  #'James Bay Lowland' #'Lost Chicken Mine' #'Lake Baikal' # Lost Chicken Mine B (beetle) Lake Baikal M (model

ANN_CYC_REQ = ANN_CYC.get(SITE)
WARM_TEMP_REQ = WARM_TEMP.get(SITE)
COLD_TEMP_REQ = COLD_TEMP.get(SITE)
main()
