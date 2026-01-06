# -*- coding: utf-8 -*-
"""
This program will do a DMC plot.  But it will be latitude vs deltSST

"""
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import iris
import sys

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def get_data():
    """
    this function willl open the file containing the data and will return 
    arrays containing:
        1. site_longitude
        2. site_latitudes
        3. T anomaly from (NOAA-ERSSTv5)
        4. standard deviation of the data
        5. number of sites
    """
    dfs = pd.read_excel(DATAFILE)
    dfs_subset = dfs[["Latitude (¡N)", "Longitude (¡E)", "NOAA_anom", "Standard dev.", "N"]]
    
    lats = dfs_subset.iloc[:,0]
    lons = dfs_subset.iloc[:,1]
    data_tanom = dfs_subset.iloc[:,2]
    data_stdev = dfs_subset.iloc[:,3]
    Npoints = dfs_subset.iloc[:,4]
    
    
    return lats, lons, data_tanom, data_stdev, Npoints
    
def plot_data(lats, Tanom, stdev):
    """
    plots the data and the errorbars
    """
    
    stdevplot = np.zeros(len(stdev))
    for i in range(0, len(stdev)):
       numeric = is_number(stdev[i])
       if numeric:
            stdevplot[i] = stdev[i]
    
    print(stdevplot)

    plt.errorbar(lats, Tanom, yerr=stdevplot, fmt='o')
    

def get_multimodel_mean(fieldname):
    """
    gets the multimodel mean and calculates a zonal average
    """
    mmm_cube = iris.load_cube(MULTIMODELMEAN, fieldname)
    zm_cube = mmm_cube.collapsed(['longitude'], iris.analysis.MEAN)
    zm_cube_max = mmm_cube.collapsed(['longitude'], iris.analysis.MAX)
    zm_cube_min = mmm_cube.collapsed(['longitude'], iris.analysis.MIN)
    
    return zm_cube, zm_cube_max, zm_cube_min
    
def plot_zm(cubemean, cubemax, cubemin, max_cubemax, min_cubemin):
    """
    latitudinal plot + zonal range of multimodel mean
    cubemean, cubemax, cubemin are the mean and the range from the multimodel mean
    max_cubemax is the maximum longitude from the model with the maximum difference
    min_cubemin is the minimum longitude from the model with the minimum difference
    """
    
    lats = cubemean.coord('latitude').points
    datamean = cubemean.data
    datamax = cubemax.data
    datamin = cubemin.data
    
    fig, ax = plt.subplots() 
    ax.plot(lats, datamean)
    ax.fill_between(lats, min_cubemin.data, max_cubemax.data, alpha=0.4)
    ax.fill_between(lats, datamin, datamax, alpha=0.4)
    ax.set_ylim(-5.0,20.0)
   
    ax.set_xlabel('latitude')
    ax.set_ylabel('SST anomaly')
    


def main():
    """
    1. get the data
    2. get multimodel mean SST 
    3. plot the data on the same file as the MMM
    """

    lats, lons, data_tanom, data_stdev, Npoints = get_data()
    #
    zonal_mean_cube, zonal_max_cube, zonal_min_cube = get_multimodel_mean('SSTmean_anomaly')
    max_zonal_mean_cube, max_zonal_max_cube, max_zonal_min_cube = get_multimodel_mean('SSTmax_anomaly')
    min_zonal_mean_cube, min_zonal_max_cube, min_zonal_min_cube = get_multimodel_mean('SSTmin_anomaly')
   
   
    plot_zm(zonal_mean_cube, zonal_max_cube, 
            zonal_min_cube, max_zonal_max_cube, min_zonal_min_cube)
    plot_data(lats, data_tanom, data_stdev)
    
    
    plt.savefig(OUTNAME + '.eps')
    plt.savefig(OUTNAME + '.pdf')
    plt.close()

DATAFILE = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/cs_mp_sst_data_30k_plusNOAA.xlsx'
MULTIMODELMEAN = '/nfs/hera1/earjcti/regridded/SST_multimodelmean.nc'
OUTNAME = '/nfs/hera1/earjcti/regridded/allplots/SST/dmc_by_latitude'

main()
