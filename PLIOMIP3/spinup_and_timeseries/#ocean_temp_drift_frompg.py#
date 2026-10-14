#!/usr/bin/env python2.7
#NAME
#    ocean_temp_drift_frompg
#PURPOSE
#    This program is loosely based on 
#  IDLPRGS/HadGEM/analysis/temperature_timeseries_frompg.pro
#
#  It should do a timeseries of the average ocean temperature
#
# search for 'main program' to find end of functions
# Julia 22/11/2016



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import iris
import cartopy.crs as ccrs
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess


def get_mean_levels(cube, weights, startlev, endlev):
    """
    get volume weighted mean between startlev and endlev
    """
    cubelevs = cube[:,startlev:endlev+1,:,:]
    weightlev = weights[:,startlev:endlev+1,: :]

    meancube = cubelevs.collapsed(['longitude','latitude','depth_1'],
                              iris.analysis.MEAN,
                              weights=weightlev)
    return meancube


def get_avg(year,weights):
    """
    gets the average temperature and potential temperature for this year
    year = year
    weights = weights for averaging.  Will be calculated here if they are unset
    """
 
    yearuse = str(year).zfill(9)
    filename=('/nfs/hera1/earjcti/um/'+exptname+'/netcdf/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    print(filename)
    cube = iris.load(filename)
    tempcube = cube.extract('insitu_T')[0]
    pottempcube = cube.extract('temp')[0]
    
    if np.ndim(weights) == 1: # this has not yet been setup 
        # get weights for calculating mean - scale by depth
        tempcube.coord('latitude').guess_bounds()
        tempcube.coord('longitude').guess_bounds()
        weights = iris.analysis.cartography.area_weights(tempcube)
        new_weights = np.zeros(np.shape(weights))
        depths = tempcube.coord('depth_1').points
        layer_depths = np.zeros(20)
        layer_depths[0] = depths[0] * 2.0
        layer_depths[19] = depths[19] - depths[18]
        for k in range(1,19):
            layer_depths[k]=((depths[k+1] - depths[k-1]) * 0.5) 
        for k in range(0,20):
            new_weights[:,k,:,:] = weights[:,k,:,:] * layer_depths[k]
    else:
        new_weights = weights

    # calculate mean and add to the 
    meantemp_cube = tempcube.collapsed(['longitude','latitude','depth_1'],
                                       iris.analysis.MEAN,
                                       weights=new_weights)
    meanpottemp_cube = pottempcube.collapsed(['longitude','latitude','depth_1'],
                                       iris.analysis.MEAN,
                                       weights=new_weights)

    # get mean for levels
    # get_mean_level accepts cube, weights, startlevel, endlevel)
    meantemp_lev1_10_cube = get_mean_levels(tempcube, new_weights, 0, 9)
    meanpottemp_lev1_10_cube = get_mean_levels(pottempcube, new_weights, 0, 9)
    meantemp_lev1_5_cube = get_mean_levels(tempcube, new_weights, 0, 4)
    meanpottemp_lev1_5_cube = get_mean_levels(pottempcube, new_weights, 0, 4)
    meantemp_lev1_cube = get_mean_levels(tempcube, new_weights, 0, 0)
    meanpottemp_lev1_cube = get_mean_levels(pottempcube, new_weights, 0, 0)
    

    return (new_weights, meantemp_cube, meanpottemp_cube,
            meantemp_lev1_10_cube, meanpottemp_lev1_10_cube,
            meantemp_lev1_5_cube, meanpottemp_lev1_5_cube,
            meantemp_lev1_cube, meanpottemp_lev1_cube)
#####################################################################
def plotdrifts(pottemp, temp, levels_plotted):
    """
    plots the drifts for potential temperature and temperature and saves to a
    file
    """

    plt.subplot(2,1,1)
    plt.plot(pottemp)
    plt.title('ocean potential temperature - ' + levels_plotted)
    plt.ylabel('degC')
    plt.xlabel('year')

    plt.subplot(2,1,2)
    plt.plot(temp)
    plt.title('ocean insitu temperature - ' + levels_plotted)
    plt.ylabel('degC')
    plt.xlabel('year')
   
    plt.tight_layout()


    fileout=('/nfs/hera1/earjcti/um/' + exptname + '/spinup/tdrift_'+exptname+'_' + levels_plotted + '.eps') 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


#######################################################
def getdrift(HadCM3,exptname,startyear,endyear):
    """
    gets the ocean temperature drift for exeriment = exptname between
    years startyear and endyear
    """

    # arrays for storing mean temperatures
    globalmean_pottemp = np.zeros(endyear-startyear+1)
    globalmean_temp = np.zeros(endyear - startyear+1)
    mean_pottemp_lev1  = np.zeros(endyear-startyear+1)
    mean_temp_lev1  = np.zeros(endyear-startyear+1)
    mean_pottemp_lev1_5  = np.zeros(endyear-startyear+1)
    mean_temp_lev1_5  = np.zeros(endyear-startyear+1)
    mean_pottemp_lev1_10  = np.zeros(endyear-startyear+1)
    mean_temp_lev1_10  = np.zeros(endyear-startyear+1)

    # obtain means for each year and store in the arrays
    weights = [0.0]
    for year in range(startyear,endyear+1):
        (weights, meantemp_cube, meanpottemp_cube,
        meantemp_lev1_10_cube, meanpottemp_lev1_10_cube,
        meantemp_lev1_5_cube, meanpottemp_lev1_5_cube,
        meantemp_lev1_cube, meanpottemp_lev1_cube) = get_avg(year,weights)
       
        globalmean_pottemp[year-startyear] = meanpottemp_cube.data
        globalmean_temp[year-startyear]  = meantemp_cube.data
        mean_pottemp_lev1[year-startyear]  = meanpottemp_lev1_cube.data
        mean_temp_lev1[year-startyear]  = meantemp_lev1_cube.data
        mean_pottemp_lev1_5[year-startyear]  = meanpottemp_lev1_5_cube.data
        mean_temp_lev1_5[year-startyear]  = meantemp_lev1_5_cube.data
        mean_pottemp_lev1_10[year-startyear]  = meanpottemp_lev1_10_cube.data
        mean_temp_lev1_10[year-startyear]  = meantemp_lev1_10_cube.data

    # plot and save
    plotdrifts(globalmean_pottemp, globalmean_temp, 'alllevs')
    plotdrifts(mean_pottemp_lev1, mean_temp_lev1, 'lev1')
    plotdrifts(mean_pottemp_lev1_5, mean_temp_lev1_5, 'lev1_5')
    plotdrifts(mean_pottemp_lev1_10, mean_temp_lev1_10, 'lev1_10')
  


################################
# main program

# annual mean
figureno=0


HadCM3='y'
exptname='xpsig'
startyear=2
endyear=3000
plt.figure(figureno)
getdrift(HadCM3,exptname,startyear,endyear)
figureno=figureno+1






sys.exit(0)

####

