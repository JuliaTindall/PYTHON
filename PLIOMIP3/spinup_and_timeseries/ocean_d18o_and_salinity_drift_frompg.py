#!/usr/bin/env python2.7
#NAME
#    ocean_d18o_and_salinity_drift_frompg
#PURPOSE
#    This program is loosely based on 
#  IDLPRGS/HadGEM/analysis/temperature_timeseries_frompg.pro
#
#  It should do a timeseries of the average ocean d18o and salinity
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


def get_avg(year,weights,incl18o):
    """
    gets the average d18o and salinity for this year
    year = year
    weights = weights for averaging.  Will be calculated here if they are unset
    """
 
    yearuse = str(year).zfill(9)
    filename=('/uolstore/home/users/earjcti/hera1/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    cubelist = iris.load(filename)
    salratiocube = cubelist.extract('salinity')[0]

    if incl18o=='y':
        ratio18ocube = cubelist.extract('otracer1')[0]
    
    if np.ndim(weights) == 1: # this has not yet been setup 
        # get weights for calculating mean - scale by depth
        salratiocube.coord('latitude').guess_bounds()
        salratiocube.coord('longitude').guess_bounds()
        weights = iris.analysis.cartography.area_weights(salratiocube)
        new_weights = np.zeros(np.shape(weights))
        depths = salratiocube.coord('depth_1').points
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
    meansalratio_cube = salratiocube.collapsed(
        ['longitude','latitude','depth_1'],
        iris.analysis.MEAN,
        weights=new_weights)

    if incl18o=='y':
        mean18oratio_cube = ratio18ocube.collapsed(
            ['longitude','latitude','depth_1'],
            iris.analysis.MEAN,
            weights=new_weights)
    else:
        mean18oratio_cube=0

    # get mean for levels
    # get_mean_level accepts cube, weights, startlevel, endlevel)
    meansalratio_lev1_10_cube = get_mean_levels(salratiocube, new_weights, 0, 9)
    meansalratio_lev1_5_cube = get_mean_levels(salratiocube, new_weights, 0, 4)
    meansalratio_lev1_cube = get_mean_levels(salratiocube, new_weights, 0, 0)

    if incl18o == 'y':
        mean18oratio_lev1_5_cube = get_mean_levels(ratio18ocube,
                                                   new_weights, 0, 4)
        mean18oratio_lev1_10_cube = get_mean_levels(ratio18ocube, 
                                                    new_weights, 0, 9)
        mean18oratio_lev1_cube = get_mean_levels(ratio18ocube, 
                                                 new_weights, 0, 0)
    else:
        mean18oratio_lev1_5_cube=0
        mean18oratio_lev1_10_cube=0
        mean18oratio_lev1_cube =0
    

    return (new_weights, meansalratio_cube, mean18oratio_cube,
            meansalratio_lev1_10_cube, mean18oratio_lev1_10_cube,
            meansalratio_lev1_5_cube, mean18oratio_lev1_5_cube,
            meansalratio_lev1_cube, mean18oratio_lev1_cube)

#####################################################################
def plotdrifts(ratio18o, ratiosal, levels_plotted, incl18o):
    """
    plots the drifts for salinity and d18o and saves to a
    file
    """

    if incl18o=='y':
        plt.subplot(2,1,1)
    else:
        plt.subplot(1,1,1)
    plt.plot((ratiosal * 1000) + 35.0)
    plt.title('ocean salinity - ' + levels_plotted)
    plt.ylabel('psu')
    plt.xlabel('year')

    if incl18o=='y':
        plt.subplot(2,1,2)
        plt.plot((ratio18o - 2005.2E-6) / 2005.2E-9)
        plt.title('ocean d18o - ' + levels_plotted)
        plt.ylabel('permille')
        plt.xlabel('year')
   
    plt.tight_layout()

    fileout=('/uolstore/home/users/earjcti/PYTHON/PLOTS/PLIOMIP3/assess_spinup/ocn_d18o_salratio/drift_'+exptname+'_' + levels_plotted + '.eps') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


#######################################################
def getdrift(HadCM3,exptname,startyear,endyear,incl18o):
    """
    gets the ocean salratioerature drift for exeriment = exptname between
    years startyear and endyear
    """

    # arrays for storing mean salratioeratures
    globalmean_18oratio = np.zeros(endyear-startyear+1)
    mean_18oratio_lev1  = np.zeros(endyear-startyear+1)
    mean_18oratio_lev1_5  = np.zeros(endyear-startyear+1)
    mean_18oratio_lev1_10  = np.zeros(endyear-startyear+1)
    globalmean_salratio = np.zeros(endyear - startyear+1)
    mean_salratio_lev1  = np.zeros(endyear-startyear+1)
    mean_salratio_lev1_5  = np.zeros(endyear-startyear+1)
    mean_salratio_lev1_10  = np.zeros(endyear-startyear+1)

    # obtain means for each year and store in the arrays
    weights = [0.0]
    for year in range(startyear,endyear+1):
        print(year)
        (weights, meansalratio_cube, mean18oratio_cube,
        meansalratio_lev1_10_cube, mean18oratio_lev1_10_cube,
        meansalratio_lev1_5_cube, mean18oratio_lev1_5_cube,
        meansalratio_lev1_cube, mean18oratio_lev1_cube) = get_avg(year,
                                                                  weights,
                                                                  incl18o)
       
        globalmean_salratio[year-startyear]  = meansalratio_cube.data
        mean_salratio_lev1[year-startyear]  = meansalratio_lev1_cube.data
        mean_salratio_lev1_5[year-startyear]  = meansalratio_lev1_5_cube.data
        mean_salratio_lev1_10[year-startyear]  = meansalratio_lev1_10_cube.data

        if incl18o=='y':
            globalmean_18oratio[year-startyear] = mean18oratio_cube.data
            mean_18oratio_lev1[year-startyear]  = mean18oratio_lev1_cube.data
            mean_18oratio_lev1_10[year-startyear]  = mean18oratio_lev1_10_cube.data
            mean_18oratio_lev1_5[year-startyear]  = mean18oratio_lev1_5_cube.data


    # plot and save
    plotdrifts(globalmean_18oratio, globalmean_salratio, 'alllevs',incl18o)
    plotdrifts(mean_18oratio_lev1, mean_salratio_lev1, 'lev1',incl18o)
    plotdrifts(mean_18oratio_lev1_5, mean_salratio_lev1_5, 'lev1_5',incl18o)
    plotdrifts(mean_18oratio_lev1_10, mean_salratio_lev1_10, 'lev1_10',incl18o)
  


################################
# main program

# annual mean
figureno=0

incl_18o='n'
HadCM3='y'
exptname='xpsie'
startyear=12  # can't start before year 12 because we aren't outputting d18o
endyear=2999
plt.figure(figureno)
getdrift(HadCM3,exptname,startyear,endyear,incl_18o)
figureno=figureno+1





sys.exit(0)

####

