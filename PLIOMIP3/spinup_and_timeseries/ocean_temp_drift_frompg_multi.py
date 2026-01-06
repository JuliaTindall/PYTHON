#!/usr/bin/env python2.7
#NAME
#    ocean_temp_drift_frompg
#PURPOSE
#    This program is loosely based on 
#  IDLPRGS/HadGEM/analysis/temperature_timeseries_frompg.pro
#
#  It should do a timeseries of the average ocean temperature.
#  It is designed for PlioMIP3 so will merge two experiments

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


def get_avg(year,weights,exptname):
    """
    gets the average temperature and potential temperature for this year
    year = year
    weights = weights for averaging.  Will be calculated here if they are unset
    """
 
    yearuse = str(year).zfill(9)
    filename=('/nfs/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
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
def plotdrifts(pottemp, temp, levels_plotted,startyear,endyear,exptname):
    """
    plots the drifts for potential temperature and temperature and saves to a
    file
    """
    years=np.arange(startyear,endyear,1)

    plt.subplot(2,1,1)
    plt.plot(years,pottemp)
    plt.title('ocean potential temperature - ' + levels_plotted)
    plt.ylabel('degC')
    plt.xlabel('year')

    plt.subplot(2,1,2)
    plt.plot(years,temp)
    plt.title('ocean insitu temperature - ' + levels_plotted)
    plt.ylabel('degC')
    plt.xlabel('year')
   
    plt.tight_layout()


    fileout=('/nfs/hera1/earjcti/um/' + exptname + '/spinup/tdrift_'+exptname+'_' + levels_plotted + '_'+str(startyear) + '_' + str(endyear) + '.eps') 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()

    filetxt=('/nfs/hera1/earjcti/um/' + exptname + '/spinup/tdrift_'+exptname+'_' + levels_plotted + '_'+str(startyear) + '_' + str(endyear) + '.txt')
    f=open(filetxt,'w')
    f.write('year,integrated temperature\n')
    for i in range(0,len(years)):
        f.write(str(years[i]) +',' + str(pottemp[i])+'\n')
    f.close()


#######################################################
def getdrift(ending):
    """
    gets the ocean temperature drift for exeriment = exptname between
    years startyear and endyear
    """
    if (ending == 'c' or ending=='d' or ending == 'e' or 
        ending == 'g' or ending=='a' or ending =='b' or 
        ending == 'h' or ending=='r'):
       yearstart_1=2
       #yearstart_1=2990
    else:
        yearstart_1=1991


    if ending == 'd':
        yearend_1=2891
    else:
        yearend_1=3001

    yearstart_2=yearend_1
    yearend_2=4001
    #yearend_2=3020
    nyears=yearend_2-yearstart_1



    # arrays for storing mean temperatures
    globalmean_pottemp = np.zeros(nyears)
    globalmean_temp = np.zeros(nyears)
    mean_pottemp_lev1  = np.zeros(nyears)
    mean_temp_lev1  = np.zeros(nyears)
    mean_pottemp_lev1_5  = np.zeros(nyears)
    mean_temp_lev1_5  = np.zeros(nyears)
    mean_pottemp_lev1_10  = np.zeros(nyears)
    mean_temp_lev1_10  = np.zeros(nyears)

    # obtain means for each year and store in the arrays
    weights = [0.0]
    for year in range(yearstart_1,yearend_1):
        expt='xpsi'+ending
        (weights, meantemp_cube, meanpottemp_cube,
        meantemp_lev1_10_cube, meanpottemp_lev1_10_cube,
        meantemp_lev1_5_cube, meanpottemp_lev1_5_cube,
        meantemp_lev1_cube, meanpottemp_lev1_cube) = get_avg(year,weights,expt)
       
        globalmean_pottemp[year-yearstart_1] = meanpottemp_cube.data
        globalmean_temp[year-yearstart_1]  = meantemp_cube.data
        mean_pottemp_lev1[year-yearstart_1]  = meanpottemp_lev1_cube.data
        mean_temp_lev1[year-yearstart_1]  = meantemp_lev1_cube.data
        mean_pottemp_lev1_5[year-yearstart_1]  = meanpottemp_lev1_5_cube.data
        mean_temp_lev1_5[year-yearstart_1]  = meantemp_lev1_5_cube.data
        mean_pottemp_lev1_10[year-yearstart_1]  = meanpottemp_lev1_10_cube.data
        mean_temp_lev1_10[year-yearstart_1]  = meantemp_lev1_10_cube.data

    for year in range(yearstart_2,yearend_2):
        expt='xqbw'+ending
        (weights, meantemp_cube, meanpottemp_cube,
        meantemp_lev1_10_cube, meanpottemp_lev1_10_cube,
        meantemp_lev1_5_cube, meanpottemp_lev1_5_cube,
        meantemp_lev1_cube, meanpottemp_lev1_cube) = get_avg(year,weights,expt)
       
        globalmean_pottemp[year-yearstart_1] = meanpottemp_cube.data
        globalmean_temp[year-yearstart_1]  = meantemp_cube.data
        mean_pottemp_lev1[year-yearstart_1]  = meanpottemp_lev1_cube.data
        mean_temp_lev1[year-yearstart_1]  = meantemp_lev1_cube.data
        mean_pottemp_lev1_5[year-yearstart_1]  = meanpottemp_lev1_5_cube.data
        mean_temp_lev1_5[year-yearstart_1]  = meantemp_lev1_5_cube.data
        mean_pottemp_lev1_10[year-yearstart_1]  = meanpottemp_lev1_10_cube.data
        mean_temp_lev1_10[year-yearstart_1]  = meantemp_lev1_10_cube.data


   
    # plot and save
    plotdrifts(globalmean_pottemp, globalmean_temp, 'alllevs',yearstart_1,
               yearend_2,expt)
    plotdrifts(mean_pottemp_lev1, mean_temp_lev1, 'lev1',yearstart_1,
               yearend_2,expt)
    plotdrifts(mean_pottemp_lev1_5, mean_temp_lev1_5, 'lev1_5',yearstart_1,
               yearend_2,expt)
    plotdrifts(mean_pottemp_lev1_10, mean_temp_lev1_10, 'lev1_10',yearstart_1,
               yearend_2,expt)
  


################################
# main program

# annual mean
figureno=0

#endings = ['d','e','g','i','j','k','l','m','n','o','p','q','r','s','t']
endings = ['t']
#endings=['e']
for ending in endings:
   # try:
        getdrift(ending)
   # except:
   #     print('failure on',ending)








sys.exit(0)

####

