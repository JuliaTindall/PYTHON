#!/usr/bin/env python2.7
#NAME
#
# This program will plot a timeseries of the NH temperature trend and the
# SH temperature trend
#


import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess

def format_lat(lat):
    direction = "N" if lat >= 0 else "S"
    return f"{abs(lat)}{direction}"

def get_avg(year,latbands):
    """
    gets average SST for the bands 30-30, 30-60N, 60-90N, 30-90N
                                          30-60S, 60-90S, 30-90S
    """
 
    yearuse = str(year).zfill(9)
    filename=('/nfs/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    cube_alllevs = iris.load_cube(filename,'TEMPERATURE (OCEAN)  DEG.C')
    cube = cube_alllevs[0,0,:,:]  # this is surface temperature
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    weights = iris.analysis.cartography.area_weights(cube) # this is in m2
  
    meantemperature=[]

    for latrange in latbands:
        latmin, latmax=latrange[0], latrange[1]
        weights_range = np.copy(weights)
        for j,lat in enumerate(cube.coord('latitude').points):
            if lat > latmax or lat < latmin:
                weights_range[j,:]=0.0
  
        meancubes = cube.collapsed(['latitude','longitude'],
                                 iris.analysis.MEAN,weights=weights_range)
        #print(meancubes)
        #print(meancubes.data)
        meantemperature.append(meancubes.data)
   
    npmeantemp = np.array(meantemperature)
  
    return (npmeantemp)

#####################################################################
def plotdrifts(temps,startyear,endyear,latbands):
    """
    plots the timeseries of temperature for each bad
    """

    plt.subplot(1,1,1)
    for i,lats in enumerate(latbands):
        
        plt.plot(temps[i,:])
    #plt.ylim(33.0,34.5)
    plt.title('temperature ' + P3name.get(exptname))
    plt.ylabel('degC')
    plt.xlabel('year')

   
    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/SST_bands__'+exptname+'_'  + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.eps') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/SST_bands_'+exptname+'_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.png') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    
    plt.close()

    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/SST_bands_'+exptname+'_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.tex') 

    f= open(fileout,'w')
    formatted_lats = ",".join(f"{format_lat(x)}-{format_lat(y)}" for x, y in latbands)
    f.write("year,"+ formatted_lats + " \n")

    for year in range(startyear,endyear):
        values = temps[:,year-startyear]
        round_vals = [f"{val:.2f}" for val in values] 
        line = f"{year}," + ",".join(round_vals)  + '\n'
        f.write(line)
    f.close()



   


#######################################################
def get_temperatures(HadCM3,exptname,startyear,endyear,latbands):
    """
    reads in the temperature and extracts for latitude bands
    """

    # arrays for storing NH, SH and tropical SST
    temps  = np.zeros((len(latbands),endyear-startyear+1))
   
    # obtain means for each year and store in the arrays
    
    for year in range(startyear,endyear+1):
        print(year)
        meantemps = get_avg(year,latbands)
        temps[:,year-startyear]=meantemps

    # plot and save
    plotdrifts(temps,startyear,endyear,latbands)
  


################################
# main program

# annual mean
figureno=0

P3name = {'xpsie' : 'EP400',    'xpsig':'EP490', 'xpsid':'LP','xpsic':'PI',
          'xpsij':'LP490'}

# bands that we calculate the temperatures over
latbands=[[-30,30], [30,60],[60,90],[30,90],[-60,-30],[-90,-60],[-90,-30]]
  


HadCM3='y'
exptname='xpsij'
startyear=1991  # can't start before year 12 because we aren't outputting d18o
endyear=2999
plt.figure(figureno)
get_temperatures(HadCM3,exptname,startyear,endyear,latbands)
figureno=figureno+1





sys.exit(0)

####

