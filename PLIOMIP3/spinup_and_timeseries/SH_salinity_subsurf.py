#!/usr/bin/env python2.7
#NAME
#
# This program will get the SH salinity for each year averaged over the 
# following
# 0-50m, 0-100m, 0-200m, 0-500m, 0-1000m and write them to a file
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



def get_avg(year, latstart):
    """
    gets average salinity for this year between latstart and -90
    """

    yearuse = str(year).zfill(9)
    filename=('/uolstore/Research/a/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    cube = iris.load_cube(filename,
                          'SALINITY (OCEAN)       (PSU-35)/1000')

    cube = iris.util.squeeze(cube)
  
    cube=cube * 1000. + 35.0
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    cube.coord('depth_1').guess_bounds()

    w_ro = iris.analysis.cartography.area_weights(cube) # this is in m2
    weights = np.copy(w_ro) 
    
    for j,lat in enumerate(cube.coord('latitude').points):
        if lat > latstart:
            weights[:,j,:]=0.0
    

    # get the weights for each of the ranges

    low_lims = [50.,100.,200.,500.,1000.]

    SH_sal_data_lims = []
    for lim in low_lims:
        print(lim)
        weights_use = np.copy(weights)
        for k,dep in enumerate(cube.coord('depth_1').points):
            if dep > lim:
                weights_use[k,:,:] = 0.0


        SH_sal_mean = cube.collapsed(['depth_1','latitude','longitude'],
                                 iris.analysis.MEAN,weights=weights_use)
        SH_sal_data_lims.append(SH_sal_mean.data)
    
  
    return (SH_sal_data_lims)

#####################################################################
def plotdrifts(salinity_0_50,salinity_0_100,salinity_0_200,salinity_0_500,salinity_0_1000,startyear,endyear,latstart):
    """
    plots the salinity
    write the timeseries of seaice area
    to a file
    """

    plt.subplot(1,1,1)
    plt.plot(salinity_0_50)
    plt.ylim(33.0,34.5)
    plt.title('Salinity ' + P3name.get(exptname,exptname))
    plt.ylabel('psu')
    plt.xlabel('year')
    if latstart < 0:
        latstartuse = str(latstart * -1.0) + 'S'
    else:
        latstartuse = str(latstart) + 'N'
    

    fileout=('/uolstore/Research/a/hera1/earjcti/um/'+exptname+'/spinup/salinity_subsurface_'+exptname+'_' + latstartuse + '-90S_' + str(int(startyear)) + '_'+ str(int(endyear)) +'.tex') 

    f= open(fileout,'w')
    f.write("year,    SH mean salinity:0-50,  0-100, 0-200, 0-500,0_1000 \n")
    for year in range(startyear,endyear):
        string = (str(year) + ','+ 
                  str(np.round(salinity_0_50[year-startyear],2)) + ',' +  
                  str(np.round(salinity_0_100[year-startyear],2)) + ','+
                  str(np.round(salinity_0_200[year-startyear],2)) + ','+
                  str(np.round(salinity_0_500[year-startyear],2)) + ','+
                  str(np.round(salinity_0_100[year-startyear],2)) + 
                  '\n')
        f.write(string)
    f.close()

    plt.show()
    plt.close()

   


#######################################################
def get_SHsalin(HadCM3,exptname,startyear,endyear,latstart):
    """
    reads in the salinity and extracts for the SH
    """

    # arrays for storing seaice
    salinity_SH_0_50  = np.zeros(endyear-startyear+1)
    salinity_SH_0_100  = np.zeros(endyear-startyear+1)
    salinity_SH_0_200  = np.zeros(endyear-startyear+1)
    salinity_SH_0_500  = np.zeros(endyear-startyear+1)
    salinity_SH_0_1000  = np.zeros(endyear-startyear+1)

    # obtain means for each year and store in the arrays
    
    for year in range(startyear,endyear+1):
        print(year)
        S_0_50,S_0_100,S_0_200,S_0_500,S_0_1000 = get_avg(year,latstart)
     
        salinity_SH_0_50[year-startyear]=S_0_50
        salinity_SH_0_100[year-startyear]=S_0_100
        salinity_SH_0_200[year-startyear]=S_0_200
        salinity_SH_0_500[year-startyear]=S_0_500
        salinity_SH_0_1000[year-startyear]=S_0_1000

    # plot and save
    plotdrifts(salinity_SH_0_50,salinity_SH_0_100,salinity_SH_0_200,
            salinity_SH_0_500,salinity_SH_0_1000,startyear,endyear,latstart)
  


################################
# main program

# annual mean
figureno=0

P3name = {'xpsie' : 'EP400',    'xpsig':'EP490', 'xpsid':'LP','xpsic':'PI',
          'xpsij':'LP490'}

latstart=-60.0   # will plot from latstart to -90.0
HadCM3='y'
exptname='xpsid'
startyear=12  # can't start before year 12 because we aren't outputting d18o
endyear=2999
plt.figure(figureno)
get_SHsalin(HadCM3,exptname,startyear,endyear,latstart)
figureno=figureno+1





sys.exit(0)

####

