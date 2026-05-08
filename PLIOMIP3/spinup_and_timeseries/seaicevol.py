#!/usr/bin/env python2.7
#NAME
#
# This program will do a sea ice volume timeseries for the NH and the SH
#
#  this is too HIGH DUE TO SOME GRIDBOXES WITH WAY TOO DEEP ICE.
#  SUGGEST USE SEAICE AREA INSTEAD
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



def get_avg(year):
    """
    gets total sea ice area for this year
    """
 
    yearuse = str(year).zfill(9)
    
    filename=('/uolstore/Research/a/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    try: 
        cube = iris.load_cube(filename,'HICE: MEAN ICE DEPTH OVER GRIDBOX  M')
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        weights = iris.analysis.cartography.area_weights(cube) # this is in m2
        SH_weights=np.copy(weights)
        NH_weights=np.copy(weights)
        
        for j,lat in enumerate(cube.coord('latitude').points):
            if lat < 0.0:
                NH_weights[0,0,j,:]=0.0
                if lat > 0.0:
                    SH_weights[0,0,j,:]=0.0
                    
        NH_ice_cube = iris.util.squeeze(cube.copy(data=cube.data * NH_weights))
        SH_ice_cube = iris.util.squeeze(cube.copy(data=cube.data * SH_weights))
            

        NH_ice = NH_ice_cube.collapsed(['latitude','longitude'],iris.analysis.SUM)
        SH_ice = SH_ice_cube.collapsed(['latitude','longitude'],iris.analysis.SUM)

        # convert from m^3 TO  KM 3
        NH_ice = NH_ice.data / (1000. * 1000. * 1000.)
        SH_ice = SH_ice.data / (1000. * 1000. * 1000.)
    except:
        print('not found',filename)
        if year > 3:
            sys.exit(0)
        NH_ice = np.nan
        SH_ice = np.nan


    return (NH_ice, SH_ice)

#####################################################################
def plotdrifts(seaicearea, hemisphere,startyear,endyear):
    """
    plots the timeseries of seaice VOLUME
    file
    """

    plt.subplot(1,1,1)
    plt.plot(seaicearea)
    plt.ylim(0.0,1000000)
  
    plt.title('seaiceVOLUME')
    plt.ylabel('km3')
    plt.xlabel('year')


    fileout=('/uolstore/home/users/earjcti/PYTHON/PLOTS/PLIOMIP3/assess_spinup/seaice/seaice_VOL_'+exptname+'_' + hemisphere + '_' + str(int(startyear)) + '_'+ str(int(endyear)) +'.eps') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    fileout=('/uolstore/home/users/earjcti/PYTHON/PLOTS/PLIOMIP3/assess_spinup/seaice/seaice_VOL_'+exptname+'_' + hemisphere + '_' + str(int(startyear)) + '_'+ str(int(endyear)) +'.png') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


    fileout=('/uolstore/home/users/earjcti/PYTHON/PLOTS/PLIOMIP3/assess_spinup/seaice/seaice_volume_'+exptname+'_' + hemisphere + '_' + str(int(startyear)) + '_'+ str(int(endyear)) +'.tex') 
    f= open(fileout,'w')
    f.write("year,    seaicearea " + hemisphere + "\n")
    for year in range(startyear,endyear):
        string = (str(year) + ','+ 
                  str(np.round(seaicearea[year-startyear],2)) + '\n')
        f.write(string)
    f.close()



#######################################################
def get_seaice(HadCM3,exptname,startyear,endyear):
    """
    reads in the sea ice volume and stores over the NH and the SH
    """

    # arrays for storing seaice
    seaicearea_NH  = np.zeros(endyear-startyear+1)
    seaicearea_SH  = np.zeros(endyear-startyear+1)

    # obtain means for each year and store in the arrays
    
    for year in range(startyear,endyear+1):
        print(year)
        (seaice_NH, seaice_SH) = get_avg(year)

        seaicearea_NH[year-startyear]=seaice_NH
        seaicearea_SH[year-startyear]=seaice_SH

    # plot and save
    plotdrifts(seaicearea_NH,'NH',startyear,endyear)
    plotdrifts(seaicearea_SH,'SH',startyear,endyear)
  


################################
# main program

# annual mean
figureno=0

incl_18o='n'
HadCM3='y'
exptname='xpsic'
startyear=0  # can't start before year 12 because we aren't outputting d18o
endyear=100
plt.figure(figureno)
get_seaice(HadCM3,exptname,startyear,endyear)
figureno=figureno+1





sys.exit(0)

####

