#!/usr/bin/env python2.7
#NAME
#
# This program will do a sea ice area timeseries for the NH and the SH
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
    filename=('/nfs/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    cube = iris.load_cube(filename,'AICE : ICE CONCENTRATION')
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

    # convert from m^2 to km^2
    NH_ice.data = NH_ice.data / (1000. * 1000.)
    SH_ice.data = SH_ice.data / (1000. * 1000.)


    return (NH_ice.data, SH_ice.data)

#####################################################################
def plotdrifts(seaicearea, hemisphere,startyear,endyear):
    """
    plots the timeseries of seaice area
    file
    """

    plt.subplot(1,1,1)
    plt.plot(seaicearea)
    plt.ylim(0.0,10000000)
  
    plt.title('seaicearea')
    plt.ylabel('km2')
    plt.xlabel('year')


    fileout=('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP3/assess_spinup/seaice/seaice_area_'+exptname+'_' + hemisphere + '_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.eps') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    fileout=('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP3/assess_spinup/seaice/seaice_area_'+exptname+'_' + hemisphere + '_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.png') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


    fileout=('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP3/assess_spinup/seaice/seaice_area_'+exptname+'_' + hemisphere + '_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.tex') 
    f= open(fileout,'w')
    f.write("year,    seaicearea " + hemisphere + "\n")
    for year in range(startyear,endyear):
        string = (np.str(year) + ','+ 
                  np.str(np.round(seaicearea[year-startyear],2)) + '\n')
        f.write(string)
    f.close()



#######################################################
def get_seaice(HadCM3,exptname,startyear,endyear):
    """
    reads in the sea ice area and stores over the NH and the SH
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
startyear=12  # can't start before year 12 because we aren't outputting d18o
endyear=2999
plt.figure(figureno)
get_seaice(HadCM3,exptname,startyear,endyear)
figureno=figureno+1





sys.exit(0)

####

