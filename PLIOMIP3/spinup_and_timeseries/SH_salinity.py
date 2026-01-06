#!/usr/bin/env python2.7
#NAME
#
# This program will plots the Sea Surface Salinity south of a given latitude
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
    filename=('/nfs/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    cube_alllevs = iris.load_cube(filename,
                                  'SALINITY (OCEAN)       (PSU-35)/1000')
    cube = cube_alllevs[0,0,:,:]  # this is surface salinity
    cube=cube * 1000. + 35.0
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    weights = iris.analysis.cartography.area_weights(cube) # this is in m2
  
    for j,lat in enumerate(cube.coord('latitude').points):
        if lat > latstart:
            weights[j,:]=0.0
      
    SH_sal_mean = cube.collapsed(['latitude','longitude'],
                                 iris.analysis.MEAN,weights=weights)

  
    return (SH_sal_mean.data)

#####################################################################
def plotdrifts(salinity,startyear,endyear,latstart):
    """
    plots the timeseries of seaice area
    file
    """

    plt.subplot(1,1,1)
    plt.plot(salinity)
    plt.ylim(33.0,34.5)
    plt.title('Salinity ' + P3name.get(exptname,exptname))
    plt.ylabel('psu')
    plt.xlabel('year')
    if latstart < 0:
        latstartuse = str(latstart * -1.0) + 'S'
    else:
        latstartuse = str(latstart) + 'N'


    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/salinity_'+exptname+'_' + latstartuse + '-90S_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.eps') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/salinity_'+exptname+'_' + latstartuse + '-90S_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.png') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    
    plt.close()

    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/salinity_'+exptname+'_' + latstartuse + '-90S_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.tex') 

    f= open(fileout,'w')
    f.write("year,    SH mean salinity: top level \n")
    for year in range(startyear,endyear):
        string = (np.str(year) + ','+ 
                  np.str(np.round(salinity[year-startyear],2)) + '\n')
        f.write(string)
    f.close()



   


#######################################################
def get_SHsalin(HadCM3,exptname,startyear,endyear,latstart):
    """
    reads in the salinity and extracts for the SH
    """

    # arrays for storing seaice
    salinity_SH  = np.zeros(endyear-startyear+1)

    # obtain means for each year and store in the arrays
    
    for year in range(startyear,endyear+1):
        print(year)
        salin = get_avg(year,latstart)
        salinity_SH[year-startyear]=salin

    # plot and save
    plotdrifts(salinity_SH,startyear,endyear,latstart)
  


################################
# main program

# annual mean
figureno=0

P3name = {'xpsie' : 'EP400',    'xpsig':'EP490', 'xpsid':'LP','xpsic':'PI',
          'xpsij':'LP490'}

latstart=-60.0   # will plot from latstart to -90.0
HadCM3='y'
exptname='xpsij'
startyear=1991  # can't start before year 12 because we aren't outputting d18o
endyear=2999
plt.figure(figureno)
get_SHsalin(HadCM3,exptname,startyear,endyear,latstart)
figureno=figureno+1





sys.exit(0)

####

