"""
#NAME
#    salinity_depth.py
#PURPOSE 
#
#  This program will plot the salinity at a given depth
#  it will also provide a difference from a control
"""

# Import necessary libraries

import os
import numpy as np
import math
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.quickplot as qplt
import iris.plot as iplt
from netCDF4 import Dataset, MFDataset
import sys
#from netCDF4 import Dataset, MFDataset

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



def get_avg(jobid, startyear):
    """
    gets the average data fpr the field
    """  
    
    allsalinity = np.zeros((144,288))
    count=0.
    field = 'SALINITY (OCEAN)       (PSU-35)/1000'

    # get template map
    cube = iris.load_cube('/nfs/hera1/earjcti/um/xqbwc/pg/' + 
                          'xqbwco#pg' + str(3900).zfill(9) + 'c1+.nc',
                          field)

    cube = iris.util.squeeze(cube)[LEVEL,:,:]
      
    for year in range(startyear, startyear+NYEARS):
        
        files = ('/nfs/hera1/earjcti/um/' + jobid + '/pg/' + 
                 jobid + 'o#pg' + str(year).zfill(9) + 'c1+.nc')
        print(files)
        f=MFDataset(files)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        depths = f.variables['depth_1'][:]
        sal=np.squeeze(f.variables['salinity'][:])
        sallev = sal[LEVEL,:,:]
        levreq=depths[LEVEL]

        allsalinity = allsalinity + sallev
        count=count+1.

    avgsaldata = ((allsalinity / count) * 1000.) + 35.0
       
    avgsalcube=cube.copy(data=avgsaldata)

    
    # plot average
    vals = np.arange(30.0, 38.5,0.5)

    qplt.contourf(avgsalcube,levels=vals,extend='both')
    titlename = (jobid + '. Years:' + str(startyear) + '-' + str(startyear + NYEARS) +  ' depth=' + np.str(levreq))
    plt.title(titlename, fontsize=10)
    plt.gca().coastlines()
      
    plt.savefig('/nfs/hera1/earjcti/um/' + jobid +  '/avgplots/salinity/salinity_' + jobid  + '_' + np.str(levreq) + '.eps')
    plt.savefig('/nfs/hera1/earjcti/um/' + jobid +  '/avgplots/salinity/salinity_' + jobid  + '_' + np.str(levreq) + '.png')
    
    plt.close()



    return avgsalcube,levreq
                   
    
#=================================================================
# MAIN PROGRAM STARTS HERE



LINUX_WIN='l'
NYEARS = 100
SEASON = 'ann'

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPT = 'xqbwe'  # xpsic PI,  xpsij-lp490  xpsik - lp560
EXPT_STARTYEAR = 3900
#EXPT = 'Eoi400_ARC4_2450-2499'

# data from good experiment
CNTL = 'xqbwc'  # xpsic pi, xpsid lp400
CNTL_STARTYEAR = 3900

LEVEL=0


expt_cube, levreq = get_avg(EXPT,EXPT_STARTYEAR)
print('got expt_cube')
cntl_cube, levreq = get_avg(CNTL,CNTL_STARTYEAR)
print('got cntl cube')
diff_cube = expt_cube - cntl_cube

print(np.shape(diff_cube))

# get means
diff_cube.coord('longitude').guess_bounds()
diff_cube.coord('latitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(diff_cube)
meandiff = diff_cube.collapsed(['longitude','latitude'],
                               iris.analysis.MEAN, weights=grid_areas)
diffchar = str(np.around(meandiff.data,2))
   
print('about to plot')
vals = np.arange(-2.0, 2.25,0.25)

qplt.contourf(diff_cube,levels=vals,extend='both',cmap='RdBu_r')
titlename = EXPT + '-' +  CNTL + '. Years:' + str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR + NYEARS) + '.Meandiff =' +  diffchar + ' depth=' + np.str(levreq)
plt.title(titlename, fontsize=10)
plt.gca().coastlines()
      
print('about to write to file')
print(levreq)
#plt.show()
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/salinity/salinity_' + EXPT + '-' + CNTL + '_' + np.str(levreq) + '.eps')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/salinity/salinity_' + EXPT + '-' + CNTL + '_' + np.str(levreq) + '.png')
plt.close()
