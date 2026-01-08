"""
#NAME
#    SAT_warming.py
#PURPOSE 
#
#  This program will show the SAT warming in the last nn years of each 
#  experiment
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

def customise_cmap2():
    """
    as customise_cmap but 19 colors only + 2 white in middle added by Julia
    """
    colors = [(84, 48, 5), (113, 70, 16), (143, 93, 27), (173, 115, 38),
              (195, 137, 60), (206, 160, 97), (216, 182, 135),
              (227, 204, 173), (238, 226, 211), (248, 248, 247),
              (212, 230, 229), (176, 212, 209), (140, 194, 190),
              (103, 176, 170), (67, 158, 150), (44, 135, 127),
              (29, 110, 100), (14, 85, 74), (0, 60, 48)]
    my_cmap = make_cmap(colors, bit=True)
    return my_cmap


def get_avg(jobid, startyear):
    """
    gets the average data fpr the field
    """  

    longfield = {'temp' : 'TEMPERATURE AT 1.5M',
                 'precip' : 'TOTAL PRECIPITATION RATE     KG/M2/S',
                 'precipmm' : 'TOTAL PRECIPITATION RATE MM/DAY',
                 'cloud_cover' : 'TOTAL CLOUD AMOUNT - RANDOM OVERLAP',
                 'mslp' : 'PRESSURE AT MEAN SEA LEVEL',
                 'mslp_hPa' : 'PRESSURE AT MEAN SEA LEVEL hPa',
                 'evapsea' : 'EVAPORATION FROM SEA (GBM)   KG/M2/S',
                 'seaiceconc' : 'AICE : ICE CONCENTRATION',
                 'icefrac' : 'SEA ICE FRACTION AFTER TIMESTEP',
                 'oceansurftemp' : 'OCN TOP-LEVEL TEMPERATURE          K',
                 'oceansurftempK' : 'OCN TOP-LEVEL TEMPERATURE K',
                 'surfsalinity': 'SALINITY (OCEAN)       (PSU-35)/1000',
                 'surfsalinitypsu': 'SALINITY (OCEAN) (PSU)',
                 'MLD' : 'MIXED LAYER DEPTH (OCEAN)          M',
                 'MLDm' : 'MIXED LAYER DEPTH (OCEAN) M',
                 'AMOC' : 'Meridional Overturning Stream Function (Atlantic)'   
                     }
    
    alltemps = np.zeros((73,96))
    count=0.

    # get template map
    cube = iris.load_cube('/nfs/hera1/earjcti/um/'+jobid+'/pcpd/' + 
                          jobid+'a#pd' + str(startyear).zfill(9) + 'dc+.nc',
                          longfield.get(field))

    cube = iris.util.squeeze(cube)
      
    for year in range(startyear, startyear+NYEARS):
        
        files = ('/nfs/hera1/earjcti/um/' + jobid + '/pcpd/' + 
                 jobid + 'a#pd' + str(year).zfill(9) + '??+.nc')
        print(files)
        f=MFDataset(files)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        temp=np.squeeze(f.variables[field][:])
        tempavg = np.mean(temp,axis=0)

        alltemps = alltemps + tempavg
        count=count+1.

       
    avgtempcube=cube.copy(data=alltemps / count)
    return avgtempcube
                   
    
#=================================================================
# MAIN PROGRAM STARTS HERE



LINUX_WIN='l'
NYEARS = 100
SEASON = 'ann'

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPTS = ['xqbwn','xqbwo']  # xpsic PI,  xpsij-lp490  xpsik - lp560
EXPT_STARTYEAR = 3900
#EXPT = 'Eoi400_ARC4_2450-2499'

# data from good experiment
CNTL = 'xqbwd'  # xpsic pi, xpsid lp400
CNTL_STARTYEAR = EXPT_STARTYEAR

field='temp'

cntl_cube = get_avg(CNTL,CNTL_STARTYEAR)
print('got cntl cube')
   
#for field in FIELDS:
for EXPT in EXPTS:
#  try:
    expt_cube = get_avg(EXPT,EXPT_STARTYEAR)
    print('got expt_cube')
    diff_cube = expt_cube - cntl_cube

    # get means
    diff_cube.coord('longitude').guess_bounds()
    diff_cube.coord('latitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(diff_cube)
    meandiff = diff_cube.collapsed(['longitude','latitude'],
                                   iris.analysis.MEAN, weights=grid_areas)
    diffchar = str(np.around(meandiff.data,2))
   
    print('about to plot')
    #mycmap = customise_cmap2()

    vals=np.arange(-4.0,4.4,0.5)
    vals = np.arange(-8.0, 9.0, 1.0)
   
    #if (CNTL == 'xpsid' or CNTL == 'xpsig' or CNTL=='xqbwd' or CNTL=='xqbwg'
    #    or CNTL == 'xqbwe' or CNTL=='xpsie'):
    #    vals = np.arange(-4.0, 4.5,0.5)
    #if EXPT == 'xpsir' or EXPT=='xqbwr':
    #vals = np.arange(-1.0, 1.1,0.1)
    #if CNTL == 'xpsic' or CNTL=='xqbwc'  or CNTL =='xqbwi':
    #    vals = np.arange(-8.0, 9.0, 1.0)
    qplt.contourf(diff_cube,levels=vals,extend='both',cmap='RdBu_r')
    titlename = EXPT + '-' +  CNTL + '. Years:' + str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR + NYEARS) + '.Meandiff =' +  diffchar
    plt.title(titlename, fontsize=10)
    plt.gca().coastlines()
      
    print('about to write to file')
    plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/' + EXPT + '-' + CNTL + '_' + field + '_' + str(EXPT_STARTYEAR) + '_' + str(EXPT_STARTYEAR+NYEARS) +  '.eps')
    plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/' + EXPT + '-' + CNTL + '_' + field + '_' + str(EXPT_STARTYEAR) + '_' + str(EXPT_STARTYEAR+NYEARS) +  '.png')
    plt.close()
#  except:
#    print('failed on',EXPT)
