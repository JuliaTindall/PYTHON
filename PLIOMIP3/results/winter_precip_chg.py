"""
#NAME
#    SAT_warming.py
#PURPOSE 
#
#  This program will show the precipitation change in the last ??? years of each 
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
import cartopy.feature as cfeature
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


def make_cmap(colors, position=None, bit=False):
    '''
    I didn't write this I found it on the web.
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mp.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap



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
    cube = iris.load_cube('/nfs/hera1/earjcti/um/xqbwd/pcpd/' + 
                          'xqbwda#pd' + str(startyear).zfill(9) + 'dc+.nc',
                          longfield.get(field))

    cube = iris.util.squeeze(cube)
      
    for year in range(startyear, startyear+NYEARS):
        dccube = iris.load_cube('/nfs/hera1/earjcti/um/' + jobid + '/pcpd/' + 
                 jobid + 'a#pd' + str(year-1).zfill(9) + 'dc+.nc',field)
        jacube = iris.load_cube('/nfs/hera1/earjcti/um/' + jobid + '/pcpd/' + 
                 jobid + 'a#pd' + str(year).zfill(9) + 'ja+.nc',field)
        fbcube = iris.load_cube('/nfs/hera1/earjcti/um/' + jobid + '/pcpd/' + 
                 jobid + 'a#pd' + str(year).zfill(9) + 'fb+.nc',field)

        tempavg = np.squeeze((dccube.data + jacube.data + fbcube.data) / 3.0)

        alltemps = alltemps + tempavg
        count=count+1.

     
    if field == 'precip': # convert to mm/day'
        alltemps = alltemps * 60. * 60. *24.
       
    avgtempcube=cube.copy(data=alltemps / count)
    avgtempcube.attributes = None
    avgtempcube.units = 'mm/day'
    return avgtempcube
                   
    
#=================================================================
# MAIN PROGRAM STARTS HERE



LINUX_WIN='l'
NYEARS = 100
#NYEARS=5
SEASON = 'djf'

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPTS = ['xqbwd']  # xpsic PI,  xpsij-lp490  xpsik - lp560
#EXPTS=['xqbwd']
EXPT_STARTYEAR = 3900
#EXPT = 'Eoi400_ARC4_2450-2499'

# data from good experiment
CNTL = 'xqbwc'  # xpsic pi, xpsid lp400
CNTL_STARTYEAR = EXPT_STARTYEAR

#FIELDS  = ['temp1.5','precip','cloud_cover','mslp',
#          'seaiceconc','oceansurftemp',
#          'surfsalinity','MLD', 'evapsea']

#FIELDS  = ['surfsalinity','MLD', 'evapsea']

#FIELDS_STEVE = ['temp1.5','precipmm','cloud_cover','mslp_hPa',
#          'icefrac','oceansurftempK',
#          'surfsalinity','MLD', 'evapsea']

#FIELDS = ['AMOC']
field = 'precip'

cntl_cube = get_avg(CNTL,CNTL_STARTYEAR)
print('got cntl cube')
 

for EXPT in EXPTS:
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
    vals=np.arange(-2.0,2.2,0.2)
#    if CNTL == 'xpsid' or CNTL == 'xpsig' or CNTL == 'xqbwd':
#        vals = np.arange(-4.0, 4.5,0.5)
#    if EXPT == 'xpsir':
#        vals = np.arange(-1.0, 1.1,0.1)
#    if CNTL == 'xpsic' or CNTL=='xqbwc'
#        vals = np.arange(-8.0, 9.0, 1.0)
        
    cmapname = 'RdBu_r'
    if field == 'precip':
        cmapname = customise_cmap2()

    print(vals)
    qplt.contourf(diff_cube,levels=vals,extend='both',cmap=cmapname)
    titlename = 'precip: ' + EXPT + '-' +  CNTL + '. Years:' + str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR + NYEARS) + '.Meandiff =' +  diffchar
    plt.title(titlename, fontsize=10)
    plt.gca().coastlines()
  
      
    print('about to write to file')
    plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + field + str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR+NYEARS) +  '.eps')
    plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + field+ str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR+NYEARS) + '.png')
    plt.close()
    
    if field == 'precip': # also plot percentage change
        ratiocube = expt_cube/cntl_cube
        ratiocube.units = None
        vals = np.arange(0.85,1.65,0.1)
        qplt.contourf(ratiocube,levels=vals,extend='both')
        titlename = EXPT + '/' +  CNTL + '. Years:' + str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR + NYEARS)
        plt.title(titlename, fontsize=10)
        plt.gca().coastlines()
      
        print('about to write to file')
        plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + field+ str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR+NYEARS) + '_ratio.eps')
        plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + field+ str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR+NYEARS) + '_ratio.png')

        # also plot north america
        plt.close()
        vals=np.arange(-0.5, 0.6, 0.1)
        print(vals)
        ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree(central_longitude=0.0))
        ax.set_extent([210., 300.,10,70], crs=ccrs.PlateCarree())
        #axplot=qplt.contourf(diff_cube,levels=vals,extend='both',
        #                     cmap=cmapname)
        axplot=qplt.contourf(diff_cube,cmap=cmapname,levels=vals,
                             extend='both')
        ax.add_feature(cfeature.STATES)
        #lons=diff_cube.coord('longitude').points
        #lats=diff_cube.coord('latitude').points
        #lons,lats=np.meshgrid(lons,lats)
             
        #ax.contourf(lons, lats,
        #             diff_cube.data,
        #             transform=ccrs.PlateCarree(),
        #             levels=vals,extend='both',cmap=cmapname)
        titlename = EXPT + '-' +  CNTL + '. Years:' + str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR + NYEARS) 
        plt.title(titlename, fontsize=10)
        plt.gca().coastlines()
        plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + field+ str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR+NYEARS) + '_NAmerica.eps')
        plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + field+ str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR+NYEARS) + '_NAmerica.png')
        plt.close()

  
