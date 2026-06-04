"""
#NAME
#    Summer_warming.py
#PURPOSE 
#
#  This program will show the SAT warming in the last 30 years of each 
#  experiment for the months JJA
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
from iris.cube import CubeList
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import LinearSegmentedColormap
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
    #bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    #if bit:
    #    for i in range(len(colors)):
    #        colors[i] = (bit_rgb[colors[i][0]],
    #                     bit_rgb[colors[i][1]],
    #                     bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mp.colors.LinearSegmentedColormap('my_colormap',cdict,len(colors))
    return cmap



def customise_cmap2(boundaries,mean):
    """
    as customise_cmap but 19 colors only + 2 white in middle added by Julia
    """
    ncols = len(boundaries)
    boundaries=boundaries+mean
    boundaries = np.round(boundaries,1)

    rdgy = mp.colormaps['RdGy_r'].resampled(ncols)
    colors = rdgy(range(ncols))

    for i,color in enumerate(colors):
        if boundaries[i]  < 0.0:
            #color[0] = 0.0
            #color[1]=0.0
            color[2]=1.0
        print(i,color,boundaries[i]+mean)

    #sys.exit(0)

    #colors = [(84, 48, 5), (113, 70, 16), (143, 93, 27), (173, 115, 38),
    #          (195, 137, 60), (206, 160, 97), (216, 182, 135),
    #          (227, 204, 173), (238, 226, 211), (248, 248, 247),
    #          (212, 230, 229), (176, 212, 209), (140, 194, 190),
    #          (103, 176, 170), (67, 158, 150), (44, 135, 127),
    #          (29, 110, 100), (14, 85, 74), (0, 60, 48)]
    my_cmap = make_cmap(colors, bit=True)

    return my_cmap,boundaries


    
def get_average(jobid, startyear, endyear):
    """
    gets the average data fpr the field
    """  

    filename = (filestart + jobid + '/database_averages/' + 
                jobid + '_Annual_Average_#pd_' + FIELD + '_' + STARTYEAR + 
                '_' + ENDYEAR + '.nc')

    longfield = {'Temperature' : 'TEMPERATURE AT 1.5M',
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

    cube = iris.load_cube(filename)
    cube = iris.util.squeeze(cube)
    
    return cube
                   
    
#=================================================================
# MAIN PROGRAM STARTS HERE

EXPTNAMES={'xqbwc':'PI', 'xqbwd':'LP','xqbwg':'EP',
           'xqbwr':'LP_alt','xqbwl':'PI$_{400}$','xqbwm':'PI$_{560}$',
           'xqbwi':'LP$_{280}$','xqbwj':'LP$_{490}$','xqbwk':'LP$_{560}$',
           'xqbwn':'LP$_{highNHorb}$','xqbwo':'LP$_{lowNHorb}$',
           'xqbwt':'PI$_{dynveg}$','xqbws':'LP$_{dynveg}$'}


LINUX_WIN='l'
NYEARS = 100
SEASON = 'ann'

filestart = '/uolstore/Research/a/hera1/earjcti/um/'
filestart = '/home/earjcti/um/'
# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPT = 'xqbwk'  # xsic PI,  xpsij-lp490  xpsik - lp560
CNTL = 'xqbwi'  # xpsic pi, xpsid lp400
STARTYEAR='3900'
ENDYEAR='4000'

FIELD = 'Temperature'

cntl_cube_ann = get_average(CNTL,STARTYEAR,ENDYEAR)
expt_cube_ann= get_average(EXPT,STARTYEAR,ENDYEAR)

diff_cube = expt_cube_ann - cntl_cube_ann

#boundaries = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
#boundaries = [0.0, 1.0, 2.0, 3.0,  4.0, 5.0, 6.0,7.0, 8.0,9.0, 10.0, 11.0, 12.0]
#cmap_use=plt.cm.get_cmap('Reds',len(boundaries))
#cmap_use.set_under('lightsteelblue')


try:
    diff_cube.coord('longitude').guess_bounds()
    diff_cube.coord('latitude').guess_bounds()
except:
    pass
grid_areas = iris.analysis.cartography.area_weights(diff_cube)
meandiff = diff_cube.collapsed(['longitude','latitude'],
                               iris.analysis.MEAN, weights=grid_areas)
diffchar = str(np.around(meandiff.data,2))

boundaries=np.arange(-6,6.1,.1)

cmap_use,boundaries = customise_cmap2(boundaries,meandiff.data) 


cs=iplt.pcolormesh(diff_cube,cmap=cmap_use,
                 norm=mp.colors.BoundaryNorm(boundaries, 
                                             ncolors=len(boundaries)-1,
                                             clip=False))
cbar=plt.colorbar(cs,orientation='horizontal',extend='both',format=lambda x, _: f"{x:.1f}")
cbar.set_ticks(boundaries)
print(len(boundaries))
#cbar.ax.set_yticklabels(['{:.0f}'.format(x) for x in boundaries])
cbar.locator = mticker.MaxNLocator(nbins=6)
cbar.update_ticks()

cbar.set_label('degC')
titlename = EXPTNAMES.get(EXPT) + '-' +  EXPTNAMES.get(CNTL) + '. Years:' + str(STARTYEAR) + '-' + str(ENDYEAR) + '. ANN. Meandiff =' +  diffchar
plt.title(titlename, fontsize=10)
plt.gca().coastlines()
plt.show()
print('about to write to file')
plt.savefig(filestart + EXPT +  '/avgplots/' + EXPT + '-' + CNTL + '_' + FIELD + '_shiftedscale..eps')
plt.savefig(filestart + EXPT +  '/avgplots/' + EXPT + '-' + CNTL + '_' + FIELD + '_shiftedscale.png')
plt.close()


