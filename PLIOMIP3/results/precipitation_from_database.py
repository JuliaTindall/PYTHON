"""
#NAME
#    Summer_warming.py
#PURPOSE 
#
#  This program will show the precipitation change of each 
#  experiment for the annual mean, djf and jja
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




def get_season(jobid, startyear, endyear):
    """
    gets the average data fpr the field
    """  

    filename = ('/nfs/hera1/earjcti/um/' + jobid + '/database_averages/' + 
                jobid + '_Monthly_Average_#pd_' + FIELD + '_' + STARTYEAR + 
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
    # if precipitation conver to mm/day
    if FIELD == 'TotalPrecipitationRate':
        print('here')
        cube.data = cube.data * 24. * 60. * 60.
   
    ann_cube = cube.collapsed(['time'],iris.analysis.MEAN)

    # get the values for the seasons
   
    djf_cubelist = CubeList([])
    djf_cubelist.append(cube[0,:,:,:])
    djf_cubelist.append(cube[1,:,:,:])
    djf_cubelist.append(cube[11,:,:,:])
    djf_3cube = djf_cubelist.merge_cube()
    djf_cube = djf_3cube.collapsed(['time'],iris.analysis.MEAN)
   
    jja_cubelist = CubeList([])
    jja_cubelist.append(cube[5,:,:,:])
    jja_cubelist.append(cube[6,:,:,:])
    jja_cubelist.append(cube[7,:,:,:])
    jja_3cube = jja_cubelist.merge_cube()
    jja_cube = jja_3cube.collapsed(['time'],iris.analysis.MEAN)
   
    return djf_cube, jja_cube, ann_cube
                   
    
#=================================================================
# MAIN PROGRAM STARTS HERE



LINUX_WIN='l'
NYEARS = 100
SEASON = 'ann'

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPT = 'xqbwi'  # xsic PI,  xpsij-lp490  xpsik - lp560
CNTL = 'xqbwc'  # xpsic pi, xpsid lp400
STARTYEAR='3900'
ENDYEAR='4000'

FIELD = 'TotalPrecipitationRate'

cntl_cube_djf, cntl_cube_jja, cntl_cube_ann = get_season(CNTL,STARTYEAR,ENDYEAR)
expt_cube_djf, expt_cube_jja, expt_cube_ann = get_season(EXPT,STARTYEAR,ENDYEAR)

diff_cube_djf = expt_cube_djf - cntl_cube_djf
diff_cube_jja = expt_cube_jja - cntl_cube_jja
diff_cube_ann = expt_cube_ann - cntl_cube_ann
print('got cntl cube')

#boundaries = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
boundaries = np.arange(-2.0,2.2,0.2)
cmapa=customise_cmap2()
cmap_use=plt.cm.get_cmap(cmapa,len(boundaries))
#cmap_use.set_under('lightsteelblue')

# plot djf
diff_cube_djf.coord('longitude').guess_bounds()
print('j1',diff_cube_djf.data)
diff_cube_djf.coord('latitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(diff_cube_djf)
meandiff_djf = diff_cube_djf.collapsed(['longitude','latitude'],
                               iris.analysis.MEAN, weights=grid_areas)
diffchar_djf = str(np.around(meandiff_djf.data,2))

print('j2',diff_cube_djf.data)
#cs=iplt.pcolormesh(diff_cube_djf,cmap=cmap_use,
#                 norm=mp.colors.BoundaryNorm(boundaries, 
#                                             ncolors=len(boundaries)-1,
#                                             clip=False))
cs=iplt.contourf(diff_cube_djf,levels=boundaries,cmap=cmap_use,extend='both')
cbar=plt.colorbar(cs,orientation='horizontal',extend='both')
#cbar.set_ticks(boundaries)
cbar.set_label('mm/day')
titlename = EXPT + '-' +  CNTL + '. Years:' + str(STARTYEAR) + '-' + str(ENDYEAR) + '. DJF. Meandiff =' +  diffchar_djf
plt.title(titlename, fontsize=10)
plt.gca().coastlines()
plt.show()

print('about to write to file')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + FIELD + '.eps')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/djf_' + EXPT + '-' + CNTL + '_' + FIELD + '.png')
plt.close()


# plot jja
diff_cube_jja.coord('longitude').guess_bounds()
diff_cube_jja.coord('latitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(diff_cube_jja)
meandiff_jja = diff_cube_jja.collapsed(['longitude','latitude'],
                               iris.analysis.MEAN, weights=grid_areas)
diffchar_jja = str(np.around(meandiff_jja.data,2))

#cs=iplt.pcolormesh(diff_cube_jja,cmap=cmap_use,
#                 norm=mp.colors.BoundaryNorm(boundaries, 
#                                             ncolors=len(boundaries)-1,
#                                             clip=False))
cs=iplt.contourf(diff_cube_jja,levels=boundaries,cmap=cmap_use,extend='both')

cbar=plt.colorbar(cs,orientation='horizontal',extend='both')
#cbar.set_ticks(boundaries)
cbar.set_label('degC')
titlename = EXPT + '-' +  CNTL + '. Years:' + str(STARTYEAR) + '-' + str(ENDYEAR) + '.  JJA. Meandiff =' +  diffchar_jja
plt.title(titlename, fontsize=10)
plt.gca().coastlines()
      
print('about to write to file')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja_' + EXPT + '-' + CNTL + '_' + FIELD + '.eps')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja_' + EXPT + '-' + CNTL + '_' + FIELD + '.png')
plt.close()
