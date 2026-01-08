"""
#NAME
#    Summer_minus_winter_temperatures.py
#PRPOSE 
#
#  This program will show the JJA-DJF temperatures both in absolute terms and 
#  relative to a control
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
   
    return djf_cube, jja_cube
                   
    
#=================================================================
# MAIN PROGRAM STARTS HERE



LINUX_WIN='l'
NYEARS = 100

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPT = 'xqbwo'  # xsic PI,  xpsij-lp490  xpsik - lp560
CNTL = 'xqbwc'  # xpsic pi, xpsid lp400
STARTYEAR='3900'
ENDYEAR='4000'

FIELD = 'Temperature'

cntl_cube_djf, cntl_cube_jja = get_season(CNTL,STARTYEAR,ENDYEAR)
expt_cube_djf, expt_cube_jja = get_season(EXPT,STARTYEAR,ENDYEAR)


cntl_summer_min_winter = cntl_cube_jja-cntl_cube_djf
expt_summer_min_winter = expt_cube_jja-expt_cube_djf
diff_cube_summ_min_winter = expt_summer_min_winter - cntl_summer_min_winter

#boundaries = [0.0, 5.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
boundaries=np.arange(-50,60,10)
cmap_use=plt.cm.get_cmap('RdBu_r',len(boundaries)-1)
#cmap_use.set_under('lightsteelblue')

# plot expt
try:
    expt_summer_min_winter.coord('longitude').guess_bounds()
    expt_summer_min_winter.coord('latitude').guess_bounds()
except:
    pass
grid_areas = iris.analysis.cartography.area_weights(expt_summer_min_winter)
#meandiff_djf = diff_cube_djf.collapsed(['longitude','latitude'],
#                               iris.analysis.MEAN, weights=grid_areas)
#diffchar_djf = str(np.around(meandiff_djf.data,2))

cs=iplt.pcolormesh(expt_summer_min_winter,cmap=cmap_use,
                 norm=mp.colors.BoundaryNorm(boundaries, 
                                             ncolors=len(boundaries)-1,
                                             clip=False))
cbar=plt.colorbar(cs,orientation='horizontal',extend='both')
cbar.set_ticks(boundaries)
cbar.set_label('degC')
titlename = EXPT + '-'  + 'JJA-DJF. Years:' + str(STARTYEAR) + '-' + str(ENDYEAR) 
plt.title(titlename, fontsize=10)
plt.gca().coastlines()
      
print('about to write to file')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja-djf_' + EXPT + '_' + FIELD + '.eps')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja-djf_' + EXPT + '_' + FIELD + '.png')
plt.close()


# plot anom
diff_cube_summ_min_winter 
try:
    diff_cube_summ_min_winter.coord('longitude').guess_bounds()
    diff_cube_summ_min_winter.coord('latitude').guess_bounds()
except:
    pass
grid_areas = iris.analysis.cartography.area_weights(diff_cube_summ_min_winter)
meandiff_jja = diff_cube_summ_min_winter.collapsed(['longitude','latitude'],
                               iris.analysis.MEAN, weights=grid_areas)
diffchar_jja = str(np.around(meandiff_jja.data,2))


#boundaries = [0.0, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
boundaries = np.arange(-18,22,4)
cmap_use=plt.cm.get_cmap('RdBu_r',len(boundaries)-1)
cs=iplt.pcolormesh(diff_cube_summ_min_winter,cmap=cmap_use,
                 norm=mp.colors.BoundaryNorm(boundaries, 
                                             ncolors=len(boundaries)-1,
                                             clip=False))
cbar=plt.colorbar(cs,orientation='horizontal',extend='both')
cbar.set_ticks(boundaries)
cbar.set_label('degC')
titlename = EXPT + '-' +  CNTL + '. Years:' + str(STARTYEAR) + '-' + str(ENDYEAR) + '.  JJA -DJF'
plt.title(titlename, fontsize=10)
plt.gca().coastlines()
      
print('about to write to file')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja-djf_' + EXPT + '-' + CNTL + '_' + FIELD + '.eps')
plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja-djf_' + EXPT + '-' + CNTL + '_' + FIELD + '.png')
plt.close()
