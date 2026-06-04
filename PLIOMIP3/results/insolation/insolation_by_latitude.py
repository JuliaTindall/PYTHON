"""
#NAME
#    insolation_by_latitude.py
#PURPOSE 
#
#  This program is wanting to look at the orbital experments in order
#  to see the insolation by latitude.  In particular I want to explain the
#  cooling over the southern ocean in the highJJA - highDJF orbit
#
#  if there is more insolation averaged over the south in the DJF orbit this 
#  will explain it
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
import sys
#from netCDF4 import Dataset, MFDataset


def get_zonal_mean_incoming_sw(expt):
    """
    gets the zonal mean annual average incoming sw radiation
    """

    filestart = '/home/earjcti/um/' + expt + '/pcpd/' + expt + 'a#pd000003900*'

    if season == 'jja' or season == 'JJA':
        months = ['jn','jl','ag']

    if season == 'djf' or season == 'DJF':
        months = ['ja','fb','dc']

    if season == 'son' or season == 'SON':
        months = ['sp','ot','nv']

    if season == 'djf' or season == 'MAM':
        months = ['mr','ar','my']

    if season == 'all' or season == 'ALL':
        months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']

    moncubelist = CubeList([])
    for mon in months:
        cube = iris.load_cube(filestart + mon + '+.nc',
                              'INCOMING SW RAD FLUX (TOA): ALL TSS')
        moncubelist.append(cube)

    iris.util.equalise_attributes(moncubelist)
    cubes = moncubelist.concatenate_cube()
    avgcube = cubes.collapsed('t',iris.analysis.MEAN)

    # it is insolation so same at every longitude hence we can just
    # extract one longitude
    zm_cube = avgcube[0,:,0]

    return zm_cube
#=================================================================
# MAIN PROGRAM STARTS HERE

EXPT1 = 'xqbwn'
EXPT2 = 'xqbwd'
EXPT3 = 'xqbwo'

season = 'MAM'   # DJF JJA all

if season == 'DJF' or season == 'JJA' or season == 'SON' or season == 'MAM':
    print('USE WITH EXTREME CAUTION')
    print('RESULTS ARE INDICATIVE ONLY')
    print('YOU SHOULD BE USING DAILY DATA')
    input('PRESS ENTER TO CHECK YOU UNDERSTAND')



incomsw_1_cube = get_zonal_mean_incoming_sw(EXPT1)
incomsw_2_cube = get_zonal_mean_incoming_sw(EXPT2)
incomsw_3_cube = get_zonal_mean_incoming_sw(EXPT3)

plt.figure(figsize=[10,5])
plt.subplot(1,2,1)
plt.plot(incomsw_1_cube.data,incomsw_1_cube.coord('latitude').points,label=EXPT1)
plt.plot(incomsw_2_cube.data,incomsw_2_cube.coord('latitude').points,label=EXPT2)
plt.plot(incomsw_3_cube.data,incomsw_3_cube.coord('latitude').points,label=EXPT3)
plt.legend()

plt.subplot(1,2,2)
plt.plot(incomsw_1_cube.data - incomsw_3_cube.data,incomsw_1_cube.coord('latitude').points,label=EXPT1 + '-' + EXPT3)
plt.plot(incomsw_1_cube.data - incomsw_2_cube.data,incomsw_1_cube.coord('latitude').points,label=EXPT1 + '-' + EXPT2)
#plt.hline(x=0)
plt.legend()
plt.title(season)
plt.show()
