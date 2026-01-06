#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sep 21 2019

@author: earjcti

This will plot the MMM precipitation and temperature anomalies.  It will hatch the
statistically robust changes based on the followin method (Mba et al 2018)

1. Roughly 80% of models must agree on the direction of the change
   (if we have 13 models then 10 must agree)
2. The ((ensemble mean change) / (the ensemble standard deviation))>1.
   (note I am not sure whether to get the ensemble standard deviation from
    pi or the mPWP maybe I will try both)

"""

import sys
import iris
import iris.quickplot as qplt
import iris.plot as iplt
#from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import netCDF4


def get_data(filereq, field, modeluse):
    """
    gets the field (field) from the file (filereq) and loads it
    into an iris cube (the model name is in modeluse)
    outputs a cube of the data that is as simple as possible
    """

    if modeluse == 'MMM':
        cube = iris.load_cube(filereq, field)
    else:
        cubes = iris.load(filereq)
        cube = cubes[0]
    cube.data = cube.data.astype('float32')

    if field == 'SST' or field == 'NearSurfaceTemperature':
        if (modeluse == 'MIROC4m' or modeluse == 'COSMOS'):
            cube.units = 'Celsius'
        else:
            cube.convert_units('Celsius')

    for coord in cube.coords():
        name = coord.standard_name
        if name != 'latitude' and name != 'longitude':
            if name is None:
                if coord.long_name is None:
                    cube.remove_coord(coord.var_name)
                else:
                    cube.remove_coord(coord.long_name)
            else:
                cube.remove_coord(name)

    if modeluse == 'EC-Earth3.1' and field == 'SST':
        cube.coord('latitude').bounds = None
        cube.coord('longitude').bounds = None

    cube.cell_methods = None

    return cube



def get_ind_model_anomaly():
    """
    gets the anomaly data for all the models puts them into a list of cubes for
    returning to the calling program
    """

    cubelist = iris.cube.CubeList([])
    for i, model in enumerate(MODELNAMES):
        cubepi = get_data(FILESTART + model + '/E280.' + FIELDNAME + '.allmean.nc',
                          FIELDNAME, model)
        cubeplio = get_data(FILESTART + model + '/EOI400.' + FIELDNAME + '.allmean.nc',
                            FIELDNAME, model)
        cubediff = cubeplio - cubepi
        cubelist.append(cubediff)


    return cubelist

def check_sign(cubelist):
    """
    we are checking that a number (NSIGN) of cubes has the same sign in a field
    for example if the field is temperature we may want to see that 10/13 models all show the
    same sign of temperature change

    input a list of cubes from models
    output a numpy array which is 1 where at least 'NSIGN' models have the same sign
           and 0 when they don't
    """

    cube = cubelist[0]
    ydim, xdim = np.shape(cube.data)
    posarr = np.zeros((ydim, xdim))
    negarr = np.zeros((ydim, xdim))

    for cubeno, cube in enumerate(cubelist):
        cubedata = cube.data
        temparr = np.ma.where(cubedata > 0, 1, 0) # temporary array
        posarr = posarr + temparr
        temparr = np.ma.where(cubedata < 0, 1, 0)
        negarr = negarr + temparr

    # if there are more than NSIGN elements in posarr or negarr than the same sign arr is set
    sign_arr = np.ma.where(posarr >= NSIGN, 1, 0)
    temparr = np.ma.where(negarr >= NSIGN, 1, 0)
    sign_arr = sign_arr + temparr
    newcube = cube.copy(data=sign_arr)

    return newcube

def check_large_anomaly(ratiocube):
    """
    if ratio is greater than 1 set to 1, otherwise set to zero
    """
    temparr = np.ma.where(ratiocube.data > 1, 1, 0)
    newcube = ratiocube.copy(data=temparr)

    return newcube


def main():
    """
    1. get the change in the field from all the models
    2. find out where > 70% of the models agree in the sign of the change (region1)
    3. get mean change in the field
    4. get the standard deviation for the control climate
    5. find out where mean change / standard deviation > 1 (region2)
    6. plot, hatching where region 1 and region2 are satisfied
    """

    namefield = {"NearSurfaceTemperature" : "SAT",
                 "TotalPrecipitation" : "Precipitation",
                 "SST" : "SST"
                 }

    modelcubelist = get_ind_model_anomaly()
    cube_sign = check_sign(modelcubelist) # checks where a given number of them are the correct sign
    meancube = iris.load_cube(FILESTART + FIELDNAME + '_multimodelmean.nc',
                              FIELDNAME + 'mean_anomaly')  # get mean
    stdevcube = iris.load_cube(FILESTART + FIELDNAME + '_multimodelmean.nc',
                               FIELDNAME + 'std_pi') # get standard deviation
    cube_mean_std_sign = check_large_anomaly(meancube / stdevcube)

    #qplt.contourf(cube_sign)
    ax = plt.axes(projection = ccrs.PlateCarree())
    if FIELDNAME == 'TotalPrecipitation':
       
        qplt.contourf(meancube, np.arange(-1.4, 1.6, 0.2), cmap='RdBu', extend='both')
        plt.figtext(0.02, 0.97,'d)',
                   horizontalalignment='left',
                   verticalalignment='top',
                   fontsize=20)
    else:
        qplt.contourf(meancube, np.arange(0, 5.5, 0.5), cmap='Reds', extend='both')
    iplt.contourf(cube_sign, 1, hatches=[None, '///'], colors='none')
    iplt.contourf(cube_mean_std_sign, 1, hatches=[None, 3 * '\\\''], colors='none')
    plt.gca().coastlines()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-', draw_labels='true')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    plt.title(namefield.get(FIELDNAME) +' anomaly: multimodel mean')

    plt.savefig(OUTSTART + FIELDNAME + 'robust_anomalies.eps')
    plt.savefig(OUTSTART + FIELDNAME + 'robust_anomalies.pdf')
    plt.close()

    OUTNC = FILESTART + 'dummy.nc'
    if FIELDNAME == 'NearSurfaceTemperature':
        OUTNC = FILESTART + 'alldata/data_for_fig2.nc'
    if FIELDNAME == 'TotalPrecipitation':
        OUTNC = FILESTART + 'alldata/data_for_fig5d.nc'
        
    print(OUTNC)
    
    cubelist = iris.cube.CubeList([meancube, cube_sign, cube_mean_std_sign])
    iris.save(cubelist, OUTNC)        
    


    return



# variable definition
LINUX_WIN = 'l'

#FIELDNAME = 'NearSurfaceTemperature'
FIELDNAME = 'SST'
UNITS = 'deg C'

#FIELDNAME = 'TotalPrecipitation'
#UNITS = 'mm/day'


if LINUX_WIN == 'l':
    FILESTART = ('/nfs/hera1/earjcti/regridded/')
    OUTSTART = FILESTART + 'allplots/' + FIELDNAME + '/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
    OUTSTART = ' '


MODELNAMES = ['CCSM4-Utr', 'COSMOS', 'CESM1.2', 'CESM2','CCSM4',
              'EC-Earth3.3', 'GISS2.1G', 'HadCM3',
              'IPSLCM6A', 'IPSLCM5A2', 'IPSLCM5A',
              'MIROC4m', 'MRI2.3',
              'NorESM-L', 'NorESM1-F',
              'CCSM4-UoT'
             ]
#ODELNAMES = ['EC-Earth3.1']
NSIGN = np.floor(len(MODELNAMES) * 0.8) # *0.8 is 80%

main()
