#NAME
#    plot temperature_and_salinity_by_basin_alt
#PURPOSE 
#
# this program will calculate temperature and salinity by basin
# using the alternative_mask file.
# These values are also obtained using more standard masks as part of the
# density calculations in plot_density_by_basin.py

# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import iris.plot as iplt
import gsw
import sys
import pandas as pd
from pathlib import Path
#from netCDF4 import Dataset, MFDataset
#from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

#============================================================================
def get_T_S(filename):
    """
    reads in temperature and salinity and converts to density (at 2000m)
    """

    # read in data
    T_cube = iris.load_cube(filename,'temp')
    #print(T_cube)
    S_cube = iris.load_cube(filename,'salinity')
    S_cube.data = (S_cube.data * 1000.) + 35.0 # convert to psu
    S_cube.attributes.pop("valid_min")
    S_cube.attributes.pop("valid_max")
    
    latitude = T_cube.coord('latitude').points
    longitude = T_cube.coord('longitude').points
    depth = T_cube.coord('depth_1').points
    latmesh,depmesh,lonmesh = np.meshgrid(latitude,depth,longitude)

   
    return T_cube,S_cube



def get_mask(filename,basin):
    """
    gets the mask from the Atlantic mask file on a v grid
    """

    cube = iris.load_cube(filename,basin + ' mask T_grid')

    return cube



####################################################################
def mask_function(cube,mask_cube):
    """
    masks the cubes
    """

    # mask V cube
    cube=iris.util.squeeze(cube)
    
    cubelev=np.copy(cube.data)
    depth_coord = cube.coord('depth_1').points
    for k in range(0,len(depth_coord)):
        cubelev[k,:,:] = np.ma.where(mask_cube.data > 0.5, cube.data[k,:,:],
                              -99999.)

    returncube = cube.copy(data=cubelev)
    returncube.data = np.ma.where(returncube.data > 1.0E20,-99999.,
                                  returncube.data)
    returncube.data.mask = np.where(returncube.data < -9999,1.0,0.0)
    #cube.data.mask = np.where(cube.data ==-99999., 1.0,0.0)
    

    return returncube


def basin_avg(cube):
    """
    calculates the basin average across the cube
    """

    cubedata = cube.data
    #print(cubedata.mask)
    zonalmean_data = np.ma.mean(cube.data,axis=2)
    #print(np.shape(cube.data))
    #print(np.shape(zonalmean_data))
    basin_avg_cube = cube.collapsed('longitude',iris.analysis.MEAN)
    basin_avg_cube = basin_avg_cube.copy(data=zonalmean_data)
    #print(zonalmean_data)
    #print(basin_avg_cube)
    #print(basin_avg_cube.data[0,90])
    #print(zonalmean_data[0,90])
    #print(cube.data[0,90,:])
    #print(cube.data.mask[0,90,:])
    #sys.exit(0)

    return basin_avg_cube




#####################################################################
def process_data(filename,basin,year):
    """
    processess all the data for each year
    """


    # file where V and dz are stored on density coordinates

    # 1. calculate  Temp, Salinity grid
    (temperature_cube,
     salinity_cube) = get_T_S(filename)

    # gets the basins over which we calculate
    mask_cube=get_mask('alternative_mask.nc',basin) # get mask on T grid

    temperature_cube = mask_function(temperature_cube,mask_cube)
    salinity_cube = mask_function(salinity_cube,mask_cube)

    # calculate stream function for the basin
    temperature_basin_avg_cube = basin_avg(temperature_cube)
    salinity_basin_avg_cube = basin_avg(salinity_cube)


  
    # save everything to a file

    fileout = (FILEINIT + 'um/' + exptname +
               '/basin_diagnostics/T_S_' + exptname + '_' +
               basin + str(year) + '.nc')

    temperature_basin_avg_cube.long_name='temperature basin'
    salinity_basin_avg_cube.long_name='salinity basin'
    iris.save([temperature_basin_avg_cube,
               salinity_basin_avg_cube],
               fileout,fill_value = -99999.)


    



#######################################################################
exptname = 'xqbwg'
startyear=3901
endyear=4000
basin='Indo_Pacific' # Indo_Pacific or Atlantic  
FILEINIT = '/uolstore/Research/a/hera1/earjcti/'
FILEINIT = '/home/earjcti/'

period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400','xpsig':'EP',
          'xpsic':'PI','xqbwd':'LP','xqbwj':'LP490','xqbwe':'EP400',
          'xqbwg':'EP',
          'xqbwc':'PI','xqbwn':'LP_warm_NH_JJA','xqbwo':'LP_warm_SH_DJF',
          'xqbwp':'EP_warm_NH_JJA','xqbwq':'EP_warm_SH_DJF'}

# get individual years diagnostics for the basin
for year in range(startyear,endyear):
    filestart = FILEINIT + 'um/' + exptname + '/pg/' + exptname 
    filename = filestart + 'o#pg00000'+str(year).zfill(4)+'c1+.nc'
    print(filename)
    process_data(filename,basin,year)

#################################################
# get the mean dianostics for the basin
sal_cubelist = CubeList([])
temp_cubelist = CubeList([])

for year in range(startyear,endyear):
     filename = (FILEINIT + 'um/' + exptname +
                 '/basin_diagnostics/T_S_' + exptname + '_' +
                 basin + str(year) + '.nc')#

     temp_cubelist.append(iris.load_cube(filename,'temperature basin'))
     sal_cubelist.append(iris.load_cube(filename,'salinity basin'))
     sal_cube = iris.load_cube(filename,'salinity basin')
   
iris.util.equalise_attributes(sal_cubelist)
iris.util.equalise_attributes(temp_cubelist)
sal_cubes = sal_cubelist.merge_cube()
temp_cubes = temp_cubelist.merge_cube()

sal_avg_cube = sal_cubes.collapsed('t',iris.analysis.MEAN)
temp_avg_cube = temp_cubes.collapsed('t',iris.analysis.MEAN)

fileout = (FILEINIT + 'um/' + exptname +
           '/basin_diagnostics/T_S_mean_' + exptname + '_' +
           basin + str(startyear) + '_' + str(endyear-1)+'.nc')
iris.save([temp_avg_cube,sal_avg_cube],
          fileout,fill_value = -99999.)



