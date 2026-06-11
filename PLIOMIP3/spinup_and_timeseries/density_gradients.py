#!/usr/bin/env python2.7
#NAME
#    density_gradients
#PURPOSE
#
#   This program will do a timeseries of the mean density southwards of 65S for
#   a) surface
#   b) top 100m
#   c) top 200m
#   d) top 1000m
#   e) 1000m-3000m
#
#  the density has previously been obtained and is in basin_diagnostics
#
#
#
# Julia 22/11/2025
#
# update Julia 30.1.2026:  we will just write out a timeseries of density
# on single levels
#
#  5m, 25m, 47.8m 203m, 447m, 666m, 995m, 1500m
#
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
from cartopy.mpl.ticker import LongitudeFormatter
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess
import matplotlib.ticker as mticker

exptname='xqbws'
startyear=3000
endyear=4000
filestart = '/uolstore/Research/a/hera1/earjcti/um/' + exptname 
latmax=-60
single_levs = 'y'  # y write out on levels
                   # n write out bands of levels
 

 
def get_density_band(year):
    """
    gets the density at each of the bands
    #   a) surface
    #   b) top 100m
    #   c) top 200m
    #   d) top 1000m
    #   e) 1000m-3000m

    """
 
    filename=(filestart+'/basin_diagnostics/'+exptname+'_Pacific'+
              str(year) +'.nc')
    print(filename)
    cubes = iris.load(filename)
    # extract 3d density
    denscube = cubes.extract('density (calculated from T and S')[0]
    # density southwards of 65S
    south_constraint = iris.Constraint(latitude=lambda v:v<latmax)
    dens_south_cube =denscube.extract(south_constraint)

    # average over area
    lat = dens_south_cube.coord('latitude')
    lon = dens_south_cube.coord('longitude')
    depth = dens_south_cube.coord('depth_1')
    if not lat.has_bounds():
        lat.guess_bounds()
    if not lon.has_bounds():
        lon.guess_bounds()
    if not depth.has_bounds():
        depth.guess_bounds()


    weights = iris.analysis.cartography.area_weights(dens_south_cube)

    # mean by depth
    dens_levels_cube = dens_south_cube.collapsed(['latitude','longitude'],
                                                 iris.analysis.MEAN,
                                                 weights=weights)

    surf_data = dens_levels_cube[0].data

    # extract cube by levels
    cube_100 = dens_levels_cube.extract(iris.Constraint(depth_1=lambda d:d<100))
    cube_200 = dens_levels_cube.extract(iris.Constraint(depth_1=lambda d:d<200))
    cube_1000 = dens_levels_cube.extract(
        iris.Constraint(depth_1=lambda d:d<1000))
    cube_1000_3000 = dens_levels_cube[13:16]
    
    
    thickness = np.diff(depth.bounds,axis=1).squeeze()
    
    avg_100=cube_100.collapsed('depth_1',iris.analysis.MEAN,
                               weights=thickness[0:len(cube_100.data)])
    avg_200=cube_200.collapsed('depth_1',iris.analysis.MEAN,
                               weights=thickness[0:len(cube_200.data)])
    avg_1000=cube_1000.collapsed('depth_1',iris.analysis.MEAN,
                               weights=thickness[0:len(cube_1000.data)])
    avg_1000_3000=cube_1000_3000.collapsed('depth_1',iris.analysis.MEAN,
                               weights=thickness[13:16])
    
    data_100=avg_100.data
    data_200=avg_200.data
    data_1000=avg_1000.data
    data_1000_3000=avg_1000_3000.data
    
  
    return (surf_data,data_100,data_200,data_1000,data_1000_3000)



def get_density_level(year):
    """
    gets the density at each of the levels
    #  5m, 25m, 47.8m 203m, 447m, 666m, 995m, 1500m

    """
 
    filename=(filestart+'/basin_diagnostics/'+exptname+'_Pacific'+
              str(year) +'.nc')
    print(filename)
    cubes = iris.load(filename)
    # extract 3d density
    denscube = cubes.extract('density (calculated from T and S')[0]
    # density southwards of 65S
    south_constraint = iris.Constraint(latitude=lambda v:v<latmax)
    dens_south_cube =denscube.extract(south_constraint)

    # average over area
    lat = dens_south_cube.coord('latitude')
    lon = dens_south_cube.coord('longitude')
    depth = dens_south_cube.coord('depth_1')
    if not lat.has_bounds():
        lat.guess_bounds()
    if not lon.has_bounds():
        lon.guess_bounds()
    if not depth.has_bounds():
        depth.guess_bounds()


    weights = iris.analysis.cartography.area_weights(dens_south_cube)

    # mean by depth
    dens_levels_cube = dens_south_cube.collapsed(['latitude','longitude'],
                                                 iris.analysis.MEAN,
                                                 weights=weights)

    surf_data = dens_levels_cube[0].data
    #print(dens_levels_cube)
    #print(dens_levels_cube.coord('depth_1').points)
   
    # extract cube by levels
 
    data_5= (dens_levels_cube.extract(iris.Constraint(depth_1=5.0))).data
    data_25= (dens_levels_cube.extract(iris.Constraint(depth_1=25.0))).data
    data_47= (dens_levels_cube.extract(iris.Constraint(depth_1=47.85))).data
    data_203= (dens_levels_cube.extract(iris.Constraint(depth_1=203.7))).data
    data_447= (dens_levels_cube.extract(iris.Constraint(depth_1=447.05))).data
    data_666= (dens_levels_cube.extract(iris.Constraint(depth_1=666.3))).data
    data_995= (dens_levels_cube.extract(iris.Constraint(depth_1=995.55))).data
    data_1500=( dens_levels_cube.extract(iris.Constraint(depth_1=1500.85))).data

    
  
    return (data_5,data_25,data_47,data_203,data_447,data_666,data_995,data_1500)




################################
# main program

all_years = []

if single_levs == 'n': # we will get density bands
    all_surf = []
    all_100 = []
    all_200 = []
    all_1000 = []
    all_1000_3000 = []

    for year in range(startyear,endyear):
        surf_data, data_100,data_200,data_1000,data_1000_3000=get_density_band(year)
        all_years.append(year)
        all_surf.append(surf_data)
        all_100.append(data_100)
        all_200.append(data_200)
        all_1000.append(data_1000)
        all_1000_3000.append(data_1000_3000)

    # Write to a text file
    fileout = (filestart+'/basin_diagnostics/'+exptname+'_Pacific_density'+
               str(startyear) +'_' + str(endyear) + '_'+str(latmax)+ '.txt')
 
    with open(fileout, "w") as f:
        string = ("Year, dens_surf,dens_0_100,dens_0_200," +
                  "dens_0_1000,dens_1000_3000\n") 
        f.write(string)

        for y, a, b, c, d, e in zip(all_years, all_surf, all_100,all_200,all_1000,all_1000_3000):
            f.write(f"{y}, {a}, {b},{c},{d},{e}\n")


    #plot
        
    plt.figure(figsize=(8, 5))
    plt.plot(all_years, all_surf, label='surface')
    plt.plot(all_years, all_100, label='0-100')
    plt.plot(all_years, all_200, label='0-200')
    plt.plot(all_years, all_1000, label='0-1000')
    plt.plot(all_years, all_1000_3000, label='1000-3000')

    # Add labels and title
    plt.xlabel('Year')
    plt.ylim(35.9,36.9)
    plt.ylabel('density')
    plt.legend()
    plt.title('Year vs density')

    # Optional: Add grid
    plt.grid(True)

    # Show the plot
    plt.show()





if single_levs == 'y': # we will get density on a range of levels
    #  5m, 25m, 47.8m 203m, 447m, 666m, 995m, 1500m

    all_5m = []
    all_25m = []
    all_47m = []
    all_203m = []
    all_447m = []
    all_666m = []
    all_995m = []
    all_1500m = []

    for year in range(startyear,endyear):
        d5,d25,d47,d203,d447,d666,d995,d1500=get_density_level(year)
        all_years.append(year)
        all_5m.append(d5)
        all_25m.append(d25)
        all_47m.append(d47)
        all_203m.append(d203)
        all_447m.append(d447)
        all_666m.append(d666)
        all_995m.append(d995)
        all_1500m.append(d1500)
       
    # Write to a text file
    fileout = (filestart+'/basin_diagnostics/'
               +exptname+'_Pacific_density_levs'+
               str(startyear) +'_' + str(endyear) + '_'+str(latmax)+ '.txt')
 
    with open(fileout, "w") as f:
        string = ("Year, 5m, 25m, 47m 203m, 447m, 666m, 995m, 1500m\n") 
        f.write(string)
        print('here',string)
        for y, a, b, c, d, e, g, h, i in zip(all_years, all_5m,
                                             all_25m,all_47m,
                                             all_203m,all_447m,
                                             all_666m,all_995m,all_1500m):
            f.write(f"{y}, {a}, {b},{c},{d},{e}, {g}, {h}, {i}\n")


    #plot
    print(all_666m)
    print(all_995m)
    plt.figure(figsize=(8, 5))
    plt.plot(all_years, all_5m, label='5m')
    plt.plot(all_years, all_47m, label='47m')
    plt.plot(all_years, all_447m, label='447m')
    plt.plot(all_years, all_666m, label='666m')
    plt.plot(all_years, all_995m, label='995m')
    
    # Add labels and title
    plt.xlabel('Year')
    plt.ylim(36,37)
    plt.ylabel('density')
    plt.legend()
    plt.title('Year vs density')
    
    # Optional: Add grid
    plt.grid(True)
    
    # Show the plot
    plt.show()


####

