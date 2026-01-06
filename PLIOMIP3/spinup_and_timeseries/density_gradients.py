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
# Julia 22/11/2016



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

exptname='xpsie'
startyear=12
endyear=1999
filestart = '/uolstore/Research/a/hera1/earjcti/um/' + exptname 
latmax=-65




def get_density(year):
    """
    gets the density at each of the levels
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


    


################################
# main program

all_years = []
all_surf = []
all_100 = []
all_200 = []
all_1000 = []
all_1000_3000 = []

for year in range(startyear,endyear):
    surf_data, data_100,data_200,data_1000,data_1000_3000=get_density(year)
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


####

