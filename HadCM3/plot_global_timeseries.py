#!/usr/bin/env python3.8
#NAME
#    TEMPERATURE_TIMESERIES
#PURPOSE
#    this program will write out the ocean temperature timeseries from 
#    the pg file and plot it.
#
# search for 'main program' to find end of functions
# Julia 14/1/2017



import sys
import numpy as np
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import iris.plot as iplt
import iris.quickplot as qplt


def get_global_avg(year):
    """
    gets the global average value of the field for the given year
    """
    cent_ind = {0: "0", 1: "1", 2: "2", 3: "3", 4: "4", 5: "5", 6: "6", 
                7: "7", 8: "8", 9: "9", 10: "a", 11: "b", 12: "c", 13: "d", 
               14: "e", 15: "f", 16: "g", 17: "h", 18: "i", 19: "j", 20: "k", 
               21: "l", 22: "m", 23: "n", 24: "o", 25: "p", 26: "q", 27: "r", 
               28: "s", 29: "t", 30: "u", 31: "v", 32: "w", 33: "x", 34: "y", 
               36: "z"}

    if NEW_VERSION:
        filename = FILESTART + '#' + 'pg' + str(year).zfill(9) + 'c1+.nc'
        print(filename)
    else:
        cent = np.floor(year/100.)
        filename = (FILESTART + '@pg' + cent_ind.get(cent) + 
                    str(np.int(year - (cent*100))).zfill(2) + 'c1.nc')
       
    cube = iris.load_cube(filename,FIELD)
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    cube.coord('depth_1').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    meanval = cube.collapsed(['depth_1','latitude','longitude'],
                             iris.analysis.MEAN, weights=grid_areas)

    if FIELD == 'salinity':
        globmean = (meanval.data * 1000.)+35.0
    else:
        globmean = meanval.data

    
    return globmean[0]

################################
# main program

EXPT = 'xpkma'
NEW_VERSION = True
STARTYEAR=1851
NYEARS=500
FIELD = 'temp'    # valid fields are 'salinity', 'temp'

FILESTART = '/nfs/hera1/earjcti/um/' + EXPT + '/pg/' + EXPT + 'o'
meandata = np.zeros(NYEARS)
for year in range(STARTYEAR, STARTYEAR+NYEARS):
    globavgval = get_global_avg(year)
    meandata[year-STARTYEAR] = globavgval


plt.plot(meandata)
plt.xlabel('year')
plt.title(FIELD)
if FIELD == 'salinity':
    plt.ylabel('psu')
if FIELD == 'temp':
    plt.ylabel('degC (potential temperature)')
plt.savefig('/nfs/hera1/earjcti/HadCM3_plots/OceanTimeseries/' + EXPT + '_' + FIELD + '.eps')
plt.savefig('/nfs/hera1/earjcti/HadCM3_plots/OceanTimeseries/' + EXPT + '_' + FIELD + '.png')
plt.close()
