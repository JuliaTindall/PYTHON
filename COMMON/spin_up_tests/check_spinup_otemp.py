#!/usr/bin/env python3.6
#NAME
#    assess_spinup_temperature.py
#PURPOSE
#    This program will plot timeseries of:
#          ocean temperature
#          SST
#          SAT
#    to assess spinup
#
# it uses data provided by based on HadGEM2/assess_spinup_temperature.py
#
# search for 'main program' to find end of functions
# Julia 2023



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris
import sys



def plotdata(filename,titlename,plotno,valmin,valmax):
    """ plots the data in a file
    requires the filename, the title of the plot
    where on the page we want the plot to go
    """
    
    
    f = open(filename, "r")
    line = f.readline()
    expt=-1
    first = 'y'
    allindex = []
    alldata = []
    alltitle = []
    

    while line:
        # test for titleline if 1st character not space use this
    
        if line[0:3] != '   ':
            expt = expt + 1
            alltitle.append(line)
        
            if first == 'n':
                allindex.append(index_list)
                alldata.append(data_list)
            
            index_list=[]
            data_list=[]
            first = 'n'
        else:
            index, data = line.split()
        
            index_list.append(np.int(index))
            data_list.append(np.float(data))
 
        
        line = f.readline()

    f.close()

    allindex.append(index_list)
    alldata.append(data_list)
            
    index_list=[]
    data_list=[]

    print(alltitle)

    ax = plt.subplot(2,2,plotno)
    colors=['black','blue','red']
    for i in range(0,len(allindex)):
        ax.plot(allindex[i],alldata[i],label=alltitle[i],color=colors[i])
    
   
    ax.set_ylabel('deg C')
    ax.set_xlabel('year')
    ax.set_ylim(valmin,valmax)
    ax.set_xlim(0,550)
    ax.set_title(titlename)
        
    ax.legend()
    

def get_global_field():
    """
    gets the mean field as a timeserie
    both averaged over the ocean and for different levels
    """
    extra = {'10':'a', '11':'b', '12':'c', '13':'d', '14':'e', '15':'f',
             '16':'g', '17':'h', '18':'i', '19':'j', '20':'k', '21':'l',
             '22':'m', '23':'n', '24':'o', '25':'p', '26':'q', '27':'r',
             '28':'s'}

    filestart = ('/nfs/hera1/earjcti/um/' + EXPTNAME + '/pg/' + EXPTNAME
                 + 'o@pg')
    ann_timeseries = np.zeros(NYEARS)
    lev_timeseries = np.zeros((NYEARS, 20))
    for year in range(STARTYEAR, STARTYEAR + NYEARS):
        cent = np.floor(year / 100.)
        yearuse = np.str(np.int(year - (cent * 100))).zfill(2)
        centind = extra.get(np.str(np.int(cent)))
   
        filename = filestart + centind + yearuse + 'c1.nc'
        # read in data
        cube = iris.load_cube(filename, FIELDNAME)

        # average
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube)
        levavg_cube = (cube.collapsed(['longitude', 'latitude'],
                                     iris.analysis.MEAN, weights = grid_areas))
        lev_timeseries[year-STARTYEAR,:] = levavg_cube.data
        
        globavg_cube = cube.collapsed(['longitude','latitude','depth_1'],
                                      iris.analysis.MEAN,
                                      weights=grid_areas)
        print(globavg_cube.data, np.mean(levavg_cube.data))
      
        #####NOTE PUT GRID AREAS2 TO BE  weighed by depth"
        sys.exit(0)

                                                   
################################
# main program

# read in data from files
FIELDNAME = 'POTENTIAL TEMPERATURE (OCEAN)  DEG.C'
STARTYEAR = 2001
NYEARS=600
EXPTNAME = 'xozzb'


field_ts, field_ts_levels = get_global_field()


fig=plt.figure(figsize=(8,8))

#hADgem
filename = '/nfs/hera1/earjcti/IDLPLOTS/HadGEM/plotall_temptimeseries.txt'
titlename = 'a) HadGEM total ocean temperature'
plotdata(filename,titlename,1,-0.1,1.0)


# HadCm3
filename = '/nfs/hera1/earjcti/IDLPLOTS/HadGEM/plotall_temptimeseries_HadCM3.txt'
titlename = 'b) HadCM3 total ocean temperature'
plotdata(filename,titlename,2,-0.1,1.0)

# HadGEM level 1
filename = '/nfs/hera1/earjcti/IDLPLOTS/HadGEM/plotall_temptimeseries_lev1.txt'
titlename = 'c) HadGEM sea surface temperature'
plotdata(filename,titlename,3,-0.5,0.5)

# HadGEM SAT
filename = '/nfs/hera1/earjcti/IDLPLOTS/HadGEM/plotall_SATtimeseries.txt'
titlename = 'd) HadGEM SAT'
plotdata(filename,titlename,4,-1.0,1.0)



fig.tight_layout()
#plt.show()
fileout='/home/earjcti/PYTHON/PLOTS/HadGEM2/temperature_trends.eps'
plt.savefig(fileout)
plt.close()


####

