#!/usr/bin/env python2.7
#NAME
#    assess_spinup_temperature.py
#PURPOSE
#    This program will plot timeseries of:
#          ocean temperature
#          SST
#          SAT
#    to assess spinup
#
# it uses data provided by based on IDLPRGS/plotall_temptimeseries.pro
#
# search for 'main program' to find end of functions
# Julia 28/08/2019



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid



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
    


################################
# main program

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

