"""
#NAME
#    Dominant_vegetation
#PURPOSE 
#
#  This program will plot the dominant vegetation.  This will be defined as any 
#  vegetation over a tunable threshold (start with 50\%)
"""

# Import necessary libraries

import os
import numpy as np
import math
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.quickplot as qplt
import iris.plot as iplt
from netCDF4 import Dataset, MFDataset
import sys
from matplotlib.colors import ListedColormap, BoundaryNorm

#from netCDF4 import Dataset, MFDataset

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def get_dominant_vegetation(file,threshold,dyn_stat_ind):
    """
    for each gridbox gets the vegetation in the file that is more than 
    threshold
    """

    if dyn_stat_ind == 's':
        fullcube = iris.load_cube(file)
        cube = fullcube[0,:,:,:]
    else:
        fullcube = iris.load_cube(file)
        cube = fullcube.collapsed('time',iris.analysis.MEAN)
        
    
    outcube = cube[0,:,:]    # a shell for storing output data

    data = cube.data
    nlats = len(cube.coord('latitude').points)
    nlons = len(cube.coord('longitude').points)
                
    outdata = np.ma.masked_all((nlats,nlons)) # dominant
    maxdata = np.ma.masked_all((nlats,nlons)) # maximum pft
    
    for i in range(0,nlons):
        for j in range(0,nlats):
            if (not data.mask[0,j,i]): # not ocean so check above threshold
                outdata[j,i]=0  # assume non found
                outdata.mask[j,i]=False
                # max value
                if np.argmax(data[:,j,i]) < 6: # remove urban
                    maxdata[j,i] = np.argmax(data[:,j,i]) + 1
                else:
                    maxdata[j,i] = np.argmax(data[:,j,i]) 
                
                # dominant pft
                for k in range(0,9):
                    if (data[k,j,i] > threshold):
                        if np.argmax(data[:,j,i]) < 6: # remove urban
                            outdata[j,i] = k + 1
                        else:
                            outdata[j,i] = k
                    

                # add in mixed forests
                if outdata[j,i] == 0:
                    if (data[0,j,i] + data[1,j,i] > threshold):
                        outdata[j,i]=9

                # add in mixed grasses and shrubs
                if outdata[j,i] == 0:
                    if (data[2,j,i] + data[3,j,i]
                        + data[4,j,i] > threshold):
                        outdata[j,i]=10

                # look at remainder
                if outdata[j,i] == 0:
                    outdata[j,i]=11
                                     
                    #sys.exit(0)
                
    outdata = np.ma.where(outdata > 1E20,-99999,outdata)
    outcube.data = outdata
    outcube.long_name = 'PFT > ' + str(threshold)
    
    maxcube = outcube.copy(data=maxdata)
    maxcube.long_name = 'most common PFT'
    
    #iris.save([outcube,maxcube],dyn_stat_ind + '_temporary.nc',fill_value =-99999)
    
    return maxcube, outcube               


def get_domveg_pi(threshold):
    """
    for each gridbox gets the vegetation in the file that is more than 
    threshold
    """

    filename = filestart + 'qrfractype.PMIP_'

    allcubes = CubeList()
    for i in range(1,10):
        vcube = iris.load_cube(filename + str(i) + '.nc')
        vcube.coord('Surface').points = [i]
        allcubes.append(vcube)

    iris.util.equalise_attributes(allcubes)
    acube = allcubes.concatenate_cube()
    cube = acube[0,:,:,:]
    
    outcube = cube[0,:,:]    # a shell for storing output data

    data = cube.data
    nlats = len(cube.coord('latitude').points)
    nlons = len(cube.coord('longitude').points)
                
    outdata = np.ma.masked_all((nlats,nlons)) # dominant
    maxdata = np.ma.masked_all((nlats,nlons)) # maximum pft
    
    for i in range(0,nlons):
        for j in range(0,nlats):
            if (not data.mask[0,j,i]): # not ocean so check above threshold
                outdata[j,i]=0  # assume non found
                outdata.mask[j,i]=False
                # max value
                if np.argmax(data[:,j,i]) < 6: # remove urban
                    maxdata[j,i] = np.argmax(data[:,j,i]) + 1
                else:
                    maxdata[j,i] = np.argmax(data[:,j,i]) 
                
                # dominant pft
                for k in range(0,9):
                    if (data[k,j,i] > threshold):
                        if np.argmax(data[:,j,i]) < 6: # remove urban
                            outdata[j,i] = k + 1
                        else:
                            outdata[j,i] = k
                    

                # add in mixed forests
                if outdata[j,i] == 0:
                    if (data[0,j,i] + data[1,j,i] > threshold):
                        outdata[j,i]=9

                # add in mixed grasses and shrubs
                if outdata[j,i] == 0:
                    if (data[2,j,i] + data[3,j,i]
                        + data[4,j,i] > threshold):
                        outdata[j,i]=10

                # look at remainder
                if outdata[j,i] == 0:
                    outdata[j,i] =11
                    
                    #sys.exit(0)
                
    outdata = np.ma.where(outdata > 1E20,-99999,outdata)
    outcube.data = outdata
    outcube.long_name = 'PFT > ' + str(threshold)
    
    maxcube = outcube.copy(data=maxdata)
    maxcube.long_name = 'most common PFT'
    
    #iris.save([outcube,maxcube],dyn_stat_ind + '_temporary.nc',fill_value =-99999)
    
    return maxcube, outcube               


def plot_cube(cube1,cube2,cube3,title1,title2,title3,fileout,type):
    """
    plots the biomes in a nice way
    """
 
    biome1 = cube1.data.copy()
    biome2 = cube2.data.copy()
    biome3 = cube3.data.copy()
    lpmask = np.isnan(biome1)

    if type == 'm':
        cmap = ListedColormap([
            'forestgreen',   # 1 broadleaf trees
            'darkgreen',     # 2 needleleaf trees
            'yellowgreen',    # 3 C3 grass
            'gold',          # 4 C4 grass
            'sienna',        # 5 shrub
            'blue',          # 6 lake
            'wheat',         # 7 baresoil
            #'lightgrey'      # 8 land ice
            'white'  # 8 land ice
        ])
        bounds = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5]
        ticklabels = ['Broadleaf','Needleleaf','C3 grass','C4 grass',
                      'Shrub','Lake','Soil','Ice']
    

    if type == 't':
        cmap = ListedColormap([
            'forestgreen',   # 1 broadleaf trees
            'darkgreen',     # 2 needleleaf trees
            'yellowgreen',    # 3 C3 grass
            'gold',          # 4 C4 grass
            'sienna',        # 5 shrub
            'blue',          # 6 lake
            'wheat',         # 7 baresoil
            #'lightgrey'      # 8 land ice
            'white',  # 8 land ice         
            'mediumseagreen',  # 9 mixed trees
            'olive',           # 10 mixed grass/shrub
            'lightgrey'        # 11 mixed
        ])
        bounds = [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5]
        ticklabels = ['Broadleaf','Needleleaf','C3 grass','C4 grass',
                      'Shrub','Lake','Soil','Ice','mixed trees',
                      'mixed grasses','mixed']
    
  
    norm = BoundaryNorm(bounds, cmap.N)

    lons = cube1.coord('longitude').points
    lats = cube1.coord('latitude').points

    fig, axes = plt.subplots(ncols=2, nrows=2,figsize=(12,8),
                             subplot_kw={'projection': ccrs.PlateCarree()},
                             constrained_layout=True)
    ax1=axes[0,0]
    ax2=axes[1,0]
    ax3=axes[1,1]
    ax4=axes[0,1]
       
  
    mesh = ax1.pcolormesh(lons, lats, biome1,
                         cmap=cmap,
                         norm=norm,
                         transform=ccrs.PlateCarree())

  #  ax1.pcolormesh(lons,lats,biome1.mask,cmap='Greys',alpha=0.2,
  #                 shading='nearest',edgecolors='k',
  #                 linewidths=0.2,transform=ccrs.PlateCarree())
    ax1.set_title(title1)
    ax1.set_facecolor('lightblue')

    
    mesh = ax2.pcolormesh(lons, lats, biome2,
                         cmap=cmap,
                         norm=norm,
                         transform=ccrs.PlateCarree())

    #ax2.coastlines()
    ax2.set_title(title2)
    ax2.set_facecolor('lightblue')


    mesh = ax3.pcolormesh(lons, lats, biome3,
                         cmap=cmap,
                         norm=norm,
                         transform=ccrs.PlateCarree())

    #ax2.coastlines()
    ax3.set_title(title3)
    ax3.set_facecolor('lightblue')

# --- Colourbar 
    cbar = fig.colorbar(mesh, ax=[ax1,ax2,ax3], orientation='horizontal',
                        fraction=0.05, ticks=np.arange(1,len(ticklabels)+1),
                        pad=0.08,aspect=40,shrink=1.0)
    cbar.ax.set_xticklabels(ticklabels,fontsize=12)

    plt.savefig(fileout)
#=================================================================
# MAIN PROGRAM STARTS HERE

threshold=0.6

filestart = '/home/earjcti/um/xqbws/database_averages/veg/'
static_file = filestart + '/P4_enh_mb_qrfrac.type.nc'
dynamic_file = filestart + '/xqbws_Monthly_Average_#pi_veg_pft_3901_4001.nc'

maxcube_p, threshcube_p = get_domveg_pi(threshold)
maxcube_d, threshcube_d = get_dominant_vegetation(dynamic_file,threshold,'d')
maxcube_s, threshcube_s = get_dominant_vegetation(static_file,threshold,'s')


plot_cube(threshcube_p,threshcube_s,threshcube_d,
          'PI static','LP static','LP dynamic',
          'above_'+str(threshold)+'.png','t')

plot_cube(maxcube_p,maxcube_s,maxcube_d,'PI static','LP static','LP dynamic',
          'most_common_veg.png','m')
