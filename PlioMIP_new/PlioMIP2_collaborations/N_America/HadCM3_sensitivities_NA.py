#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 05.11.2021 by Julia

This program will plot the HadCM3 sensitivities over the N. America region
"""

import os
import sys
import numpy as np
#import matplotlib as mp
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#import netCDF4
#from netCDF4 import Dataset, MFDataset
import iris
import iris.analysis.cartography
import iris.quickplot as qplt
import iris.coord_categorisation
from iris.experimental.equalise_cubes import equalise_attributes

#os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid
import warnings
warnings.filterwarnings("ignore")



def plotmap(i, titlename, datatoplot, longitudes, latitudes, fig):
        """
        will plot the data in a map format

        """

        xplot = 2
        yplot = 2


        plotpos = np.mod(i, xplot * yplot) + 1
        plt.subplot(xplot, yplot, plotpos)
        lons, lats = np.meshgrid(longitudes, latitudes)

      #  map = Basemap(llcrnrlon=-130.0, urcrnrlon=-105.0,
      #                llcrnrlat=30.0, urcrnrlat=45.0,
      #                projection='cyl', resolution='l')

        map = Basemap(llcrnrlon=-170.0, urcrnrlon=-60.0,
                      llcrnrlat=10.0, urcrnrlat=70.0,
                      projection='cyl', resolution='l')

        #map.drawmapboundary
        x, y = map(lons, lats)
        map.drawcoastlines(linewidth=1.0)
        map.drawstates(linewidth=0.5)

        V = np.arange(-1., 1.1, 0.1)
        cs = map.contourf(x, y, datatoplot, V, cmap='RdBu',
                          extend='both')
        plt.title(titlename)
       # x, y = map(-117, 36)
       # plt.plot(x, y, 'ok', markersize=5)


        if plotpos == (xplot * yplot) or (i + 1) == len(expts):
             # Shrink current axis by 20% and put a legend to the right
            plt.subplots_adjust(left=0.05, bottom=0.1, right=0.82, top=0.9,
                                wspace=0.1, hspace=0.0)

            cb_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
           
            cbar = fig.colorbar(cs, cax=cb_ax, orientation='vertical')
            #cbar = plt.colorbar(fig, orientation='horizontal')
            #fig.colorbar(fix, ax=axs[:, col], shrink=0.6)
            cbar.ax.tick_params(labelsize=15)
            cbar.set_label(UNITS, fontsize=15)
            print('plotted colorbar')
            #plt.show()
            #plt.tight_layout()
            fileout = (FILEOUTSTART + np.str(np.int(np.ceil(i/8)))
                       + '.eps')
            print(fileout)
            plt.savefig(fileout, bbox_inches='tight')

            fileout = (FILEOUTSTART + np.str(np.int(np.ceil(i/8)))
                       + '.pdf')

            plt.savefig(fileout, bbox_inches='tight')
            plt.close()

def get_mean(exptname):
    """
    get's the mean of the field of interest from the experiment (exptname)
    """
    print(exptname)
    filestart = ('/nfs/hera1/earjcti/um/' + exptname + '/netcdf/' +
                 exptname + 'a@pd' + EXTRA)

    ncubes = 0
    cubelist = iris.cube.CubeList([])

    for year in range (STARTYEAR, ENDYEAR):
        for month in MONTHNAMES:
            filename = filestart + np.str(year) + month + '.nc'
            cubet = iris.load_cube(filename, FIELDNAME)
            cube = iris.util.squeeze(cubet)
            cubelist.append(cube)
            ncubes = ncubes + 1

    equalise_attributes(cubelist)
    allcubes = cubelist.merge_cube()
    mean_cube = result = allcubes.collapsed('t', iris.analysis.MEAN)
    mean_cube = mean_cube * 60. * 60. * 24.
    mean_cube.units = 'mm/day'
   
    return mean_cube
##########################################################
# main program
# set up variable information

FIELDNAME = 'TotalPrecipitation'
FILEOUTSTART = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
                'Collaborators/N_America/' + FIELDNAME + '_HadCM3')
MONTHNAMES = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
EXTRA='k'
STARTYEAR=40
ENDYEAR=60

cntl = 'xozza'
expts = ['xozzb','xozzh','xozzi']
UNITS = 'mm/day'
LINUX_WIN = 'l'
FIELDNAME = 'TOTAL PRECIPITATION RATE     KG/M2/S'
title = {'xozzb' : 'EOI400', 'xozzh' : 'P1 N_America', 'xozzi': 'PI N_america'}

cntlmean = get_mean(cntl)
fig = plt.figure(figsize=(6.0, 4.0))
       
for i, expt in enumerate(expts):
   meandata = get_mean(expt)
   anomcube = meandata - cntlmean
   datatoplot, longitudes = (shiftgrid(180., anomcube.data,
                                       anomcube.coord('longitude').points,
                                       start=False))
            
   plotmap(i, title.get(expt), datatoplot, longitudes,
           anomcube.coord('latitude').points, fig)
