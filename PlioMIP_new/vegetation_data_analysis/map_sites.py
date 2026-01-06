#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on December 2020
# 
#  This program will plot a map showing the amplitude of the seasonal cycle 
#  over Land for the pliocene and the preindustrial
#
#
#import os
import numpy as np
import pandas as pd
#import scipy as sp
#import cf
import iris
#import iris.util
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import netCDF4
from mpl_toolkits.basemap import Basemap, shiftgrid
#from netCDF4 import Dataset, MFDataset
#import iris.analysis.cartography
#import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
#import cf_units as unit
#from iris.experimental.equalise_cubes import equalise_attributes
import cartopy
import cartopy.crs as ccrs
#import matplotlib.ticker as mticker
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from mpl_toolkits.basemap import Basemap

import sys

 



    

  
def main():
    """
    plot map with all sites on.
    """

    sitedata = []
    # site data is
    # sitename, sitelat, sitelon, 
    sitedata.append(['MI', 77.5, 261])
    sitedata.append(['BP', 79, 278])
    sitedata.append(['LE', 67, 172])
    sitedata.append(['LCM', 64, -142])
    sitedata.append(['Lake Baikal', 56, 108])
             

    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_global()
    for data in sitedata:
        plt.plot(data[2],data[1], color='red',marker='o', 
                    linewidth=2, transform = ccrs.Geodetic())
  
    
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/map_sites.eps')
   
    plt.savefig(fileout)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'vegetation/map_sites.png')
   
    plt.savefig(fileout)
    plt.close()
  

##########################################################
# main program

main()
