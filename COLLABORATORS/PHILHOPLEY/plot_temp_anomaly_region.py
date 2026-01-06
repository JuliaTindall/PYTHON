#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Monday June 3rd 2019

#@author: earjcti
#
# phil has asked if I can difference Xiaofangs orbital sensitivity
# experiments over africa
#
########################################################
# other notes are

import os
import numpy as np
import scipy as sp
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import netCDF4
#from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
#from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
from iris.cube import CubeList
import iris.quickplot as qplt
import cartopy as cart
import cartopy.crs as ccrs
#import cf_units as unit
#from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans
import sys

def make_cmap(colors, position=None, bit=False):
    '''
    I didn't write this I found it on the web.
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mp.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap


# functions start her

def customise_cmap2():
    """
    precipitation colormap
    """
    colors = [(84, 48, 5), (113, 70, 16), (143, 93, 27), (173, 115, 38),
              (195, 137, 60), (206, 160, 97), (216, 182, 135),
              (227, 204, 173), (238, 226, 211), (248, 248, 247),
              (212, 230, 229), (176, 212, 209), (140, 194, 190),
              (103, 176, 170), (67, 158, 150), (44, 135, 127),
              (29, 110, 100), (14, 85, 74), (0, 60, 48)]
    my_cmap = make_cmap(colors, bit=True)
    return my_cmap


def get_mean_annual(expt,field,extra,startyear,endyear):
    """
    gets the mean annual averaged over 50 years
    """

    if field == 'temp':
        long_field='TEMPERATURE AT 1.5M'
    if field == 'precip':
        long_field='TOTAL PRECIPITATION RATE     KG/M2/S'

    filestart = ('/nfs/hera1/earjcti/um/' + expt + '/netcdf/' 
                 + expt + 'a@pd' + extra)
    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']

    cubes = CubeList([])
    for year in range(startyear,endyear):
        print(year)
        for mon in monthnames:
            filename = filestart + str(year) + mon + '.nc'
            cube = iris.load_cube(filename,long_field)
            cubes.append(cube)
       
    iris.util.equalise_attributes(cubes)
    fieldcube = cubes.concatenate_cube()

    avgcube = fieldcube.collapsed('t',iris.analysis.MEAN)
    avgcube = iris.util.squeeze(avgcube)

    if field == 'precip':
        avgcube.data = avgcube.data * 60. * 60. * 24.
        avgcube.units='mm/day'

    return avgcube


def plotdata(cube,field):
    """
    plots the data
    """

    ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree(central_longitude=0.0))
    # set for Africa
    ax.set_extent([-30,60,-40,40],ccrs.PlateCarree())

    # get data to overplot
    site_lons = [24.6,33.9,35.75, 36.0]
    site_lats = [-27.6, -10.0, 4, 5]

    print(field)
    if field == 'temp':
        print('here')
        V=np.arange(-6,7,1)
        cbar='RdBu_r'
        titlename='Temperature Anomaly'
        cube.units='degC'
    if field == 'precip':
        V = np.arange(-3.0,3.5,0.5)
        cbar=customise_cmap2()
        titlename = 'Precipitation Anomaly'

    # plot map
    qplt.contourf(cube,levels=V,cmap=cbar,extend='both')
    # overplot sites
    plt.scatter(site_lons,site_lats,c='black',marker='x',s=20,
                transform=ccrs.Geodetic())
    plt.text(site_lons[0]+1.0,site_lats[0],'T')
    plt.text(site_lons[1]-7.0,site_lats[1],'KB')
    plt.text(site_lons[2]-7.0,site_lats[2]-2.0,'NF')
    plt.text(site_lons[3]+1.0,site_lats[3],'SF')



    # plot stuff to make it look nice
    plt.gca().coastlines()
    plt.title(titlename)
    plt.savefig(field + '.png')
    plt.show()
    plt.close()
    

    

#end def plotdata

##########################################################
# main program


expt='xogzc'  #xogzc-EmaxPminOmax max-sh_summer insolation
cntl='xogzf'  #xogzf-EmaxPmaxOmax max-nh_summer insolation
field='temp'
extra='w'
startyear=50
endyear=100

expt_cube = get_mean_annual(expt,field,extra,startyear,endyear)
cntl_cube = get_mean_annual(cntl,field,extra,startyear,endyear)


cube_anom=expt_cube - cntl_cube
plotdata(cube_anom,field)

sys.exit(0)
