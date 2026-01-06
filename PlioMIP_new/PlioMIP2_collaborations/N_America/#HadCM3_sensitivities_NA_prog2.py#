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

def simplify_cube(cube):
    """
    gets cube and makes sure dimensions are longitude, latitude surface and
    t. 
    """    
    for coord in cube.coords():
        if coord.var_name == 'level275':
            coord.var_name = 'surface'
    
    cube.coord('surface').points = 0.0
    cube.coord('surface').units = 'm'
    cube.coord('surface').attributes = None
    
    cube.data = np.where(cube.data > 1.0E10, 0., cube.data)
    return cube


def plotmap(titlename, datatoplot, longitudes, latitudes):
        """
        will plot the data in a map format

        """

      
        lons, lats = np.meshgrid(longitudes, latitudes)

      #  map = Basemap(llcrnrlon=-130.0, urcrnrlon=-105.0,
      #                llcrnrlat=30.0, urcrnrlat=45.0,
      #                projection='cyl', resolution='l')

        map = Basemap(llcrnrlon=-150.0, urcrnrlon=-80.0,
                      llcrnrlat=15.0, urcrnrlat=55.0,
                      projection='cyl', resolution='l')

        #map.drawmapboundary
        x, y = map(lons, lats)
        map.drawcoastlines(linewidth=1.0)
        map.drawstates(linewidth=0.5)

        if FIELDNAME == 'TotalPrecipitation':
           V = np.arange(-1., 1.1, 0.1)
           cmapname = 'RdBu'
        if FIELDNAME == 'NSAT':
           V = np.arange(0,8.5,0.5)
           cmapname = 'Reds'
        if FIELDNAME == 'pmine':
           V = np.arange(-1,1.1,0.1)
           cmapname = 'RdBu'
        cs = map.contourf(x, y, datatoplot, V, cmap=cmapname,
                          extend='both')
        plt.title(titlename)
       # x, y = map(-117, 36)
       # plt.plot(x, y, 'ok', markersize=5)


       # if plotpos == (xplot * yplot) or (i + 1) == len(expts):
       #      # Shrink current axis by 20% and put a legend to the right
       #     plt.subplots_adjust(left=0.05, bottom=0.1, right=0.82, top=0.9,
       #                         wspace=0.1, hspace=0.0)#
       #
       # cb_ax = plt.add_axes([0.85, 0.15, 0.02, 0.7])
        cbar = plt.colorbar(cs, orientation='horizontal')
       #     #fig.colorbar(fix, ax=axs[:, col], shrink=0.6)
        cbar.ax.tick_params(labelsize=15)
        cbar.set_label(UNITS, fontsize=15)
       #     print('plotted colorbar')
       #plt.show()
       #plt.tight_layout()
        fileout = (FILEOUTSTART+'.eps')
        plt.savefig(fileout, bbox_inches='tight')

        fileout = (FILEOUTSTART + '.png')
        plt.savefig(fileout, bbox_inches='tight')
        plt.close()

def get_pmine(filename_):
    """
    will add up all the fluxes that make evaporation and returns
    p-e within a cube
    The fluxes are:
      evaporation from canopy
      evaporation from sea
      transpiration (this is exactly the same as evaporation from soil
                     I have checked some examples)
      sublim from surface
    """

    varnames_sec = ["EVAPORATION FROM SEA (GBM)   KG/M2/S",
                "TRANSPIRATION RATE           KG/M2/S"]
    
    varnames_ts = ["EVAP FROM CANOPY - AMOUNT   KG/M2/TS",
                "SUBLIM. FROM SURFACE (GBM)  KG/M2/TS"]
    
    precipcube = iris.load_cube(filename_,
				"TOTAL PRECIPITATION RATE     KG/M2/S")
    precipcube = simplify_cube(precipcube)
    precipcube.units = 'kg m-2 s-1'

    for i, var in enumerate(varnames_sec):
        cube = iris.load_cube(filename_,var)
        cube = simplify_cube(cube)
        if i == 0:
            cubetot = cube
        else:
            cubetot = cubetot + cube
        
        
    for i, var in enumerate(varnames_ts):
        cube = iris.load_cube(filename_,var)
        cube = simplify_cube(cube)
        cube.data = cube.data / (30. * 60.)
        cube.units = 'kg m-2 s-1'
        cubetot = cubetot + cube
    
    pminecube = precipcube - cubetot

    return pminecube


def get_mean(exptname, field):
    """
    get's the mean of the field of interest from the experiment (exptname)
    """
    filestart = ('/nfs/hera1/earjcti/um/' + exptname + '/pd/' +
                 exptname + 'a@pd' + EXTRA.get(exptname))

    print('getting mean',exptname)
    ncubes = 0
    cubelist = iris.cube.CubeList([])

    for year in range (STARTYEAR, ENDYEAR):
        print(year)
        for i,month in enumerate(MONTHNAMES):
            filename = filestart + np.str(year).zfill(2) + month + '.nc'
        #    print(filename,FIELDNAME)
            if FIELDNAME == 'pmine':
               cubet = get_pmine(filename)
            else:
                cubet = iris.load_cube(filename, field)
	    
            cubet.coord('t').attributes = None
            cube = iris.util.squeeze(cubet)
            cubelist.append(cube)
            ncubes = ncubes + 1

#    for cube in cubelist:
#	    print(cube)
#	    print(cube.coord('t'))
#	    print(cube.coord('latitude'))
    equalise_attributes(cubelist)
    iris.util.unify_time_units(cubelist)
    allcubes = cubelist.merge_cube()
    mean_cube = result = allcubes.collapsed('t', iris.analysis.MEAN)
    if FIELDNAME == 'TotalPrecipitation' or FIELDNAME == 'pmine':
	    mean_cube = mean_cube * 60. * 60. * 24.
	    mean_cube.units = 'mm/day'
   
    return mean_cube
##########################################################
# main program
# set up variable information

desc = {'tenvo':'E280', 'tenvj':'EOI400','tenvq': 'E400',
	'xozza':'PI', 'xozzc' : 'K1_(3.060Ma)', 'xozzd' : 'G17_(2.950Ma)',
	'xozze' : 'Km3_(3.155Ma)','xozzf':'3.053Ma',
	'tenvk':'EOI350','tenvl':'EOI450','xoorb':'EOI400_dynveg','xoora':'PI'}
EXTRA = {'tenvo' : 't', 'tenvj': 'o','tenvq' : 't',
	 'xozza':'o', 'xozzc':'o',  'xozzd':'o','xozze':'o','xozzf':'o',
	 'tenvk': 'o','tenvl' :'o','xoorb' : 't'}
longfield = {'TotalPrecipitation' : "TOTAL PRECIPITATION RATE     KG/M2/S",
	     'NSAT' : 'TEMPERATURE AT 1.5M'} 
fieldunits = {"TotalPrecipitation" : "mm/day",
	      "NSAT" : "degC","pmine" : "mm/day"} 
expt= 'xoorb'  
cntl = 'tenvo'
#expts = ['xozzb','xozzh','xozzi']

FIELDNAME = 'pmine'
#FIELDNAME = 'TotalPrecipitation'
#FIELDNAME = 'NSAT'

SENSNAME = desc.get(expt) + '-' + desc.get(cntl)
FILEOUTSTART = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
                'Collaborators/N_America/HadCM3_sensitivities/'
		+ FIELDNAME + '_HadCM3_' + SENSNAME)
MONTHNAMES = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
STARTYEAR=0
ENDYEAR=100

UNITS = fieldunits.get(FIELDNAME)
LINUX_WIN = 'l'

cntlmean = get_mean(cntl,longfield.get(FIELDNAME))
meandata = get_mean(expt, longfield.get(FIELDNAME))
anomcube = meandata - cntlmean
print(anomcube.data)
datatoplot, longitudes = (shiftgrid(180., anomcube.data,
                                      anomcube.coord('longitude').points,
                                      start=False))
            
plotmap(SENSNAME + ' ' + FIELDNAME, datatoplot, longitudes,
           anomcube.coord('latitude').points)
