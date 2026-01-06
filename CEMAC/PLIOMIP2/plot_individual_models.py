#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.09.2019 by Julia

This program will plot a given field from the individual models
for either the Pliocene or the preindustrail or the difference between them

It will subtract the multimodel mean so that the differences are very clear
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
import iris.coord_categorisation

#os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid


def getmeanfield(period):
    """
    get the mean values from the mean value file

    inputs: fieldname (probably NearSurfaceTemperature or TotalPrecipitation)
            period (probably mPWP PI or anomaly)

    returns the mean value from the multuimodelmean.nc file
    """

    meanfile = (FILESTART + 'regridded/' + FIELDNAME +
                '_multimodelmean.nc')
    meanfield = FIELDNAME + 'mean_' + period

    cube = iris.load_cube(meanfile, meanfield)


    return cube


def getmodelfield(modelname, period):
    """
    get the mean values from the model data
    inputs: modelname (ie HadCM3)
            period (likely EOI400 or E280)
    returns:  a cube contatining the mean data from the model

    """

    modfile = (FILESTART + 'regridded/' + modelname + '/' +
               period + '.' + FIELDNAME + '.allmean.nc')

    tempcube = iris.load(modfile)
    cube = tempcube[0]
    cube.units = UNITS

    #this will make all the dimensions of all the cubes match.


    for coord in cube.coords():
        name = coord.standard_name
        if name not in ['latitude', 'longitude']:
            if name is None:
                if coord.long_name is None:
                    cube.remove_coord(coord.var_name)
                else:
                    cube.remove_coord(coord.long_name)
            else:
                cube.remove_coord(name)

    for coord in cube.coords():   # now this will be longitude or latitude
        coord.points = coord.points.astype('float32')
        coord.var_name = coord.standard_name
        coord.long_name = coord.standard_name

    return cube


class Plotalldata:
    """
    This will plot the data from the timeperiod (ie mpwp or pi)
    """
    def __init__(self, timeperiod, anom_cubes):
        self.nmodels = len(MODELNAMES)
        self.filestart = (FILESTART + '/regridded/allplots/' +
                          FIELDNAME + '/' + timeperiod + '_individual')
        self.timeperiod = timeperiod
        self.anom_cubes = anom_cubes

        if (FIELDNAME == 'NearSurfaceTemperature'
            or FIELDNAME == 'SST'):
                self.valmin = -5.
                self.valmax = 6.
                self.diff = 1.
                self.colormap = 'RdBu_r'

        if FIELDNAME == 'TotalPrecipitation':
            self.valmin = -2.
            self.valmax = 2.1
            self.diff = 0.2
            self.colormap = 'RdBu'


    def plotdata(self):
        """
        this will plot all the cubes to a .eps or .png file
        input anom_cubes : a list of cubes containing the anomalies from the mean
        """


        fig = plt.figure(figsize=(11.0, 8.5))
        for i in range(0, self.nmodels):

            cubedata = self.anom_cubes[i].data
            latitudes = self.anom_cubes[i].coord('latitude').points
            lon = self.anom_cubes[i].coord('longitude').points
            datatoplot, longitudes = (shiftgrid(180., cubedata,
                                                lon, start=False))
            #if (np.mod(i, 8) + 1) == 1:
            #    title_ = (MODELNAMES[i] + ':' +
            #              self.timeperiod + ' (model - MMM)')
            #else:
            #    title_ = (MODELNAMES[i])

            title_ = (MODELNAMES[i])
            self.plotmap(i, title_,
                         datatoplot, longitudes, latitudes, fig)


        return

    def plotmap(self, i, titlename, datatoplot, longitudes, latitudes, fig):
        """
        will plot the data in a map format

        """

        xplot = 4
        yplot = 4


        plotpos = np.mod(i, xplot * yplot) + 1
        plt.subplot(xplot, yplot, plotpos)
        lons, lats = np.meshgrid(longitudes, latitudes)

        map = Basemap(llcrnrlon=-180.0, urcrnrlon=180.0,
                      llcrnrlat=-90.0, urcrnrlat=90.0,
                      projection='cyl', resolution='l')

        #map.drawmapboundary
        x, y = map(lons, lats)
        map.drawcoastlines(linewidth=0.5)

        V = np.arange(self.valmin, self.valmax, self.diff)
        cs = map.contourf(x, y, datatoplot, V, cmap=self.colormap,
                          extend='both')
        plt.title(titlename)


        if plotpos == (xplot * yplot) or (i + 1) == self.nmodels:
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
            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                       + '.eps')
            plt.savefig(fileout, bbox_inches='tight')

            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                       + '.pdf')

            plt.savefig(fileout, bbox_inches='tight')
            plt.close()


##########################################################
# main program
# set up variable information

#FIELDNAME = 'NearSurfaceTemperature'
#UNITS = 'Celsius'
FIELDNAME = 'SST'
UNITS = 'Celsius'
#FIELDNAME = 'TotalPrecipitation'
#UNITS = 'mm/day'
#FIELDNAME = 'SST'
LINUX_WIN = 'l'

if LINUX_WIN == 'l':
    FILESTART = '/nfs/hera1/earjcti/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'

MODELNAMES = ['CESM2', 'IPSLCM6A', 'COSMOS', 
            'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
            'MIROC4m', 'IPSLCM5A2', 'HadCM3',
            'GISS2.1G', 'CCSM4', 
            'CCSM4-Utr', 'CCSM4-UoT', 
            'NorESM-L', 'MRI2.3', 'NorESM1-F'
             ]

#MODELNAMES = ['NorESM-L']



# set up cubelists to store data
mpwp_anom_cubes = iris.cube.CubeList([])
pi_anom_cubes = iris.cube.CubeList([])
anom_anom_cubes = iris.cube.CubeList([])


#################################################
# get mean data
mean_plio_cube = getmeanfield('mPWP')
mean_pi_cube = getmeanfield('pi')
mean_anom_cube = getmeanfield('anomaly')


for model, modelname in enumerate(MODELNAMES):
    model_plio_cube = getmodelfield(modelname, 'EOI400')
    model_pi_cube = getmodelfield(modelname, 'E280')
    print(modelname)
    if modelname == 'EC-Earth3.1' and FIELDNAME == 'SST':
       model_pi_cube.coord('latitude').bounds = None
       model_pi_cube.coord('longitude').bounds = None

    model_anom_cube = model_plio_cube - model_pi_cube

    mpwp_anom_cubes.append(model_plio_cube - mean_plio_cube)
    pi_anom_cubes.append(model_pi_cube - mean_pi_cube)
    anom_anom_cubes.append(model_anom_cube - mean_anom_cube)

##################################################
# plot the cubes for the model anomalies relative to the mean

obj = Plotalldata('mPWP', mpwp_anom_cubes)
obj.plotdata()

obj = Plotalldata('PI', pi_anom_cubes)
obj.plotdata()

obj = Plotalldata('mPWP-PI', anom_anom_cubes)
obj.plotdata()

#
