#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.09.2019 by Julia

This program will plot a given field from the individual models 
for either the Pliocene or the preindustrail or the difference between them

It will subtract the multimodel mean so that the differences are very clear
"""

import os
import numpy as np
import scipy as sp
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from mpl_toolkits.basemap import Basemap, shiftgrid
#import Basemap
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
import cf_units as unit
from iris.experimental.equalise_cubes import equalise_attributes
import sys


def getmeanfield(fieldname, period):
    """
    get the mean values from the mean value file 
    
    inputs: fieldname (probably NearSurfaceTemperature or TotalPrecipitation)
            period (probably mPWP PI or anomaly)

    returns the mean value from the multuimodelmean.nc file
    """
    
    meanfile = ('/nfs/hera1/earjcti/regridded/' + fieldname + 
                '_multimodelmean.nc')
    meanfield = fieldname + 'mean_' + period
    
    cube = iris.load_cube(meanfile, meanfield)

    return cube


def getmodelfield(modelname, fieldname, period):
     """
     get the mean values from the model data
     inputs: modelname (ie HadCM3)
             fieldname (probably NearSurfaceTemperature or TotalPrecipitation)
             period (likely EOI400 or E280)
     returns:  a cube contatining the mean data from the model
     
     """
     
     modfile = ('/nfs/hera1/earjcti/regridded/' + modelname + '/' + 
                period + '.' + fieldname + '.allmean.nc')
    
     tempcube = iris.load(modfile)
     cube = tempcube[0]
     cube.units = UNITS
     
     return cube

class Plotalldata:   
    def __init(self):
        self.data = []
        
        
    def plotdata(self):
        """
        this will plot all the cubes to a .eps or .png file
        input anom_cubes : a list of cubes containing the anomalies from the mean
              modelnames : the names of all the models
        """
        self.nmodels = len(modelnames)
        self.filestart = ('/nfs/hera1/earjcti/regridded/allplots/' +
                          fieldname + '/individual')
    
        for i in range(0, self.nmodels):
            
            cubedata = anom_cubes[i].data
            self.latitudes = anom_cubes[i].coord('latitude').points
            lon = anom_cubes[i].coord('longitude').points
            self.datatoplot, self.longitudes = (shiftgrid(180.,
                                                          cubedata,
                                                          lon,
                                                          start=False))
            self.model = modelnames[i]
            self.plotmap(i)
        

        return      
    
    def plotmap(self, i):
        """
        will plot the data in a map format
        
        """
        
        plotpos = np.mod(i, 8) + 1
        plt.subplot(3, 3, plotpos)
        lons, lats = np.meshgrid(self.longitudes, self.latitudes)
        
        map=Basemap(llcrnrlon=-180.0, urcrnrlon=180.0, 
                    llcrnrlat=-90.0, urcrnrlat=90.0, 
                    projection='cyl',resolution='l')
   
        map.drawmapboundary
        x, y = map(lons, lats)
        map.drawcoastlines()
    
        V = np.arange(-5.0, 6.0, 1.0)    
        cs = map.contourf(x, y, self.datatoplot, V, cmap='RdBu_r',
                          extend='both')
        plt.title(self.model)
        
        print(i,self.nmodels)
        if plotpos == 8 or (i + 1) == self.nmodels:
            plt.subplot(3, 3, 9)
            plt.gca().set_visible(False)
            cbar = plt.colorbar(cs, orientation='horizontal')
            cbar.set_label(UNITS)
            
            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                        + '.eps')
            plt.savefig(fileout, bbox_inches='tight')
            
            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                    + '.png')
            
            plt.savefig(fileout, bbox_inches='tight')
        
        
##########################################################
# main program
# set up variable information

fieldname='NearSurfaceTemperature'
filename=' '
linux_win='l'


modelnames=['CESM1.0.5','COSMOS','EC-Earth3.1','GISS','HadCM3',
            'IPSLCM6A','IPSLCM5A2','IPSLCM5A',
            'MIROC4m','MRI-CGCM2.3',
            'NorESM-L','NorESM1-F',
            'UofT',
            ]

#modelnames=['GISS']
UNITS = 'Celsius'

#fieldnames=['TotalPrecipitation']
#units=['mm/day']

# set up cubelists to store data
mpwp_anom_cubes=iris.cube.CubeList([])
pi_anom_cubes=iris.cube.CubeList([])
anom_anom_cubes=iris.cube.CubeList([])


#################################################
# get mean data
mean_plio_cube = getmeanfield(fieldname, 'mPWP')
mean_pi_cube = getmeanfield(fieldname, 'pi')
mean_anom_cube = getmeanfield(fieldname, 'anomaly')


for model in range(0,len(modelnames)):
    model_plio_cube = getmodelfield(modelnames[model], fieldname, 'EOI400')
    model_pi_cube = getmodelfield(modelnames[model], fieldname, 'E280')
    model_anom_cube = model_plio_cube - model_pi_cube

    mpwp_anom_cubes.append(model_plio_cube - mean_plio_cube)
    pi_anom_cubes.append(model_pi_cube - mean_pi_cube)
    anom_anom_cubes.append(model_anom_cube - model_pi_cube)
    
##################################################
# plot the cubes for the model anomalies relative to the mean

obj = Plotalldata()
obj.plotdata(fieldname, mpwp_anom_cubes, modelnames)

#sys.exit(0)
