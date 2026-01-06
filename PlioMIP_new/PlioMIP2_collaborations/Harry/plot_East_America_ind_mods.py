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
import iris.quickplot as qplt
import iris.coord_categorisation
import cartopy
import cartopy.crs as ccrs

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
    monthnumber = {'January' : 1.0,
                   'August' : 8.0} 
   

    meanfile = (FILESTART + 'regridded/' + FIELDNAME +
                '_multimodelmean_month.nc')
    if period == 'plio - pi':
        meanfield = FIELDNAME + period
    else:
        meanfield = FIELDNAME + 'mean_' + period

    print(meanfile, meanfield)
    cube = iris.load_cube(meanfile, meanfield)
    if MONTH == 'Annual':
        cubeperiod = cube.collapsed('time', iris.analysis.MEAN)
    else:
        cubeperiod = cube.extract(iris.Constraint(time=monthnumber.get(MONTH)))
    

    return cubeperiod


def getmodelfield(modelname, period):
    """
    get the mean values from the model data
    inputs: modelname (ie HadCM3)
            period (likely EOI400 or E280)
    returns:  a cube contatining the mean data from the model

    """
    monthnumber = {'January' : 1.0,
                   'August' : 8.0} 

    print(modelname)
    if MONTH == 'Annual':
        modfile = (FILESTART + 'regridded/' + modelname + '/' +
                   period + '.' + FIELDNAME + '.allmean.nc')

        tempcube = iris.load(modfile)
        cube = tempcube[0]
    else:
        modfile = (FILESTART + 'regridded/' + modelname + '/' +
                   period + '.' + FIELDNAME + '.mean_month.nc')
        cube = iris.load_cube(modfile)
        if ((modelname == 'NorESM-L' or modelname == 'NorESM1-F') 
            and FIELDNAME == 'SST'):
            cubetemp = iris.load(modfile)
            cubetemp2=cubetemp[0]
            cube = cubetemp2[np.int(monthnumber.get(MONTH))-1, :, :]
        else:
            monthreq = iris.Constraint(month = monthnumber.get(MONTH))
            cube = iris.load_cube(modfile, monthreq)
       
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
        self.filestart = ('/nfs/see-fs-02_users/earjcti' + 
                          '/PYTHON/PLOTS/PLIOMIP2/Collaborators/HJD_OMZ/' +
                          FIELDNAME + '_' + timeperiod + '_' + MONTH)
        self.timeperiod = timeperiod
        self.anom_cubes = anom_cubes

        if timeperiod == 'mPWP-PI':
            if FIELDNAME == 'NearSurfaceTemperature':
                self.valmin = 0.
                self.valmax = 8.
                self.diff = 0.5
                self.colormap = 'Reds'

            if FIELDNAME == 'SST':
                self.valmin = 0.
                self.valmax = 4.
                self.diff = 0.25
               
                if MONTH == 'August':
                    self.valmax=6.5
                    self.diff=0.5
                self.colormap = 'Reds'

                
             
            if FIELDNAME == 'TotalPrecipitation':
                self.valmin = -2.
                self.valmax = 2.1
                self.diff = 0.2
                if MONTH == 'August':
                    self.valmin = -3.
                    self.valmax = 3.2
               
                self.colormap = 'RdBu'

        else:
             if ((FIELDNAME == 'NearSurfaceTemperature' or FIELDNAME == 'SST')
                 and MONTH == 'Annual'):
                self.valmin = 15.
                self.valmax = 30.
                self.diff = 1.
                self.colormap = 'Spectral_r'
                
             if (FIELDNAME == 'NearSurfaceTemperature'and MONTH == 'January'):
                self.valmin = 0.
                self.valmax = 25.
                self.diff = 1.
                self.colormap = 'Spectral_r'
            
             if (FIELDNAME == 'SST'and MONTH == 'January'):
                self.valmin = 15.
                self.valmax = 30.
                self.diff = 1.
                self.colormap = 'Spectral_r'
            
             if ((FIELDNAME == 'NearSurfaceTemperature' or FIELDNAME == 'SST')
                 and MONTH == 'August'):
                self.valmin = 24.
                self.valmax = 36.
                self.diff = 1.
                self.colormap = 'Spectral_r'

             if FIELDNAME == 'TotalPrecipitation':
                self.valmin = 0.
                self.valmax = 10.
                self.diff = 1.0
                self.colormap = 'winter'



    def plotdata(self):
        """
        this will plot all the cubes to a .eps or .png file
        input anom_cubes : a list of cubes containing the anomalies from the mean
        """


        fig = plt.figure(figsize=(8.5, 9.5))
        xplot=3
        yplot=3
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
            plotpos = np.mod(i, xplot * yplot) + 1
       
            self.plotmap(i, plotpos, title_,
                         datatoplot, longitudes, latitudes, fig, xplot, yplot)
            if plotpos == (xplot * yplot):
                # new figure
                 fig = plt.figure(figsize=(8.5, 9.5))
       

        return

    def plotmap(self, i, plotpos, titlename, datatoplot, 
                longitudes, latitudes, fig,
                xplot,yplot):
        """
        will plot the data in a map format

        """


        plt.subplot(xplot, yplot, plotpos)
        lons, lats = np.meshgrid(longitudes, latitudes)

        map = Basemap(llcrnrlon=-85.0, urcrnrlon=-70.0,
                      llcrnrlat=25.0, urcrnrlat=40.0,
                      projection='cyl', resolution='l')

        #map.drawmapboundary
        x, y = map(lons, lats)
        map.drawcoastlines(linewidth=0.5)

        V = np.arange(self.valmin, self.valmax, self.diff)
        cs = map.contourf(x, y, datatoplot, V, cmap=self.colormap,
                          extend='both')
        plt.title(titlename)

        print(i+1, self.nmodels)
        if plotpos == (xplot * yplot) or (i+1)  == self.nmodels:
             # Shrink current axis by 20% and put a legend to the right
            plt.subplots_adjust(left=0.05, bottom=0.1, right=0.82, top=0.9,
                                wspace=0.1, hspace=0.0)

            cb_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
           
            cbar = fig.colorbar(cs, cax=cb_ax, orientation='vertical')
            #cbar = plt.colorbar(fig, orientation='horizontal')
            #fig.colorbar(fix, ax=axs[:, col], shrink=0.6)
            cbar.ax.tick_params(labelsize=15)
            cbar.set_label(UNITS, fontsize=15)
            print('plotted colorbar', i)
            #plt.show()
            #plt.tight_layout()
            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                       + '.eps')
            plt.savefig(fileout, bbox_inches='tight')

            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                       + '.png')

            plt.savefig(fileout, bbox_inches='tight')
            plt.close()
       
def plot_mmm(mean_plio_cube, mean_pi_cube, mean_anom_cube):
    """
    plots the MMM to a file
    """

    vmin = {'mPWPNearSurfaceTemperatureJanuary' : 5.0, 
            'anomNearSurfaceTemperatureJanuary' : 1.0,
            'mPWPNearSurfaceTemperatureAugust' : 27.0, 
            'anomNearSurfaceTemperatureAugust' : 1.0,
            'mPWPNearSurfaceTemperatureAnnual' : 15.0, 
            'anomNearSurfaceTemperatureAnnual' : 2.0,
            'mPWPSSTJanuary' : 15.0, 'anomSSTJanuary' : 1.0,
            'mPWPSSTAugust' : 27.0, 'anomSSTAugust' : 1.0,
            'mPWPSSTAnnual' : 20.0, 'anomSSTAnnual' : 1.6,
            'mPWPTotalPrecipitationJanuary' : 0.0, 
            'anomTotalPrecipitationJanuary' : -1.0,
            'mPWPTotalPrecipitationAugust' : 0.0, 
            'anomTotalPrecipitationAugust' : -2.0,
            'mPWPTotalPrecipitationAnnual' : 0.0, 
            'anomTotalPrecipitationAnnual' : -1.0,
            
            }
    vmax = {'mPWPNearSurfaceTemperatureJanuary' : 26.0, 
            'anomNearSurfaceTemperatureJanuary' : 3.2, 
            'mPWPNearSurfaceTemperatureAugust' : 32.5, 
            'anomNearSurfaceTemperatureAugust' : 5.0,
            'mPWPNearSurfaceTemperatureAnnual' : 30.0, 
            'anomNearSurfaceTemperatureAnnual' : 5.0,
            'mPWPSSTJanuary' : 30.0, 'anomSSTJanuary' : 3.2, 
            'mPWPSSTAugust' : 32.5, 'anomSSTAugust' : 4.0,
            'mPWPSSTAnnual' : 30.0, 'anomSSTAnnual' : 3.8,
            'mPWPTotalPrecipitationJanuary' : 6.5, 
            'anomTotalPrecipitationJanuary' : 1.1, 
            'mPWPTotalPrecipitationAugust' : 6.5, 
            'anomTotalPrecipitationAugust' : 2.2,
            'mPWPTotalPrecipitationAnnual' : 6.5, 
            'anomTotalPrecipitationAnnual' : 1.1,
            }

    vdiff = {'mPWPNearSurfaceTemperatureJanuary' : 1.0, 
             'anomNearSurfaceTemperatureJanuary' : 0.2,
             'mPWPNearSurfaceTemperatureAugust' : 0.5, 
             'anomNearSurfaceTemperatureAugust' : 0.5,
             'mPWPNearSurfaceTemperatureAnnual' : 1.0, 
             'anomNearSurfaceTemperatureAnnual' : 0.25,
             'mPWPSSTJanuary' : 1.0, 'anomSSTJanuary' : 0.2,
             'mPWPSSTAugust' : 0.5, 'anomSSTAugust' : 0.5,
             'mPWPSSTAnnual' : 0.5, 'anomSSTAnnual' : 0.2,
             'mPWPTotalPrecipitationJanuary' : 0.5, 
             'anomTotalPrecipitationJanuary' : 0.1,
             'mPWPTotalPrecipitationAugust' : 0.5, 
             'anomTotalPrecipitationAugust' : 0.2,
             'mPWPTotalPrecipitationAnnual' : 0.5, 
             'anomTotalPrecipitationAnnual' : 0.1,
             
             }
    
    cmap = {'mPWPSST' : 'Spectral_r', 'anomSST' : 'Reds',
            'mPWPNearSurfaceTemperature' : 'Spectral_r', 
            'anomNearSurfaceTemperature' : 'Reds',
            'mPWPTotalPrecipitation' : 'winter', 
            'anomTotalPrecipitation' : 'RdBu'}

  
    # plot mPWP
    ax1 = plt.subplot(1, 2, 1, projection = ccrs.PlateCarree())
    ax1.coastlines()
    ax1.set_extent([-85., -70., 25., 40.], ccrs.PlateCarree())
    V = np.arange(vmin.get('mPWP'+FIELDNAME + MONTH), 
                  vmax.get('mPWP'+FIELDNAME + MONTH),
                  vdiff.get('mPWP'+FIELDNAME + MONTH))
    qplt.contourf(mean_plio_cube, V, cmap=cmap.get('mPWP'+FIELDNAME), 
                  extend='both')
    plt.title(MONTH + ' MMM ' + FIELDNAME + ':mPWP')
    
    # plt mPWP- PI
    
    ax2 = plt.subplot(1, 2, 2, projection = ccrs.PlateCarree())
    ax2.coastlines()
    ax2.set_extent([-85., -70., 25., 40.], ccrs.PlateCarree())
    V = np.arange(vmin.get('anom'+FIELDNAME + MONTH), 
                  vmax.get('anom'+FIELDNAME + MONTH),
                  vdiff.get('anom'+FIELDNAME + MONTH))  
    qplt.contourf(mean_anom_cube, V, cmap=cmap.get('anom'+FIELDNAME),
                  extend='both')
    plt.title('mPWP - PI')
   

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'Collaborators/HJD_OMZ/' + FIELDNAME + '_MMM_' + MONTH + '.png')
    plt.savefig(fileout)
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/' + 
               'Collaborators/HJD_OMZ/' + FIELDNAME + '_MMM_' + MONTH + '.eps')
    plt.savefig(fileout)
    plt.show()
    plt.close()
   

##########################################################
# main program
# set up variable information

#FIELDNAME = 'NearSurfaceTemperature'
#UNITS = 'Celsius'
FIELDNAME = 'SST'
UNITS = 'Celsius'
FIELDNAME = 'TotalPrecipitation'
UNITS = 'mm/day'
#FIELDNAME = 'SST'
LINUX_WIN = 'l'
#MONTH = 'January'
MONTH = 'August'
#MONTH = 'Annual'

if LINUX_WIN == 'l':
    FILESTART = '/nfs/hera1/earjcti/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'

MODELNAMES = ['CESM2', 'HadGEM3', 'IPSLCM6A', 'COSMOS', 
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
mean_plio_cube = getmeanfield('plio')
mean_pi_cube = getmeanfield('pi')
mean_anom_cube = getmeanfield('plio - pi')

plot_mmm(mean_plio_cube, mean_pi_cube, mean_anom_cube)

for model, modelname in enumerate(MODELNAMES):
    model_plio_cube = getmodelfield(modelname, 'EOI400')
    model_pi_cube = getmodelfield(modelname, 'E280')
    if modelname == 'EC-Earth3.1' and FIELDNAME == 'SST':
       model_pi_cube.coord('latitude').bounds = None
       model_pi_cube.coord('longitude').bounds = None

    model_anom_cube = model_plio_cube - model_pi_cube

    mpwp_anom_cubes.append(model_plio_cube)
    pi_anom_cubes.append(model_pi_cube)
    anom_anom_cubes.append(model_anom_cube)

##################################################
# plot the cubes for the model anomalies relative to the mean

obj = Plotalldata('mPWP', mpwp_anom_cubes)
obj.plotdata()

obj = Plotalldata('PI', pi_anom_cubes)
obj.plotdata()

obj = Plotalldata('mPWP-PI', anom_anom_cubes)
obj.plotdata()

#
