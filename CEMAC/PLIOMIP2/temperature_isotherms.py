#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:43:50 2019

@author: earjcti

This will plot the multimodel mean near surface air temperature with
isotherms (at a given level) marked out.  It is to understand which regions
are likely to support thermofrost

"""

import sys
import os
import iris
import numpy as np
import matplotlib.pyplot as plt

#os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid


def get_data(filereq, field, modeluse):
    """
    gets the field (field) from the file (filereq) and loads it
    into an iris cube (the model name is in modeluse)
    outputs a cube of the data that is as simple as possible
    """

    if modeluse == 'MMM':
        cube = iris.load_cube(filereq, field)
    else:
        cubes = iris.load(filereq)
        cube = cubes[0]
    cube.data = cube.data.astype('float32')

    if (modeluse == 'MIROC4m' or modeluse == 'COSMOS'):
        cube.units = 'Celsius'
    else:
        cube.convert_units('Celsius')

    for coord in cube.coords():
        name = coord.standard_name
        if name != 'latitude' and name != 'longitude':
            if name is None:
                if coord.long_name is None:
                    cube.remove_coord(coord.var_name)
                else:
                    cube.remove_coord(coord.long_name)
            else:
                cube.remove_coord(name)

    cube.cell_methods = None

    return cube

def plotmean_newaxis(cube, modelno_):
    """
    this will add a new axis to the cube which contains the model number
    this is needed for concatenation
    """

    tempcube_ = iris.util.new_axis(cube)
    tempcube_.add_dim_coord(iris.coords.DimCoord(modelno_,
                                                 standard_name='model_level_number',
                                                 long_name='model',
                                                 var_name='model',
                                                 units=None,
                                                 bounds=None,
                                                 coord_system=None,
                                                 circular=False), 0)
    return tempcube_

class Plotalldata:
    """
    This will plot the data from the timeperiod (ie mpwp or pi)
    """
    def __init__(self, cubelist, timeperiod_, fieldname, modelname,
                 outname):

        """
        inputs are:
                    modelname - like 'HadCM3' or 'MMM' multimodelmean
                    cube - a single cube containing the data
        """


        self.titlename = timeperiod_ + ' ' + modelname
        self.cubelist = cubelist
        self.outname = outname

        if fieldname == 'NearSurfaceTemperature':
            self.valmin = -30.
            self.valmax = 30.
            self.diff = 5.
            self.colormap = 'RdBu_r'
            self.use_cbar = 'y'
            
        if fieldname == 'MMM':
            self.valmin = 0.
            self.valmax = 30.
            self.diff = 5.
            self.colormap = 'gist_ncar_r'
            self.use_cbar = 'n'

     


    def plotdata(self):
        """
        this will plot all the cubes to a .eps or .png file
        input anom_cubes : a list of cubes containing the anomalies from the mean
        """

        cube = self.cubelist[0]
        latitudes = cube.coord('latitude').points
        lon = cube.coord('longitude').points

        plt.subplot(1,1,1)
        self.plotmap(lon, latitudes, 30., 180.)
        plt.savefig(self.outname + '_Asia.eps', bbox_inches='tight')
        plt.savefig(self.outname + '_Asia.png', bbox_inches='tight')
        plt.close()
        
        plt.subplot(1,1,1)
        self.plotmap(lon, latitudes, 180., 360.)
        plt.savefig(self.outname + '_America.eps', bbox_inches='tight')
        plt.savefig(self.outname + '_America.png', bbox_inches='tight')
        plt.close()
        
        plt.subplot(1,1,1)
        self.plotmap(lon, latitudes, 0., 360.)
        plt.savefig(self.outname + '_Globe.eps', bbox_inches='tight')
        plt.savefig(self.outname + '_Globe.png', bbox_inches='tight')
        plt.close()


        return

    def plotmap(self, longitudes, latitudes, left, right):
        """
        will plot the data in a map format

        """

        lons, lats = np.meshgrid(longitudes, latitudes)

        map = Basemap(llcrnrlon=left, urcrnrlon=right,
                      llcrnrlat=0.0, urcrnrlat=90.0,
                      projection='cyl', resolution='l')

        #map.drawmapboundary
        lonmap, latmap = map(lons, lats)
        map.drawcoastlines()

        for i, cube in enumerate(self.cubelist):
            cubedata = cube.data
            contourplot = map.contour(lonmap, latmap, cube.data, 
                                      ISOTHERM_NEEDED,
                                   colors=COLORUSE[i], linewidths=1.5,
                                   linestyles='solid')

        plt.title(self.titlename)
    



def main_model_ind():
    """
    get the data from the individual models
    plot the multimodel mean and the individual models on the same 
    figure
    """
    
    # get a mean value of all the cubes which are within +/- 2deg
    cubelist = iris.cube.CubeList([])
    for i, model in enumerate(MODELNAMES):
        filename = (FILESTART + model + '/' + exptname.get(timeperiod) +
                    '.' + FIELDNAME + '.allmean.nc')

        mod_cube = get_data(filename, FIELDNAME, model)
        mod_cube_2deg = get_within_pm2deg(mod_cube)

        if i == 0:
            tempcube_orig = mod_cube_2deg.copy(data=mod_cube_2deg.data)
        else:
            tempcube_orig.data = mod_cube_2deg.data
        tempcube = plotmean_newaxis(tempcube_orig, i)
        cubelist.append(tempcube)


    twodeg_cube = cubelist.concatenate_cube()
    any_in_range_cube = twodeg_cube.collapsed(['model_level_number'],
                                              iris.analysis.MEAN)


# now we want to plot in blue if multimodel mean is in the range and
# in red if only some of the models are in range

    data_mmm = permafrost_mean_cube.data
    data_any = any_in_range_cube.data
    ysize, xsize = np.shape(data_mmm)
    data_both = np.zeros(np.shape(data_mmm))
    for i in range(0, xsize):
        for j in range(0, ysize):
            if np.ma.is_masked(data_mmm[j, i]):
                data_both[j, i] = data_any[j, i]
            else:
                data_both[j, i] = -100.

    all_mean_cube = any_in_range_cube.copy(data=data_both)

    plioobj = Plotalldata(all_mean_cube, timeperiod, 'TEMP +/- 2DEG',
                          '+/- 2deg: any (red) MMM (blue)',
                          OUTSTART + timeperiod +'_within_2deg')
    plioobj.plotdata()
    
    
    
def main_time(timeperiod):
    """
    the main routine for a single timeperiod (likely mPWP or PI)
    toplot: y if we want to plot, n if we don't
    """

    exptname = {"pi" : "E280",
                "mPWP" : "EOI400"}

    filename = (FILESTART  + FIELDNAME + '_multimodelmean.nc')
    print(filename)

    mean_cube = get_data(filename, FIELDNAME + 'mean_' + timeperiod, 
                         'MMM')

# now we need to get all the models and see if any of the values
# are within +/-degC
    #if len(MODELNAMES) > 1:
    #    main_model_ind()


    return mean_cube

def main():
    """
    the main routine that will split the timeperiod up if appropriate
    if there are two timeperiods it will plot them on the same figure
    """
    
    cubelist_times = iris.cube.CubeList([])
    if len(TIMEPERIODS) ==1 : # just do main_timeperiod
        cube = main_time(TIMEPERIODS[0])
        cubelist_times.append(cube)
        
    if len(TIMEPERIODS) == 2:  # plotting two timeperiods on same figure
    
        cube1 = main_time(TIMEPERIODS[0])
        cubelist_times.append(cube1)
        cube2 = main_time(TIMEPERIODS[1])
        cubelist_times.append(cube2)
        
        data1 = cube1.data
        data2 = cube2.data
        
        
       
    # plot isotherm contour
        title = np.str(ISOTHERM_NEEDED) + 'deg: mPWP (red), PI (blue)'
        plotobj = Plotalldata(cubelist_times, '', 'MMM',
                          title,
                          OUTSTART + 'MMM_' + np.str(ISOTHERM_NEEDED) + 'deg')
        plotobj.plotdata()
       
        

# variable definition
LINUX_WIN = 'l'
ISOTHERM_NEEDED = 0.0

if LINUX_WIN == 'l':
    FILESTART = '/nfs/hera1/earjcti/regridded/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
    
FIELDNAME = 'NearSurfaceTemperature'
OUTSTART = FILESTART + 'allplots/' + FIELDNAME + '/'
UNITS = 'deg C'
#TIMEPERIODS = ['pi']
#TIMEPERIODS = ['mPWP']
TIMEPERIODS = ['pi','mPWP']
COLORUSE = ['blue','red']

#if len(TIMEPERIODS) > 1:
#    MODELNAMES = []
#else:
#    MODELNAMES = ['CESM1.0.5', 'COSMOS', 'EC-Earth3.1', 'GISS', 'HadCM3',
#                  'IPSLCM6A', 'IPSLCM5A2', 'IPSLCM5A',
#                  'MIROC4m', 'MRI-CGCM2.3',
#                  'NorESM-L', 'NorESM1-F',
#                  'UofT'
#                  ]

main()
