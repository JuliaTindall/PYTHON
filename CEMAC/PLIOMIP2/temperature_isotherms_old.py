#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:43:50 2019

@author: earjcti

This will plot the multimodel mean near surface air temperature with
+/-2 degree isotherms marked out.  It is to understand which regions
are likely to support thermofrost

"""

import sys
import iris
from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.pyplot as plt


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

def get_within_pm2deg(cube):
    """
    This routine will mask out all the regions that are not within
    + or - 2degrees
    input: a temperature cube
    output: a new cube with all data not within +/-deg masked out.
    """

    cubedata = cube.data

    newdata = np.ma.masked_outside(cubedata, -2.0, 2.0, copy=True)

    newcube = cube.copy(data=newdata)

    return newcube

class Plotalldata:
    """
    This will plot the data from the timeperiod (ie mpwp or pi)
    """
    def __init__(self, cube, timeperiod_, fieldname, modelname,
                 outname):

        """
        inputs are:
                    modelname - like 'HadCM3' or 'MMM' multimodelmean
                    cube - a single cube containing the data
        """


        self.titlename = timeperiod_ + ' ' + modelname
        self.cube = cube
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

        if fieldname == 'TEMP +/- 2DEG':
            self.valmin = -10.
            self.valmax = -5.0
            self.diff = 1.0
            self.colormap = 'jet'
            self.use_cbar = 'n'



    def plotdata(self):
        """
        this will plot all the cubes to a .eps or .png file
        input anom_cubes : a list of cubes containing the anomalies from the mean
        """


        cubedata = self.cube.data
        latitudes = self.cube.coord('latitude').points
        lon = self.cube.coord('longitude').points
        datatoplot, longitudes = (shiftgrid(180., cubedata,
                                            lon, start=False))


        self.plotmap(datatoplot, longitudes, latitudes)


        return

    def plotmap(self, datatoplot, longitudes, latitudes):
        """
        will plot the data in a map format

        """

        lons, lats = np.meshgrid(longitudes, latitudes)

        map = Basemap(llcrnrlon=-180.0, urcrnrlon=180.0,
                      llcrnrlat=0.0, urcrnrlat=90.0,
                      projection='cyl', resolution='l')

        #map.drawmapboundary
        lonmap, latmap = map(lons, lats)
        map.drawcoastlines()

        values = np.arange(self.valmin, self.valmax, self.diff)
        contourplot = map.contourf(lonmap, latmap, datatoplot, values,
                                   cmap=self.colormap,
                                   extend='both')
        plt.title(self.titlename)

        if self.use_cbar == 'y':
            cbar = plt.colorbar(contourplot, orientation='horizontal')
            cbar.set_label(UNITS, size=10)
            cbar.ax.tick_params(labelsize=8, labelrotation=60)


        fileout = (self.outname + '.eps')
        plt.savefig(fileout, bbox_inches='tight')

        fileout = (self.outname + '.png')#

        plt.savefig(fileout, bbox_inches='tight')
        plt.close()

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

    filename = ('/nfs/hera1/earjcti/regridded/'
                + FIELDNAME + '_multimodelmean.nc')

    mean_cube = get_data(filename, FIELDNAME + 'mean_' + timeperiod, 
                         'MMM')

    permafrost_mean_cube = get_within_pm2deg(mean_cube)

# now we need to get all the models and see if any of the values
# are within +/-degC
    if len(MODELNAMES) > 1:
        main_model_ind()


    return permafrost_mean_cube

def main():
    """
    the main routine that will split the timeperiod up if appropriate
    if there are two timeperiods it will plot them on the same figure
    """
    
    if len(TIMEPERIODS) ==1 : # just do main_timeperiod
        dummy = main_time(TIMEPERIODS[0])
        
    if len(TIMEPERIODS) == 2:  # plotting two timeperiods on same figure
        cube1 = main_time(TIMEPERIODS[0])
        cube2 = main_time(TIMEPERIODS[1])
        
        data1 = cube1.data
        data2 = cube2.data
        ysize, xsize = np.shape(data1)
        databoth = np.zeros((ysize, xsize))
      
        for i in range(0, xsize):
            for j in range(0, ysize):
                if np.ma.is_masked(data1[j, i]):
                    pass
                else:
                    data1[j, i] = 10
                    databoth[j, i] = databoth[j, i] + data1[j, i]
                if np.ma.is_masked(data2[j, i]):
                    pass
                else:
                    data2[j, i] = 20
                    databoth[j, i] = databoth[j, i] + data2[j, i]
       
        databoth = np.where(databoth == 0, np.nan, databoth) 
        print(databoth[:,0])

        cubeboth = cube1.copy(data=databoth)
        
        plotobj = Plotalldata(cubeboth, '', 'MMM',
                              '+/- 2deg: mPWP (green), PI (orange), both (blue)',
                              OUTSTART + 'MMM_within_2deg')
        plotobj.plotdata()
        plt.show()
        
        

# variable definition
FILESTART = ('/nfs/hera1/earjcti/regridded/')
FIELDNAME = 'NearSurfaceTemperature'
OUTSTART = FILESTART + 'allplots/' + FIELDNAME + '/'
UNITS = 'deg C'
#TIMEPERIODS = ['pi']
#TIMEPERIODS = ['mPWP']
TIMEPERIODS = ['pi','mPWP']

if len(TIMEPERIODS) > 1:
    MODELNAMES = []
else:
    MODELNAMES = ['CESM1.0.5', 'COSMOS', 'EC-Earth3.1', 'GISS', 'HadCM3',
                  'IPSLCM6A', 'IPSLCM5A2', 'IPSLCM5A',
                  'MIROC4m', 'MRI-CGCM2.3',
                  'NorESM-L', 'NorESM1-F',
                  'UofT'
                  ]

main()
