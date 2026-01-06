#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

#This program will plot the fields for Jonathan Holmes from the 
#averaged data (data_6ka - data_11ka)
#
#We will do a map plot in the vicinity of Ireland
#We will plot the absolute value at 6ka, and an anomaly for all the other
#slices
#
#The fields we intend to plot are JJA:
#    d18op, SAT, precip amount, circulation patterns, SST, d18osw
#Created on Sat Mar  7 13:58:35 2020
#
#@author: julia

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import iris.quickplot as qplt

#os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
#os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid

def get_allcubes(field):
    """
    gets the data from all the timeslices and puts them in a list of cubes
    returns the list of cubes
    """
    allcubes = iris.cube.CubeList([])
    
    for i, slice in enumerate(SLICES):
        cube = get_data(slice, field)
        if i == 0:
            # we will plot the raw data
            allcubes.append(cube)
        else:
            # we will plot the anomaly from the first timeslice
            anom_data = cube.data - allcubes[0].data
            newcube = cube.copy(data=anom_data)
            allcubes.append(newcube)

    return allcubes

def get_abs_plotvals():
    """
    sets values for plot range for if we are doing a non_anomaly plot
    """
    
    if FIELDREQ == 'd18o':
        vmin = -10.0
        vmax = 1.0
        vdiff =1.0 
        
    if FIELDREQ == 'TEMPERATURE AT 1.5M':
        vmin = 0.0
        vmax = 35.0
        vdiff = 5.0
    print(FIELDREQ)
   
    print(FIELDREQ)
    if (FIELDREQ == 'TOTAL PRECIPITATION RATE    MM/DAY'
    or FIELDREQ == 'evap'):
        vmin=0.0
        vmax=3.0
        vdiff=0.1
   
        
    return vmin, vmax, vdiff

def get_anom_plotvals():
    """
    sets values for plot range for if we are doing ananomaly plot
    """
    if FIELDREQ == 'd18o':
        vmin = -10
        vmax = 11
        vdiff = 1

    if FIELDREQ == 'TEMPERATURE AT 1.5M':
        vmin = -3.0
        vmax = 3.2
        vdiff = 0.2
        
    if (FIELDREQ == 'TOTAL PRECIPITATION RATE    MM/DAY'
    or FIELDREQ == 'evap'):
        vmin=-10
        vmax=11
        vdiff=1.0
        
    return vmin, vmax, vdiff

def jja_weight_d18o_by_precip_amt(d18ocube, filename):
    """
    if we are doing d18o we need to weight by precipitation amount to calculate
    the three month average
    """
    # get precipitation amount
    allprecipcube = iris.load_cube(filename, 'TOTAL PRECIPITATION RATE    MM/DAY')
  
    # get jja fields
    precip_jja_cube = allprecipcube[5:8, :, :]
    d18o_jja_cube = d18ocube[5:8, :, :]
    
    # weight d18o by precipitaion amount and average
    newdata = ((precip_jja_cube.data * d18o_jja_cube.data) / 
               np.mean(precip_jja_cube.data, axis=0))
    newcube = precip_jja_cube.copy(data=newdata)
    weighted_cube = newcube.collapsed('t', iris.analysis.MEAN)
    
    
    return weighted_cube

def jja_avg(cube, filename):
    """
    simple 3 month average over jja of the cube
    """
    #
    # get jja fields
    jja_cube = cube[5:8, :, :]
    
    # weight d18o by precipitaion amount and average
    avg_cube = jja_cube.collapsed('t', iris.analysis.MEAN)
    
    if FIELDREQ == 'TEMPERATURE AT 1.5M':
        avg_cube.data = avg_cube.data - 273.15
    
    return avg_cube

def sum_evap(filename_):
    """
    will add up all the fluxes that make evaporation and returns
    total evaporation within a cube
    The fluxes are:
      evaporation from canopy
      evaporation from sea
      transpiration
      sublim from surface
    """

    varnames_mm = ["EVAPORATION FROM SEA    MM/DAY",
                "TRANSPIRATION RATE    MM/DAY",
                "EVAP FROM CANPOPY    MM/DAY",
                ]
    
    varnames_ts = ["SUBLIM. FROM SURFACE (GBM)  KG/M2/TS"]
    
    
    for i, var in enumerate(varnames_mm):
        print(filename_, var)
        cube = iris.load_cube(filename_,var)
        cube.data = np.where(np.isnan(cube.data), 0, cube.data)
        cube.data = np.where(cube.data > 1.0E10, 0, cube.data)
        #cube = simplify_cube(cube)
        if i == 0:
            cubetot = cube
        else:
            cubetot = cubetot + cube
        print(var,cube.data)
            
        
    for i, var in enumerate(varnames_ts):
        cube = iris.load_cube(filename_,var)
        cube.data = np.where(np.isnan(cube.data), 0, cube.data)
        cube.data = np.where(cube.data > 1.0E10, 0, cube.data)
        #cube = simplify_cube(cube)
        cube.data = cube.data * 48.0
        cube.units = 'mm/day'
        cubetot = cubetot + cube
     
   
    return cubetot

def get_data(timeslice, field):
    """
    will get the data for the given timeslice and field
    and return in a cube
    """
    filename = FILESTART + timeslice + '.nc'
    print(filename, field)
    if field == 'evap':
        fieldcube = sum_evap(filename)
    else:
        print(filename)
        print(field)
        fieldcube = iris.load_cube(filename, field)
    
    if field == 'd18o':
        newcube = jja_weight_d18o_by_precip_amt(fieldcube, filename)
    else:
        newcube = jja_avg(fieldcube, filename)
        
    return newcube
        
def plot_data(cubelist):
    """
    plots all the cubes : one from each timeslice
    """
  
    for i, cube in enumerate(cubelist):
        print(i)
        plt.subplot(2, 3, i+1)
    
        map = Basemap(llcrnrlon=-20.0, urcrnrlon=20.0,
                      llcrnrlat=45.0, urcrnrlat=70.0,
                      projection='cyl', resolution='l')

        latitudes = cube.coord('latitude').points
        longitudes = cube.coord('longitude').points
        
        data_shift, lons_shift = shiftgrid(
                180, cube.data, longitudes, start=False)
        
        lons, lats = np.meshgrid(lons_shift, latitudes)
        x, y = map(lons, lats)
        map.drawcoastlines(linewidth=0.5)
        
        plt.rcParams['text.latex.preamble']=[r"\usepackage{wasysym}"]
        if i==0:
            cmap_j = 'rainbow'
            vmin, vmax, vdiff = get_abs_plotvals()
            titlename = SLICES[i]
            cs = map.contourf(x, y, data_shift, 
                          levels =np.arange(vmin,vmax,vdiff),
                          cmap = cmap_j, extend='both')
            plt.title(titlename,fontsize=10)
            cbar = plt.colorbar(cs, orientation='horizontal')
            cbar.ax.tick_params(labelsize='small')
            cbar.set_label(u'\u2030',#horizontalalignment='left',
                           fontsize=8,labelpad=-30)

       
        else:
            cmap_j = 'RdBu_r'
            vmin, vmax, vdiff = get_anom_plotvals()
            titlename = SLICES[i] + '-' + SLICES[0]
            cs = map.contourf(x, y, data_shift * 10., 
                          levels =np.arange(vmin,vmax,vdiff),
                          cmap = cmap_j, extend='both')
            plt.title(titlename,fontsize=10)
            cbar = plt.colorbar(cs, orientation='horizontal')
            cbar.ax.tick_params(labelsize='small')
            cbar.set_label('x 10' + u'\u2030', 
                           #horizontalalignment='left',
                           fontsize=8,labelpad=-30)



      
        
    plt.savefig(FILEOUTSTART + '.png')
    plt.savefig(FILEOUTSTART + '.eps')
    plt.savefig(FILEOUTSTART + '.pdf')
    plt.close()


def plot_winds(u_cubes, v_cubes, lev):
    """
    plots all the cubes : one from each timeslice
    """

      
    #fig = plt.figure(figsize=[12.8, 9.6])
    fig = plt.figure(figsize=[50, 40],dpi=200)
    for i, ucube in enumerate(u_cubes):
        vcube = v_cubes[i]
        
        subplotno = np.int('23' + np.str(i+1))
        ax = fig.add_subplot(subplotno)
        print(subplotno)
        
        map = Basemap(llcrnrlon=-20.0, urcrnrlon=20.0,
                      llcrnrlat=30.0, urcrnrlat=70.0,
                      projection='cyl', resolution='l')
        
        #map = Basemap(projection='cyl', resolution='l')
       

        latitudes = ucube.coord('latitude').points
        longitudes = ucube.coord('longitude').points
        
        u_shift, lons_shift = shiftgrid(
                180, ucube.data, longitudes, start=False)

        v_shift, lons_shift = shiftgrid(
                180, vcube.data, longitudes, start=False)
        
        lons, lats = np.meshgrid(lons_shift, latitudes)
        x, y = map(lons, lats)
        map.drawcoastlines()
       

        if i == 0: 
           titlename = SLICES[i]
           scalesize = 200
        else:
           scalesize = 50
           titlename = SLICES[i]  + '-' + SLICES[0]
        if lev == '850':
            scalesize = scalesize/3.0
            
        n=1 # plot every nth arrow
        Q = map.quiver(x[::n, ::n], y[::n, ::n], 
                       u_shift[::n, ::n], v_shift[::n, ::n], scale=scalesize,
                       headwidth=7, headlength=9)
        plt.title(titlename, fontsize=30, loc='left')
        qk = ax.quiverkey(Q, 0.7, 1.05, 10, ' 1 m/s', labelpos='E',
                          fontproperties={'size':30})
       
       
    
    plt.savefig(FILEOUTSTART + '.png')
    plt.savefig(FILEOUTSTART + '.pdf')
    plt.close()


def main():
    """
    Control: 1. See whether we are plotting winds
             2. read in data 
             3. setup data for plot
             4. plot : currenaly can only do a 6 panel figure
    """

    if FIELDREQ[0:2] == 'UV': # plotting winds
        winds = 'y'
        level = FIELDREQ[2:5]
        field1 = FIELDREQ[0:1] + ' COMPNT OF WIND AT ' + level + 'hPa'
        field2 = FIELDREQ[1:2] + ' COMPNT OF WIND AT ' + level + 'hPa'
        cubelist_u = get_allcubes(field1)            
        cubelist_v = get_allcubes(field2)            
        plot_winds(cubelist_u, cubelist_v, level)
    else:
        winds = 'n'
        field1 = FIELDREQ
        cubelist = get_allcubes(field1)            
        plot_data(cubelist)


############################################################
print('julia start')
SLICES = ['11ka', '10ka', '9ka', '8ka', '7ka', '6ka']
FIELDREQ = 'd18osw'
                    # 'UV200' 'UV850'
                    # 'TOTAL PRECIPITATION RATE    MM/DAY'
                    #'TEMPERATURE AT 1.5M'   
                    # evap
                    # d18o
LINUX_WIN = 'l'
NSLICES = len(SLICES)

shortfield = {'TOTAL PRECIPITATION RATE    MM/DAY' : 'precip',
              'TEMPERATURE AT 1.5M' : 'temp',
              'd18o' : 'd18o_p',
              'UV850' : 'winds_at_850hPa',
              'UV200' : 'winds_at_200hPa',
              'evap' : 'evap',
              'd18osw': 'd18o_sw'}
              
if FIELDREQ[0:2] == 'UV':
    filetype = 'winds_'
else:
    filetype = 'data_'

#if LINUX_WIN == 'w':
#   FILESTART = ('C:\\Users\\julia\\OneDrive\\WORK\\COLLABORATORS\\'
#                + 'JONATHAN_HOLMES\\data_slices\\' + filetype)
#   FILEOUTSTART = ('C:\\Users\\julia\\OneDrive\\WORK\\COLLABORATORS\\'
#                + 'JONATHAN_HOLMES\\plots\\' + FIELDREQ)
#else:
FILESTART = 'modeloutput/' + filetype                
FILEOUTSTART = 'modeloutput/plots/' + shortfield.get(FIELDREQ)

main()

