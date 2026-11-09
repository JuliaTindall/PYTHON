#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 29th


#@author: earjcti
#
# This program will plot warming vs increased precipitation

import os
import numpy as np
import scipy as sp
#import cf
import iris
import iris.util
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
#from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
import cf_units as unit
from iris.cube import CubeList
#from iris.experimental.equalise_cubes import equalise_attributes
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import sys


#filestart = '/uolstore/Research/a/hera1/earjcti/'
filestart ='/home/earjcti/'

def plotmean_newaxis(cube, modelno_):
     tempcube=iris.util.new_axis(cube)
     tempcube.add_dim_coord(iris.coords.DimCoord(modelno_, 
            standard_name='model_level_number', long_name='model', 
            var_name='model', 
            units=None,
            bounds=None,
            coord_system=None, circular=False),0) 
     return tempcube
 

def resort_coords(cube,levelno):
    """
    this will make all the dimensions of the cube match.  They will all be
    longitude, latitude, level-no (ie 1 for first model, 2 for second model...)
    
    input is the cube and the level number
    output is the cube with the new dimensions
    here we will also equalise the units
    """
    if fieldnames[field] == 'SeaIceConc':
         cube.units = None
    
    for coord in cube.coords():  
        print(coord)
        name=coord.standard_name
        if name !='latitude' and name!='longitude':
             
            cube.remove_coord(coord)
                
    for coord in cube.coords():   # now this will be longitude or latitude
        coord.points=coord.points.astype('float32') 
        coord.var_name = coord.standard_name
        coord.long_name = coord.standard_name
         
    newcube = plotmean_newaxis(cube, [levelno])
    # this will make sure cell_methods match and that cubes can
    # be concatenated
    newcube.cell_methods = None
    
        
    return newcube
 
def get_pliomip1_data(fieldreq):
    """
    we will get the data from pliomip1
    returns an array of the mean data and the min, max mean, of the seasonal cycle
    """
    
    PLIOMIP1_FILE = (FILESTART[:-10] + 'PLIOMIP1/means_for_' 
                     + fieldreq + '.txt')
    
    print('julia',PLIOMIP1_FILE)
    sys.exit(0)
    f1 = open(PLIOMIP1_FILE)
    
    lines = f1.readlines()
    lines[:] = [line.rstrip('\n') for line in lines]
    
    # means are at the start of the file
    allanoms_list = []
    for i in range(1, len(lines)):
        line = lines[i]
        if line[0:9] == 'modelname':
            break    # we have now got all the means
        modname, eoi400, e280, anom = line.split(',')
        if modname == 'MEAN':
            meananom=anom
        else:
            allanoms_list.append(anom)
            
    # find line which contains 'jan feb mar ' which is the start of the seasonal cycle
    string = 'jan feb mar'
    min_seas_cyc = np.zeros(12) + 1000.
    max_seas_cyc = np.zeros(12) - 1000.
    for i, line in enumerate(lines):
        if string in line:
            index = i
            
    for i in range(index + 1, len(lines)):
        line = lines[i]
        if 'modelname' in line:
            break   
        modname, anomstr = line.split(',')
        anom_arr = np.array(anomstr.strip('[]').split(), dtype=float)
        if modname == 'MEAN':
            mean_seas_cyc = anom_arr
        else:
            for j in range(0,12):
                min_seas_cyc[j] = np.min([anom_arr[j], min_seas_cyc[j]]) 
                max_seas_cyc[j] = np.max([anom_arr[j], max_seas_cyc[j]]) 
    

    allanoms = np.asarray(allanoms_list, dtype=float)
    
        
    return allanoms, min_seas_cyc, max_seas_cyc, mean_seas_cyc
  

def calc_polar_amp(cube_expt, cube_cntl):
    """
    calculates the polar amplification between experiment cube and control cube
    """
    cube_anom = cube_expt - cube_cntl

    cube_anom.coord('latitude').guess_bounds()
    cube_anom.coord('longitude').guess_bounds()
  
    grid_areas = iris.analysis.cartography.area_weights(cube_anom)
    nlat=len(cube_anom.coord('latitude').points)
    lats=cube_anom.coord('latitude').points
        
    grid_areas_nh_polar=np.zeros(grid_areas.shape)
    grid_areas_sh_polar=np.zeros(grid_areas.shape)
    polarval=60.0
    for j in range(0,nlat):
         if lats[j] <= -1.0*polarval:
              grid_areas_sh_polar[j,:]=grid_areas[j,:]
         if lats[j] >= polarval:
              grid_areas_nh_polar[j,:]=grid_areas[j,:]
           
    nh_polar = cube_anom.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_nh_polar)
    sh_polar = cube_anom.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sh_polar)
    mean_cube = cube_anom.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas)

    #print(nh_polar.data)
    #print(sh_polar.data)
    #print(mean_cube.data)

    NH_amp = nh_polar.data / mean_cube.data
    SH_amp = sh_polar.data / mean_cube.data
      

    return NH_amp, SH_amp, mean_cube.data

#####################################
def plot_polaramp(modelnames,modelsource,exptname,cntlname):
   
    
    namefield = {"NearSurfaceTemperature" : "SAT",
                 "TotalPrecipitation" : "Precipitation",
                 "SST" : "SST",
                 "SeaIceConc" : 'Sea ice fraction'
                 }
    p2expt = {'LP':'EOI400','PI':'E280'}
    
    # set up lists to store all values for both experiment and control    
    exptcubes = CubeList([])
    cntlcubes = CubeList([])
    NH_all_amp = []
    SH_all_amp = []
    mean_warm_all = []
    
    for modelno in range(0,len(modelnames)):
        if modelsource[modelno] == 'P3':
             FILESTART = FILESTARTP3
             exptuse = exptname
             cntluse = cntlname
        if modelsource[modelno] == 'P2':
             FILESTART = FILESTARTP2
             exptuse = p2expt.get(exptname)
             cntluse = p2expt.get(cntlname)
        modeluse=modelnames[modelno]
        print(modeluse)
        exptfiles=[]
        cntlfiles=[]
        if modeluse == 'HadCM3-min':
             exptfiles.append(FILESTART+'HadCM3/LP_min_LSM_change')
             cntlfiles.append(FILESTART+'HadCM3/'+cntlname)
        elif modeluse == 'MIROC4m-min':
             exptfiles.append(FILESTART+'MIROC4m/LP_min_LSM_change')
             cntlfiles.append(FILESTART+'MIROC4m/'+ cntluse)
        else:
             exptfiles.append(FILESTART+modeluse+'/'+exptuse)
             cntlfiles.append(FILESTART+modeluse+'/'+cntluse)
       
                        
        for fileno,exptfile in enumerate(exptfiles):

            file_e = exptfile + ".NearSurfaceTemperature.allmean.nc"
            file_c = cntlfiles[fileno]+".NearSurfaceTemperature.allmean.nc"

            cube_e = iris.load_cube(file_e)
            cube_c = iris.load_cube(file_c)
            
            NH_amp, SH_amp, meanwarm = calc_polar_amp(cube_e,cube_c)
            NH_all_amp.append(NH_amp)
            SH_all_amp.append(SH_amp)
            mean_warm_all.append(meanwarm)
            
 
    
    #######################################################
    # we are going to plot anomaly vs Pliocene temperature
    # and use a symbol for the model
    #fig,ax=plt.subplots(2,1,1)
    
    markers = [
    'o', 's', '^', 'P', 'v', '>',     # circle, square, triangles
    'D', 'd',                         # diamonds
    'H', '<',                         # pentagons / filled plus
    'p', 'h', '*',                    # star, hexagons
    'X']                    # X, x, plus

    colors = plt.cm.tab20.colors # 20 distinct colours
    
    fig, ax = plt.subplots(figsize=[8,8])

    for i, model in enumerate(modelnames):
         color = colors[i % len(colors)]
         if modelsource[i] == 'P3':
              facecolor = color
              edgecolor = color
         else: # P2
              facecolor = 'white'
              edgecolor = color
         ax.scatter(
              SH_all_amp[i], # x-value
              NH_all_amp[i], # y-value
              marker=markers[i % len(markers)],
              s=80,
              facecolors=facecolor,
              edgecolors=edgecolor,
              linewidth=1.5,
              label=model
         )

  
    ax.set_ylabel('NH polar amplification')
    ax.set_xlabel('SH polar amplification')
    
    box = ax.get_position()
   # ax.set_position([box.x0, box.y0+(0.15*box.height), box.width * 0.8, box.height])
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0,0.5))


    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()

    lo = min(xmin, ymin)
    hi = max(xmax, ymax)

    ax.plot([lo, hi], [lo, hi], 'k--', linewidth=1)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    
    titlename='Polar amplification (NH vs SH)'
    ax.set_title(titlename)
    plt.show()
    sys.exit(0)
    plt.savefig('plots/polar_amplification/NH_vs_SH.png')
    plt.close()

    #####################################
    # plot NH polar amp vs mean warming
    
    fig, ax = plt.subplots(figsize=[8,8])

    for i, model in enumerate(modelnames):
         color = colors[i % len(colors)]
         if modelsource[i] == 'P3':
              facecolor = color
              edgecolor = color
         else: # P2
              facecolor = 'white'
              edgecolor = color
         ax.scatter( mean_warm_all[i],
              NH_all_amp[i],
              marker=markers[i % len(markers)],
              s=80,
              facecolors=facecolor,
              edgecolors=edgecolor,
              linewidth=1.5,
              label=model
         )

  
    ax.set_xlabel('mean warming')
    ax.set_ylabel('NH polar amplification')
    
    box = ax.get_position()
   # ax.set_position([box.x0, box.y0+(0.15*box.height), box.width * 0.8, box.height])
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
   
    
    titlename='Polar amplification (NH vs mean warming)'
    ax.set_title(titlename)
    plt.savefig('plots/polar_amplification/NH_vs_mean.png')
    plt.close()

    #####################################
    # plot SH polar amp vs mean warming
    
    fig, ax = plt.subplots(figsize=[8,8])

    for i, model in enumerate(modelnames):
         color = colors[i % len(colors)]
         if modelsource[i] == 'P3':
              facecolor = color
              edgecolor = color
         else: # P2
              facecolor = 'white'
              edgecolor = color
         ax.scatter( mean_warm_all[i],
              SH_all_amp[i],
              marker=markers[i % len(markers)],
              s=80,
              facecolors=facecolor,
              edgecolors=edgecolor,
              linewidth=1.5,
              label=model
         )

  
    ax.set_xlabel('mean warming')
    ax.set_ylabel('SH polar amplification')
    
    box = ax.get_position()
   # ax.set_position([box.x0, box.y0+(0.15*box.height), box.width * 0.8, box.height])
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
   
    
    titlename='Polar amplification (SH vs mean warming)'
    ax.set_title(titlename)
    plt.savefig('plots/polar_amplification/SH_vs_mean.png')
    plt.close()

    
    sys.exit(0)
    plt.close()

##########################################################
# main program
        
filename=' '
linux_win='l'


P3modelnames=['ACCESS', 'AWI-CM2', 
            'AWI-ESM2', 'COSMOS', 'HadCM3', 'HadCM3-min','MIROC4m','MIROC4m-min']
P2modelnames = [  'CCSM4',  'CCSM4-UoT','CCSM4-Utr',  
                  'CESM2','COSMOS', 'CESM1.2',
                  'EC-Earth3.3', 'GISS2.1G','HadCM3',
                  'HadGEM3','IPSLCM5A','IPSLCM5A2',
                  'IPSLCM6A', 'MIROC4m','MRI2.3', 
                  'NorESM-L',  'NorESM1-F']
            
modelnames = P3modelnames + P2modelnames
modelsource = ['P3'] * len(P3modelnames) + ['P2'] * len(P2modelnames)
PLIOMIP1 = 'n'

exptname='LP'
cntlname='PI'

if linux_win=='w':
   FILESTART='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded/'
else:
   FILESTARTP3=filestart + '/regridded100_P3/'
   FILESTARTP2=filestart + '/regridded100/'

   
   plot_polaramp(modelnames,modelsource,exptname,cntlname)

#sys.exit(0)

