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
  
   

#####################################
def plot_warm_vs_precip(modelnames,modelsource,exptname,cntlname):
   
    
    namefield = {"NearSurfaceTemperature" : "SAT",
                 "TotalPrecipitation" : "Precipitation",
                 "SST" : "SST",
                 "SeaIceConc" : 'Sea ice fraction'
                 }
    p2expt = {'LP':'EOI400','PI':'E280'}
    
    # set up lists to store all values for both experiment and control    
    temp_global_mean=[]
    temp_global_sd=[]
    precip_global_mean=[]
    precip_global_sd=[]
    
    # store months in a numpy array (model, experiment,month)
    monmeans=np.zeros((len(modelnames),2,12))   
    monsd=np.zeros((len(modelnames),2,12))
    
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
        filenames=[]
        if modeluse == 'HadCM3-min':
             filenames.append(FILESTART+'HadCM3/LP_min_LSM_change')
             filenames.append(FILESTART+'HadCM3/'+cntlname)
        elif modeluse == 'MIROC4m-min':
             filenames.append(FILESTART+'MIROC4m/LP_min_LSM_change')
             filenames.append(FILESTART+'MIROC4m/'+ cntluse)
        else:
             filenames.append(FILESTART+modeluse+'/'+exptuse)
             filenames.append(FILESTART+modeluse+'/'+cntluse)
       
        # get temperature
        means=[]
        sds=[]
        nmon=12
                        
        for fileno in range(0,len(filenames)):

            f=open(filenames[fileno]+".NearSurfaceTemperature.data.txt","r")
            f1=f.readlines()
            f2 = [x.replace('\n', '') for x in f1]
            
            # get the means according to their position in the file
            all_mean_sd=f2[2]
            all_mon_mean_sd=f2[5:5+12]
            all_lat_mean_sd=f2[20:20+180]
           
            # extract global mean
            mean,sd=all_mean_sd.split(',')
            means.append(mean)
            sds.append(sd)
            
        temp_global_mean.append(means)
        temp_global_sd.append(sds)

        # get preciptation
        means=[]
        sds=[]
        nmon=12
                        
        for fileno in range(0,len(filenames)):

            f=open(filenames[fileno]+".TotalPrecipitation.data.txt","r")
            f1=f.readlines()
            f2 = [x.replace('\n', '') for x in f1]
            
            # get the means according to their position in the file
            all_mean_sd=f2[2]
            all_mon_mean_sd=f2[5:5+12]
            all_lat_mean_sd=f2[20:20+180]
           
            # extract global mean
            mean,sd=all_mean_sd.split(',')
            means.append(mean)
            sds.append(sd)
            
        precip_global_mean.append(means)
        precip_global_sd.append(sds)
       
    
  

    #############################################################
    # plot warming vs precip
   
    expt_temp_mean=[float(item[0]) for item in temp_global_mean]
    expt_temp_2sigma=([float(item[0])*2.0 for item in temp_global_sd])

    cntl_temp_mean=[float(item[1]) for item in temp_global_mean]
    cntl_temp_2sigma=([float(item[1])*2.0 for item in temp_global_sd])

    anom_temp_mean = [e - c for e, c in zip(expt_temp_mean, cntl_temp_mean)]

    expt_precip_mean=[float(item[0]) for item in precip_global_mean]
    expt_precip_2sigma=([float(item[0])*2.0 for item in precip_global_sd])

    cntl_precip_mean=[float(item[1]) for item in precip_global_mean]
    cntl_precip_2sigma=([float(item[1])*2.0 for item in precip_global_sd])

    anom_precip_mean = [e - c for e, c in zip(expt_precip_mean,
                                              cntl_precip_mean)]
 
    
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
              anom_temp_mean[i], # x-value
              anom_precip_mean[i], # y-value
              marker=markers[i % len(markers)],
              s=80,
              facecolors=facecolor,
              edgecolors=edgecolor,
              linewidth=1.5,
              label=model
         )

  
    ax.set_xlabel('LP - PI temperature (°C)')
    ax.set_ylabel('LP - PI precipitatino (mm day)')
    #ax.axhline(0, color='k', linestyle='--', linewidth=0.5)
    
    box = ax.get_position()
   # ax.set_position([box.x0, box.y0+(0.15*box.height), box.width * 0.8, box.height])
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0,0.5))
   
    
    titlename='temperature anomaly vs precipitation anomaly'
    ax.set_title(titlename)
    #plt.show()
    #sys.exit(0)
    plt.savefig('plots/warming_vs_precip.png')
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

   
   plot_warm_vs_precip(modelnames,modelsource,exptname,cntlname)

#sys.exit(0)

