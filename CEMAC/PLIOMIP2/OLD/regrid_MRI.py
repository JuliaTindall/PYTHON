#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 14 14:13:50 2019

#@author: earjcti
#
# This program will regrid the MRI data that is needed for PLIOMIP2.  
# We will put 100 year average fields onto a 1deg X 1deg standard grid
#

import os
import numpy as np
import scipy as sp
#import cf
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys

#####################################
def regrid_data(fieldnamein,fieldnameout,exptnamein,exptnameout,plotno):

    # get filenames
    filestart='/nfs/hera1/pliomip2/data/MRI-CGCM2.3/'
    filename=(filestart+'/'+fieldnamein+'/'+exptnamein+'.'+fieldnamein+'.*.nc')
    outname=('/nfs/hera1/earjcti/regridded/MRI-CGCM2.3/'+exptnameout+'_'+fieldnameout+'.nc')

    #filename='C:\\Users\\julia\\OneDrive\\DATA\\MRI_DATA\\*.nc'

    # load data into iris cubes        
    cubes=iris.load(filename)
    ncubes=len(cubes)
    
    
    # we are going to average the cubes manually as merge doesn't seem to work
    count=0
    for i in range(0,ncubes):  
        if i ==0:
            mean_cubes=cubes[i]
        else:
            mean_cubes=mean_cubes+cubes[i]
        count=count+1
        #print(i,mean_cubes.data[:,0,0,0])
    
    mean_cubes=mean_cubes/count
    mean_cubes.long_name=cubes[i].long_name
    
    # now regrid the cube onto a 1X1 grid
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'
    
    cubegrid=iris.load_cube('one_lev_one_deg.nc')
    regrid_mean=mean_cubes.regrid(cubegrid,iris.analysis.Linear())
    
    
    print(outname)
    iris.save(regrid_mean,outname,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
       
    
     ########################################################
    # check that we have averaged properly.  To do this we are 
    # going to plot the annual cycle of the global mean field
    # for the regridded dataset and also for each year
   
    #global mean
    plt.subplot(1,2,plotno) # global mean from each year
    subcube=cubes[0]
    for i in range(0,ncubes):
        subcube.data=cubes[i].data
        if i ==0:
            subcube.coord('latitude').guess_bounds() 
            subcube.coord('longitude').guess_bounds()
            grid_areas = iris.analysis.cartography.area_weights(subcube)
        subcube_mean = subcube.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
        plt.plot(subcube_mean.data,color='r')
        
    # global mean from average
   
    regrid_mean.coord('latitude').guess_bounds() 
    regrid_mean.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(regrid_mean)
    temporal_mean = regrid_mean.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
   
   
    plt.plot(temporal_mean.data,color='b',label='avg')
    plt.title('global average by month')
    plt.legend()
    plt.show()
    
   



##########################################################
# main program

# dictionaries for MRI.py
exptname = {
        "e280" : "E280",
        "eoi400" : "EOI400"}

fieldname = {
        "pr" : "TotalPrecipitation",
        "evap" : "TotalEvap",
        "psl" : "MSLP",
        "hfss" : "SensibleHeatFlux",
        "clt" : "TotalCloudCover",
        "tas" :"SurfaceTemperature",
        "ts" : "SurfaceSkinTemperature",
        "rsdt" : "DownSWTOA",
        "rsut" : "UpSWTOA",
        "rlut" : "OLR_TOA",
        "rsus" : "SurfSWUp", 
        "rlus" : "SurfLWUp",
        "rsds" : "SurfSWDown",
        "rlds" : "SurfLWDown",
        "rsuscs" : "ClearSkySurfSWUp",
        "rsdscs" : "ClearSkySurfSWDown", 
        "rldscs" :  "ClearSkySurfLWDown",
        "sic" : "SeaIceConcentration",
        "tos" : "SeaSurfaceTemp"
            }


exptnamein=['eoi400','e280']
fieldnamein=(['pr','evap','psl','hfss','clt','tas','ts',
              'rsdt','rsut','rlut','rsus','rlus','rsds','rlds',
              'rsuscs','rsdscs','rldscs','sic','tos'] )

for i in range(0,len(exptnamein)):
    for field in range(0,len(fieldnamein)):

        fieldnameout=fieldname.get(fieldnamein[field])
        exptnameout=exptname.get(exptnamein[i])

        regrid_data(fieldnamein[field],fieldnameout,exptnamein[i],exptnameout,i+1)

#sys.exit(0)