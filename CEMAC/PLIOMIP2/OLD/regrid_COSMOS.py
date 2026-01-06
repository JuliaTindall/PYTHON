#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 14 14:13:50 2019

#@author: earjcti
#
# This program will regrid the COSMOS data that is needed for PLIOMIP2.  
# We will put 100 year average fields onto a 1deg X 1deg standard grid.
# Basically this data was uploaded with 100 years of monthly data in a single
# file
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
import iris.analysis.cartography
import sys

#####################################
def regrid_data(fieldnamein,fieldnameout,exptnamein,exptnameout,filename):

   
    outname=('/nfs/hera1/earjcti/regridded/COSMOS/'+exptnameout+'_'+fieldnameout+'.nc')

    cube=iris.load_cube(filename)
    
    nt,ny,nx=np.shape(cube)
    
    # get the annual mean data in a numpy array
    nmonths=12
    nyears=nt/nmonths
    #print(nt,nmonths,nyears)
    
    reform_data=np.zeros((nmonths,nyears,ny,nx))
    for year in range(0,nyears):
        for month in range(0,nmonths):
            t=(year*12)+month
            reform_data[month,year,:,:]=cube.data[t,:,:]
            
    mean_data=np.mean(reform_data,axis=1)
    #print(np.shape(mean_data))
    
    #######################################################################
    # we are going to try and change the original cube to one we want 
    # in order to preserve metadata
    
    # change the time coordinate to be an ascending sequence
    # set up a subcube with just 12 time points
    times=np.arange(0,nt)
    cube.coord('time').points=times
    subcube=cube[0:12,:,:]
    
    # replace subcube data with meandata
    subcube.data=mean_data
   
    # now regrid the cube onto a 1X1 grid
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'
    
    cubegrid=iris.load_cube('one_lev_one_deg.nc')
    regrid_mean=subcube.regrid(cubegrid,iris.analysis.Linear())
    
   
    print(outname)
    iris.save(regrid_mean,outname,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
     
    
    ########################################################
    # check that we have averaged properly.  To do this we are 
    # going to plot the annual cycle of the global mean field
    # for the regridded dataset and also for each year
   
    #global mean
    plt.subplot(1,2,1) # global mean from each year
    for i in range(0,100):
        subcube.data=reform_data[:,i,:,:]
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
    plt.title('globavg '+fieldnamein)
    plt.legend()
    
    # check at 30N
    plt.subplot(1,2,2)
    for i in range(0,100): # mean at 30N from each year
        subcube.data=reform_data[:,i,:,:]
        slice_30N_subcube= subcube.extract(iris.Constraint(latitude=32))
        mean_30N=slice_30N_subcube.collapsed(['longitude'], iris.analysis.MEAN)
        plt.plot(mean_30N.data,color='r')
    
    #mean at 30N
    slice_30N= regrid_mean.extract(iris.Constraint(latitude=32))
    mean_30N=slice_30N.collapsed(['longitude'], iris.analysis.MEAN)
    
   
    plt.plot(mean_30N.data,color='b',label='avg')
    plt.title('average at 30N by month')
    plt.legend()
    plt.show()


##########################################################
# main program

# dictionaries for COSMOS data
exptname = {
        "E280" : "E280",
        "Eoi400" : "EOI400"}

fieldname = {
        "pr" : "TotalPrecipitation",
        "ts" : "SurfaceTemperature",
        "sic" : "SeaIceConcentration"
        }

COSMOS_FIELDS ={"pr" : "TotalPrecip",
        "ts" : "SurfaceTemp",
        "sic" : "SeaIceAreaFraction"
        }

# this is regridding where all results are in a single file
fieldnamein=['pr','ts','sic']
exptnamein=['Eoi400','E280']

for field in range(0,len(fieldnamein)):
    for expt in range(0,len(exptnamein)):
        filename=('/nfs/hera1/pliomip2/data/AWI/COSMOS/'+exptnamein[expt]+
                  '/'+exptnamein[expt]+
                  '.'+COSMOS_FIELDS.get(fieldnamein[field])+
                  '_CMIP6_name_'+fieldnamein[field]+
                  '_2650-2749_monthly_mean_time_series.nc')
#filename='C:\\Users\\julia\\OneDrive\\DATA\\COSMOS_DATA\\E280.TotalPrecip_CMIP6_name_pr_2650-2749_monthly_mean_time_series.nc'

        fieldnameout=fieldname.get(fieldnamein[field])
        exptnameout=exptname.get(exptnamein[expt])
        regrid_data(fieldnamein[field],fieldnameout,exptnamein[expt],exptnameout,filename,)

#sys.exit(0)