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

   
    #outname=('/nfs/hera1/earjcti/regridded/MIROC/'+exptnameout+'_'+fieldnameout+'.nc')
    outstart=exptnameout+'_'+fieldnameout
    
    cube=iris.load_cube(filename)
    cube_12=cube[0:12,:,:]
    
   
    #print(cube)
    #print(cube.coord('latitude').bounds)
    #print(cube.coord('lat_bnds'))
    #sys.exit(0)
    
    
  
    # now regrid the cube onto a 1X1 grid (we will first try regridding the raw data)
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'
    
    cubegrid=iris.load_cube('one_lev_one_deg.nc')
    regridded_cube=cube.regrid(cubegrid,iris.analysis.Linear())
    #print(regridded_cube)
    #print(np.shape(cube))
    #print(np.shape(regridded_cube))
   

    # reshape the timeseries into months and years
    nt,nyorig,nxorig=np.shape(cube)
    nt,ny,nx=np.shape(regridded_cube)
    # get the annual mean data in a numpy array
    nmonths=12
    nyears=nt/nmonths
     
    reform_data=np.zeros((nmonths,nyears,ny,nx))   
    reform_orig=np.zeros((nmonths,nyears,nyorig,nxorig))
    for year in range(0,nyears):
        for month in range(0,nmonths):
            t=(year*12)+month
            reform_data[month,year,:,:]=regridded_cube.data[t,:,:]  
            reform_orig[month,year,:,:]=cube.data[t,:,:]
            
    
    # then get mean and standard deviation
    mean_month_regridded_data=np.mean(reform_data,axis=1)
    mean_regridded_data=np.mean(mean_month_regridded_data,axis=0)
    sd_month_regridded_data=np.std(reform_data,axis=1)
    sd_regridded_data=np.std(np.mean(reform_data,axis=0),axis=0)
    #print(np.shape(mean_data))
    
    #######################################################################
    # we are going to try and change the original cube to one we want 
    # in order to preserve metadata
    
    # change the time coordinate to be an ascending sequence
    # set up a subcube with just 12 time points
    times=np.arange(0,nt)
    regridded_cube.coord('time').points=times
    subcube_mean_mon=regridded_cube[0:12,:,:]
    subcube_sd_mon=regridded_cube[0:12,:,:]
    subcube_mean=regridded_cube[0:1,:,:]
    subcube_sd=regridded_cube[0:1,:,:]
    
    # replace subcube data with meandata
   
    subcube_mean_mon.data=mean_month_regridded_data
    subcube_sd_mon.data=sd_month_regridded_data
    subcube_mean.data[0]=mean_regridded_data
    subcube_sd.data[0]=sd_regridded_data
   
   
    print(outstart)
    outmeanmon=outstart+'mean_mon.nc'
    outsdmon=outstart+'sd_mon.nc'
    outmean=outstart+'mean.nc'
    outsd=outstart+'sd.nc'
    iris.save(subcube_mean_mon,outmeanmon,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
    iris.save(subcube_sd_mon,outsdmon,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
    iris.save(subcube_mean,outmean,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
    iris.save(subcube_sd,outsd,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
    
    print('mean1',subcube_mean_mon.data[0,0,0])
    print('sd1',subcube_sd_mon.data[0,0,0])
   
    ########################################################
    # check that we have averaged properly.  To do this we are 
    # going to plot the annual cycle of the global mean field
    # for the regridded dataset and also for each year
    
   
   
    #global mean
    plt.subplot(2,2,1) # global mean from each year
    #subcube=subcube_mean_mon.copy(data=)  # set up structure of subcube
    
    print('mean2',subcube_mean_mon.data[0,0,0])
    print('sd2',subcube_sd_mon.data[0,0,0])
    for i in range(0,100):
        subcube=cube_12.copy(data=reform_orig[:,i,:,:])
        subcube.coord('latitude').bounds
        subcube.coord('longitude').bounds
        grid_areas = iris.analysis.cartography.area_weights(subcube)
        newdata_plot = subcube.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
        plt.plot(newdata_plot.data,color='r')
        
    print('mean3',subcube_mean_mon.data[0,0,0])
    print('sd3',subcube_sd_mon.data[0,0,0]) 
    # global mean from average
    subcube_mean_mon.coord('latitude').guess_bounds()
    subcube_mean_mon.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(subcube_mean_mon)
    temporal_mean = subcube_mean_mon.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
   
    plt.plot(temporal_mean.data,color='b',label='avg')
    plt.title('globavg '+fieldnamein)
    plt.legend()
    print('mean5',subcube_mean_mon.data[0,0,0])
    print('sd5',subcube_sd_mon.data[0,0,0])
    
    # check at 30N
    plt.subplot(2,2,2)
    for i in range(0,100): # mean at 30N from each year
        subcube=cube_12.copy(data=reform_orig[:,i,:,:])
        slice_30N_subcube= subcube.extract(iris.Constraint(latitude=30))
        mean_30N=slice_30N_subcube.collapsed(['longitude'], iris.analysis.MEAN)
        plt.plot(mean_30N.data,color='r')
    
    #mean at 30N
    slice_30N= subcube_mean_mon.extract(iris.Constraint(latitude=30))
    mean_30N=slice_30N.collapsed(['longitude'], iris.analysis.MEAN)
    
   
    plt.plot(mean_30N.data,color='b',label='avg')
    plt.title('average at 30N by month')
    plt.legend()
    

#   standard deviation
    plt.subplot(2,2,3) # global st from each year
    subcube=cube_12.copy(data=np.std(reform_orig[:,:,:,:],axis=1))
    subcube.coord('latitude').bounds
    subcube.coord('longitude').bounds
    grid_areas = iris.analysis.cartography.area_weights(subcube)
    newdata_plot = subcube.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
    plt.plot(newdata_plot.data,color='r',label='orig')
        
    # global sd from average
    subcube_sd_mon.coord('latitude').guess_bounds()
    subcube_sd_mon.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(subcube_sd_mon)
    temporal_sd = subcube_sd_mon.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
   
    print('mean9',subcube_mean_mon.data[0,0,0])
    print('sd9',subcube_sd_mon.data[0,0,0])
    plt.plot(temporal_sd.data,color='b',linestyle='--',label='regridded')
    plt.title('globavg '+fieldnamein)
    plt.legend()
    
    # check at 30N
    plt.subplot(2,2,4)
    subcube=cube_12.copy(data=np.std(reform_orig[:,:,:,:],axis=1))
    slice_30N_subcube= subcube.extract(iris.Constraint(latitude=30))
    sd_30N=slice_30N_subcube.collapsed(['longitude'], iris.analysis.MEAN)
    plt.plot(sd_30N.data,color='r',label='orig')
    
    
    #sd at 30N
    slice_30N= subcube_sd_mon.extract(iris.Constraint(latitude=30))
    sd_30N=slice_30N.collapsed(['longitude'], iris.analysis.MEAN)
    
   
    plt.plot(sd_30N.data,color='b',label='regridded')
    plt.title('average at 30N by month')
    plt.legend()
    plt.show()



##########################################################
# main program

# dictionaries for MIROC data
exptname = {
        "E280" : "E280",
        "Eoi400" : "EOI400"}

fieldname = {
        "pr" : "TotalPrecipitation",
        "ts" : "SurfaceTemperature",
        "sic" : "SeaIceConcentration"
        }

MIROC_FIELDS ={"pr" : "pr",
        "ts" : "tas",
        "sic" : "SeaIceAreaFraction"
        }

# this is regridding where all results are in a single file
#fieldnamein=['pr','ts','sic']
#exptnamein=['Eoi400','E280']

fieldnamein=['pr','ts']
exptnamein=['Eoi400']

for field in range(0,len(fieldnamein)):
    for expt in range(0,len(exptnamein)):
        #filename=('/nfs/hera1/pliomip2/data/AWI/COSMOS/'+exptnamein[expt]+
        #          '/'+exptnamein[expt]+
        #          '.'+COSMOS_FIELDS.get(fieldnamein[field])+
        #          '_CMIP6_name_'+fieldnamein[field]+
        #          '_2650-2749_monthly_mean_time_series.nc')
        filename=('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\MIROC_DATA\\'+
            'MIROC4m_'+exptnamein[expt]+'_Amon_'+
            MIROC_FIELDS.get(fieldnamein[field])+'.nc')

        fieldnameout=fieldname.get(fieldnamein[field])
        exptnameout=exptname.get(exptnamein[expt])
        regrid_data(fieldnamein[field],fieldnameout,exptnamein[expt],exptnameout,filename,)

#sys.exit(0)