#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 14 14:13:50 2019

#@author: earjcti
#
# This program will regrid some of the data that is needed for PLIOMIP2.  
# We will put 100 year average fields onto a 1deg X 1deg standard grid
# it can be used where experiments have been uploaded as one file per yer
#
# currently it can do: HadCM3 and MRI
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
from iris.experimental.equalise_cubes import equalise_attributes
import sys



#####################################
def regrid_data(fieldnamein,fieldnameout,exptnamein,exptnameout,plotno,modelname):

    # get filenames and setup model specfiic variables
    
  
    filestart='/nfs/hera1/pliomip2/data/'
    
    if modelname=='HadCM3':
        modeluse='LEEDS/HadCM3'
        zname='surface'
        xname='longitude'
        yname='latitude'
        tname='t'
        filestart=filestart+modeluse+'/'+exptnamein+'/'
    
    if modelname=='MRI-CGCM2.3':
        modeluse=modelname    
        zname='pressure level'
        xname='lon'
        yname='lat'
        tname='time'
        filestart=filestart+modelname+'/'
 
        
        
   
    filename=(filestart+'/'+fieldnamein+'/'+exptnamein+'.'+fieldnamein+'.*.nc')

    #filestart='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'+modelname+'_DATA\\'
    #filename=filestart+exptnamein+'.'+fieldnamein+'*.nc'
    outstart=('/nfs/hera1/earjcti/regridded/'+modelname+'/'+exptnameout+'.'+
        fieldnameout+'.')

    # load data into iris cubes        
    cubes=iris.load(filename)
    ncubes=len(cubes)
    print('there are ',ncubes,' cubes')
    print(cubes[0].coord)
    # change z dimension to year
    for i in range(0,len(cubes)):
        # put the year on the surface coordinate
        cubes[i].coord(zname).points=i
        
    # now regrid the cubes onto a 1X1 grid
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'
    
    cubegrid=iris.load_cube('one_lev_one_deg.nc')
    
    # setup a cubelist to store regridded temperature
    regrid_temp=iris.cube.CubeList([])
    for i in range(0,ncubes):
        regrid=cubes[i].regrid(cubegrid,iris.analysis.Linear())
        regrid_temp.append(regrid)
    
    
        
    #make sure the metadata on all cubes are the same
    equalise_attributes(regrid_temp)
    regrid_cubes=regrid_temp.concatenate_cube()
    regrid_cubes.coord(zname).rename('year')
    regrid_cubes.coord(tname).rename('month')
    regrid_cubes.coord('month').points=(np.arange(0,12))

       
      # extract monthly data from cubes
    print(regrid_cubes.coord('month'))
    sys.exit(0)
    jan_slice = regrid_cubes.extract(iris.Constraint(month='Jan'))
    feb_slice = regrid_cubes.extract(iris.Constraint(month='Feb'))
    mar_slice = regrid_cubes.extract(iris.Constraint(month='Mar'))
    apr_slice = regrid_cubes.extract(iris.Constraint(month='Apr'))        
    may_slice = regrid_cubes.extract(iris.Constraint(month='May'))
    jun_slice = regrid_cubes.extract(iris.Constraint(month='Jun'))
    jul_slice = regrid_cubes.extract(iris.Constraint(month='Jul'))
    aug_slice = regrid_cubes.extract(iris.Constraint(month='Aug'))
    sep_slice = regrid_cubes.extract(iris.Constraint(month='Sep'))
    oct_slice = regrid_cubes.extract(iris.Constraint(month='Oct'))
    nov_slice = regrid_cubes.extract(iris.Constraint(month='Nov'))
    dec_slice = regrid_cubes.extract(iris.Constraint(month='Dec'))
 
  
    # will give a mean annual cycle
    mean_mon_data=regrid_cubes.collapsed('year',iris.analysis.MEAN)
    sd_mon_data=regrid_cubes.collapsed('year',iris.analysis.STD_DEV)
    
    # will put the annual average for each year in the 'surface' column
    mean_year_data=regrid_cubes.collapsed('month',iris.analysis.MEAN)
    
    # find mean and standard deviation for the 100 years
    mean_data=mean_mon_data.collapsed('month',iris.analysis.MEAN)
    sd_data=mean_year_data.collapsed('year',iris.analysis.STD_DEV)
    
    
    # write the cubes out to a file
    
    outfile=outstart+'mean_month.nc'
    iris.save(mean_mon_data,outfile,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
       
    outfile=outstart+'sd_month.nc'
    iris.save(sd_mon_data,outfile,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
       
    outfile=outstart+'allmean.nc'
    iris.save(mean_data,outfile,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
       
    outfile=outstart+'allstdev.nc'
    iris.save(sd_data,outfile,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
    
     ##########################################################################
    # get the global mean and standard deviation and write them all out to a file
    #
    textout=outstart+'data.txt'
    
    file1= open(textout,"w") 
    
    # get mean field for cube
            
    mean_data.coord('latitude').guess_bounds()
    mean_data.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mean_data)
    tempcube=mean_data.collapsed(['latitude','longitude'], 
                                iris.analysis.MEAN,weights=grid_areas)
    meanann=tempcube.data
    
    # get standard deviation
    # 1. mean for each year
    
    mean_year_data.coord('latitude').guess_bounds()
    mean_year_data.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mean_year_data)
    tempcube=mean_year_data.collapsed(['latitude','longitude'], 
                                iris.analysis.MEAN,weights=grid_areas)
    stdevcube=tempcube.collapsed(['year'],iris.analysis.STD_DEV)
    stdevann=stdevcube.data
    
    print('meanstdev',meanann,stdevann)
    plt.plot(tempcube.data)
    plt.plot([0,100],[meanann,meanann])
    plt.plot([0,100],[meanann+stdevann+stdevann,meanann+stdevann+stdevann])
    plt.plot([0,100],[meanann-stdevann-stdevann,meanann-stdevann-stdevann])
    plt.title('data and data+/-2sd')

    # write out to a file
    file1.write('global annual mean and standard deviation\n')
    file1.write(np.str(np.round(meanann,2))+' '+np.str(np.round(stdevann,3))+'\n')
    
    # get monthly means and standard deviation
    file1.write('monthly means and standard deviations \n')
    file1.write('month    mean    sd  \n')
    
    mean_mon_data.coord('latitude').guess_bounds()
    mean_mon_data.coord('longitude').guess_bounds()
    grid_areas2 = iris.analysis.cartography.area_weights(mean_mon_data)
    tempcube=mean_mon_data.collapsed(['latitude','longitude'], 
                                iris.analysis.MEAN,weights=grid_areas2)
    meanmon=tempcube.data
    
    # get monthly average using grid areas from year average
    # to calculate standard deviation
    jan_avg=jan_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    feb_avg=feb_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    mar_avg=mar_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    apr_avg=apr_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    may_avg=may_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    jun_avg=jun_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    jul_avg=jul_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    aug_avg=aug_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    sep_avg=sep_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    oct_avg=oct_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    nov_avg=nov_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
    dec_avg=dec_slice.collapsed(['latitude','longitude'],iris.analysis.MEAN,weights=grid_areas)
   
    stdevmon=np.zeros(12)
   
    stdevcube=jan_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[0]=stdevcube.data
    stdevcube=feb_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[1]=stdevcube.data
    stdevcube=mar_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[2]=stdevcube.data
    stdevcube=apr_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[3]=stdevcube.data
    stdevcube=may_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[4]=stdevcube.data
    stdevcube=jun_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[5]=stdevcube.data
    stdevcube=jul_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[6]=stdevcube.data
    stdevcube=aug_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[7]=stdevcube.data
    stdevcube=sep_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[8]=stdevcube.data
    stdevcube=oct_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[9]=stdevcube.data 
    stdevcube=nov_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[10]=stdevcube.data 
    stdevcube=dec_avg.collapsed(['time'],iris.analysis.STD_DEV)
    stdevmon[11]=stdevcube.data
    
    for i in range(0,12):
        file1.write(np.str(i)+' '+np.str(np.round(meanmon[i],2))+' '+np.str(np.round(stdevmon[i],3))+'\n')
  
    
    file1.close()
    
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
   
    mean_mon_data.coord('latitude').guess_bounds() 
    mean_mon_data.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(mean_mon_data)
    temporal_mean = mean_mon_data.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
   
   
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




modelname='HadCM3'   # HadCM3 or MRI-CGCM2.3
exptnamein=['e280','eoi400']
#fieldnamein=(['pr','evap','psl','hfss','clt','tas','ts',
#              'rsdt','rsut','rlut','rsus','rlus','rsds','rlds',
#              'rsuscs','rsdscs','rldscs','sic','tos'] )
fieldnamein=['pr','tas']
#exptnamein=['e280']

for i in range(0,len(exptnamein)):
    for field in range(0,len(fieldnamein)):

        fieldnameout=fieldname.get(fieldnamein[field])
        exptnameout=exptname.get(exptnamein[i])

        if modelname=='HadCM3':
            fieldin=fieldname.get(fieldnamein[field])
        else:
            fieldin=fieldnamein[field]
                
        regrid_data(fieldin,fieldnameout,exptnamein[i],exptnameout,i+1,modelname)

#sys.exit(0)
