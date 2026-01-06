#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 3rd 2019


#@author: earjcti
#
# This program will regrid some of the data that is needed for PLIOMIP2.  
# We will put 100 year average fields onto a 1deg X 1deg standard grid
# it can be used where experiments have been uploaded with 100 years in 
# one file
# 
# it can currently do MIROC4 and COSMOS

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
import iris.coord_categorisation
import cf_units as unit
import sys



#####################################
def regrid_data(fieldnamein,fieldnameout,exptnamein,exptnameout,filename,modelname,linux_win,fielduse):

   
    
    print('moodelname is',modelname)
    print('filename is',filename)
        
        
    # outfile
    if linux_win=='l':
        outstart=('/nfs/hera1/earjcti/regridded/'+modelname+'/'+exptnameout+'.'+
        fieldnameout+'.')
    else:
        outstart=('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
              +modelname+'\\'+exptnameout+'.'+fieldnameout+'.')
   
    print(modelname)
    if modelname=='EC-Earth3.1': # all fields in one file
        allcube=iris.load(filename)
        ncubes=len(allcube)
        for i in range(0,ncubes):
            if allcube[i].var_name==fielduse:
                cube=allcube[i]
       
    else: # one field in a file
        cube=iris.load_cube(filename)
   
    # now regrid the cube onto a 1X1 grid (we will first try regridding the raw data)
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'
    
    # do not need to regrid UofTdata
    if modelname =='UofT-CCSM4':
        regridded_cube=cube
    else:
        print('regridding')
        cubegrid=iris.load_cube('one_lev_one_deg.nc')
        regridded_cube=cube.regrid(cubegrid,iris.analysis.Linear())
    
   
    refdate='days since 0800-01-01 00:00:00'
    # for cosmos
    if modelname=='COSMOS':
        # cosmos data is in a strange time coordinate line yyyymmdd
        # we need to convert it to days since reference time
        origpoints=regridded_cube.coord('time').points
        npoints=len(origpoints)
        print(npoints)
       
        print(origpoints)
        yeararr=np.zeros(npoints)
        montharr=np.zeros(npoints)
        dayarr=np.zeros(npoints)
        daydecimal=np.zeros(npoints)
        dayssinceref=np.zeros(npoints)
        for i in range(0,npoints):
            origstr=str(origpoints[i])
            yeararr[i]=origstr[:][0:4]
            montharr[i]=origstr[:][4:6]
            dayarr[i]=origstr[:][6:8]
            daydecimal[i]=origstr[:][8:]
            dayssinceref[i]=dayssinceref[i-1]+dayarr[i]+daydecimal[i]-daydecimal[i-1]
        # subtract 1 from days since reference date (as reference date will be 1st Jan)
        dayssinceref=dayssinceref-1
       
        
        regridded_cube.coord('time').points=dayssinceref
        #  end of COSMOS loop
    
    # for EC-Earth3.1
    if modelname=='EC-Earth3.1':
        # convert from hours to days
        origpoints=regridded_cube.coord('time').points
        newpoints=origpoints/24.
        regridded_cube.coord('time').points=newpoints
        refdate='days since 2390-01-01 00:00:00'
        # regrid to mm/day
        if fieldnamein=='pr':
            regridded_cube.data=regridded_cube.data * 60. *60. *24. 
            cube.data=cube.data* 60. *60. *24. 
            regridded_cube.name='Total precipitation'
            regridded_cube.long_name='Total precipitation'
            regridded_cube.units='mm/day'
        # convert to celcius
        if fieldnamein=='ts':
            regridded_cube.convert_units('Celsius')
            cube.convert_units('Celsius')
        
    if modelname=='UofT':
        # if precipitation convert from m/s to mm/day
        if fieldnamein=='pr':
            regridded_cube.data=regridded_cube.data * 60. *60. *24. *1000.
            cube.data=cube.data* 60. *60. *24. *1000.
            regridded_cube.name='Total precipitation'
            regridded_cube.long_name='Total precipitation'
            regridded_cube.units='mm/day'
        # convert to celcius
        if fieldnamein=='ts':
            regridded_cube.convert_units('Celsius')
            cube.convert_units('Celsius')
           
        # we need to add the missing time coordinate
        points=(np.arange(0,1200)*30)+15. # go for middle of month
        u = unit.Unit('days since 0800-01-01 00:00:00',
               calendar=unit.CALENDAR_360_DAY) # put as 360 day calendar because of the way
                                               # the data was sent.
       
        regridded_cube.add_dim_coord(iris.coords.DimCoord(points, 
                standard_name='time', long_name='time', 
                var_name='time', 
                #units=iris.unit.Unit('days since 0800-01-01 00:00:00',calendar='360_day'),
                units=u,
                bounds=None,
                coord_system=None, circular=False),0)
    else: 
        
         regridded_cube.coord('time').units=refdate
         
       # end of Uof T loop
   
    
  
    # add auxillary coordinates month and year

    
    iris.coord_categorisation.add_month(regridded_cube, 'time', name='month')
    iris.coord_categorisation.add_year(regridded_cube, 'time', name='year')
    
   
     # mean the cubes over the month and year dimension
     # standard deviation for each month
    mean_mon_data=regridded_cube.aggregated_by('month',iris.analysis.MEAN)
    mean_mon_data.long_name=fieldnameout
   
    sd_mon_data=regridded_cube.aggregated_by('month',iris.analysis.STD_DEV)
    sd_mon_data.long_name=fieldnameout
    
    mean_year_data=regridded_cube.aggregated_by('year',iris.analysis.MEAN)
    mean_year_data.long_name=fieldnameout
    
    print('j2 mean mon',np.shape(mean_mon_data))
    print('j3 mean year',np.shape(mean_year_data))
    
    mean_data=mean_mon_data.collapsed('time',iris.analysis.MEAN)
    mean_data.long_name=fieldnameout
    sd_data=mean_year_data.collapsed('time',iris.analysis.STD_DEV)
    sd_data.long_name=fieldnameout
    
    print('j4 all mon',np.shape(mean_data))
    print('j5 all year',np.shape(sd_data))
    
    # extract monthly data from cubes
 
    jan_slice = regridded_cube.extract(iris.Constraint(month='Jan'))
    feb_slice = regridded_cube.extract(iris.Constraint(month='Feb'))
    mar_slice = regridded_cube.extract(iris.Constraint(month='Mar'))
    apr_slice = regridded_cube.extract(iris.Constraint(month='Apr'))        
    may_slice = regridded_cube.extract(iris.Constraint(month='May'))
    jun_slice = regridded_cube.extract(iris.Constraint(month='Jun'))
    jul_slice = regridded_cube.extract(iris.Constraint(month='Jul'))
    aug_slice = regridded_cube.extract(iris.Constraint(month='Aug'))
    sep_slice = regridded_cube.extract(iris.Constraint(month='Sep'))
    oct_slice = regridded_cube.extract(iris.Constraint(month='Oct'))
    nov_slice = regridded_cube.extract(iris.Constraint(month='Nov'))
    dec_slice = regridded_cube.extract(iris.Constraint(month='Dec'))
    
    
   
    
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
    stdevcube=tempcube.collapsed(['time'],iris.analysis.STD_DEV)
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
    plt.subplot(2,2,1) # global mean from each year
    #subcube=subcube_mean_mon.copy(data=)  # set up structure of subcube
   
    if modelname=='COSMOS' or modelname=='UofT' or modelname=='EC-Earth3.1':
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
    else:
        cube.coord('latitude').bounds
        cube.coord('longitude').bounds
    grid_areas = iris.analysis.cartography.area_weights(cube)
    newcube=cube.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
    nt=len(newcube.data)
    nyears=nt/12
    
    for i in range(0,nyears):
        tstart=i*12
        tend=(i+1)*12
        plotdata=newcube.data[tstart:tend]
        plt.plot(plotdata,color='r')
        

  
    # global mean from average
   
    grid_areas = iris.analysis.cartography.area_weights(mean_mon_data)
    temporal_mean = mean_mon_data.collapsed(['latitude','longitude'], iris.analysis.MEAN,weights=grid_areas)
   
    plt.plot(temporal_mean.data,color='b',label='avg')
    plt.title('globavg '+fieldnamein)
    plt.legend()
    
  
    
    
    # check at 30N
    plt.subplot(2,2,2)
    
    bounds=cube.coord('latitude').bounds
    nbounds=cube.coord('latitude').nbounds
    nbounds,dummy=np.shape(bounds)
    for i in range(0,nbounds):
        if (bounds[i,0]>=32. >=bounds[i,1] or bounds[i,0]<=32. <bounds[i,1]): 
            index=i
   
    subcube=cube[:,index,:]
    cube_avg_30N=subcube.collapsed(['longitude'],iris.analysis.MEAN)
    
    for i in range(0,nyears):
        tstart=i*12
        tend=(i+1)*12
        plotdata=cube_avg_30N.data[tstart:tend]
        plt.plot(plotdata,color='r')
    
    #mean at 30N
    slice_30N= mean_mon_data.extract(iris.Constraint(latitude=32))
    mean_30N=slice_30N.collapsed(['longitude'], iris.analysis.MEAN)
    
   
    plt.plot(mean_30N.data,color='b',label='avg')
    plt.title('average at 30N by month')
    plt.legend()
    plt.show()
    plt.close()
    

  


##########################################################
# main program

filename=' '
linux_win='l'
modelname='COSMOS' # MIROC4m  COSMOS UofT EC-Earth3.1

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

COSMOS_FIELDS ={"pr" : "TotalPrecip",
        "ts" : "SurfaceTemp",
        "sic" : "SeaIceAreaFraction"
        }

ECearth_FIELDS ={"pr" : "totp",
        "ts" : "tas",
        "sic" : "SeaIceAreaFraction"
        }

ECearth_EXPT={"Eoi400": "mPlio",
              "E280":"PI"
              }
# this is regridding where all results are in a single file
#fieldnamein=['pr','ts','sic']
#exptnamein=['Eoi400','E280']

fieldnamein=['pr','ts']
exptnamein=['Eoi400','E280']

#fieldnamein=['ts']
#exptnamein=['Eoi400','E280']
if linux_win=='l':
    filestart='/nfs/hera1/pliomip2/data/'
else:
    filestart='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'

for expt in range(0,len(exptnamein)):
    for field in range(0,len(fieldnamein)):

        fieldnameout=fieldname.get(fieldnamein[field])
        exptnameout=exptname.get(exptnamein[expt])

        if modelname == 'MIROC4m':
            filename=filestart+modelname+'/'
            fielduse=MIROC_FIELDS.get(fieldnamein[field])
            print(fielduse)
            filename=(filename+fielduse+
                      '/MIROC4m_'+exptnamein[expt]+'_Amon_'+fielduse+'.nc')
        if modelname == 'COSMOS':
            if linux_win=='l':
                filename=filestart+'/AWI/COSMOS/'
                filename=filename+exptnamein[expt]+'/'
            else:
                filename=filestart+'/COSMOS/'
            fielduse=COSMOS_FIELDS.get(fieldnamein[field])
            filename=(filename+exptnamein[expt]+'.'+fielduse+
                      '_CMIP6_name_'+fieldnamein[field]+
                      '_2650-2749_monthly_mean_time_series.nc')
        if modelname == 'UofT':
            if linux_win=='l':
                filename=filestart+modelname+'/'
                filename=filename+'UofT-CCSM4/'+exptnamein[expt]+'/Amon/'
            else:
                filename=filestart+'UofT-CCSM4\\'+exptnamein[expt]+'\\'
            fielduse=MIROC_FIELDS.get(fieldnamein[field])
            filename=(filename+fielduse+
                      '_Amon_'+exptnamein[expt]+
                      '_'+modelname+'-CCSM4_gr.nc')
        if modelname=='EC-Earth3.1':
            filename=(filestart+modelname+'/'+
                      ECearth_EXPT.get(exptnamein[expt])+
                      '.'+modelname+'.surface.nc')
            fielduse=ECearth_FIELDS.get(fieldnamein[field])
            
       
        
        
        regrid_data(fieldnamein[field],fieldnameout,exptnamein[expt],exptnameout,filename,modelname,linux_win,fielduse)

#sys.exit(0)