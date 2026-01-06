#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Monday June 3rd 2019

#@author: earjcti
#
# Manuel has asked if I can extract some fields from Kanhu's experiments so 
# that he can do some averages on them.  He would like data from the 
# xkrax experiment which has time varying vegetation, greenhouse gases varying 
#
# Because he is comparing to GNIP he would like data from years 1957-2014 
# as these are the years where we have observations.  This is j57-k14 I think
#
#
########################################################
# other notes are

import os
import numpy as np
import scipy as sp
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import cf_units as unit
import sys

#####################################
def extract_fields(filestart,fileoutstart,startyear,endyear,varnamein):



    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    
    
    # loop over months
    print(varnamein)
    
    moncube=iris.cube.CubeList([])
    for mon in range(0,len(monthnames)):
    #for mon in range(9,12):
  
        allcubes=iris.cube.CubeList([])
        #loop over years
        for year in range(startyear,endyear+1):
            yearuse=year-1800
           
            if year < 1800:
                extrause='??'
            else:
                if year < 1900:
                    extrause='i'
                else:
                    if year < 2000:
                        extrause='j'
                        yearuse=yearuse-100
                    else:
                        if year < 2100:
                            extrause='k'
                            yearuse=yearuse-200
                   
           
            stringyear=np.str(yearuse).zfill(2)
      
            
            filename=filestart+extrause+stringyear+monthnames[mon]+'.nc'
            print(filename)
           
           
            cube=iris.load_cube(filename,varnamein)
            
            
          
            u = unit.Unit('months since 0800-12-01 00:00:00',
                                  calendar=unit.CALENDAR_360_DAY)
            cube.coord('t').units=u
           
            cube.coord('t').points=year-startyear+1
           
            if varnamein=='Stash code = 338':
               
                cube.coord('level-1').rename('zdim')
               
              
                cube16o= (cube.extract(iris.Constraint(zdim=1.))+
                          cube.extract(iris.Constraint(zdim=4.))+
                          cube.extract(iris.Constraint(zdim=7.))+
                          cube.extract(iris.Constraint(zdim=10.)))
                cube18o= (cube.extract(iris.Constraint(zdim=2.))+
                          cube.extract(iris.Constraint(zdim=5.))+
                          cube.extract(iris.Constraint(zdim=8.))+
                          cube.extract(iris.Constraint(zdim=11.)))
                cube=((cube18o / cube16o)-2005.2E-6)/2005.2E-9
                
                cube.rename('d18o')
                cube.units='unknown'
                
             
            
            if varnamein=='TOTAL PRECIPITATION RATE     KG/M2/S':
                cube=(cube.extract(iris.Constraint(surface=0.000000)))
                
            if varnamein=='TEMPERATURE AT 1.5M':
                cube=(cube.extract(iris.Constraint(ht=-1.000000)))
               
           
            allcubes.append(cube)
           
        
           
        #make sure the metadata on all cubes are the same
        equalise_attributes(allcubes)
        unify_time_units(allcubes)
        for i in range(1,len(allcubes)):
            allcubes[i].coord('t').attributes=allcubes[0].coord('t').attributes
        
        
        
      
        catcube=allcubes.concatenate_cube()
       
       
        nc,ny,nx=np.shape(catcube)
        if nc != endyear-startyear+1:
            print('the cube has not been concatenated correctly')
            print('we should have',endyear-startyear+1,'times')
            print('we have',nc,'times')
            sys.exit(0)
         
        # find mean across cube dimension
        tempcube = catcube.collapsed('t', iris.analysis.MEAN)
        meancube=iris.util.new_axis(tempcube, 't')
        meancube.coord('t').points=mon+1
        meancube.coord('t').bounds=None
       
       
        # append to cube containing all months
        moncube.append(meancube)
       
       
    # unifty attributes on cubes
    equalise_attributes(moncube)
    unify_time_units(moncube)
    for i in range(1,len(moncube)):
            moncube[i].coord('t').attributes=moncube[0].coord('t').attributes
   
    
    outcube=moncube.concatenate_cube()
    
    
    return(outcube)
  
    
   

##########################################################
# main program

# this is regridding where all results are in a single file
# create a dictionary with the long field names in and the field names we want
# we are also using dictionaries so that we only have to change timeperiod name
# when rerunning

linuxwin='l'
            	
filestart={"l":'/nfs/hera1/earjcti/um/netcdf/xkrax_netcdf/xkraxa@pd',
           "w":'C:/Users/julia/OneDrive/WORK/DATA/TEMPORARY/xkraxa@pd'}

startyear=1958
endyear=2014

fileoutstart='output/'

fieldnames=['TOTAL PRECIPITATION RATE     KG/M2/S','TEMPERATURE AT 1.5M','Stash code = 338']
#fieldnames=['TOTAL PRECIPITATION RATE     KG/M2/S']

allcubes=[]

cube1=extract_fields(filestart.get(linuxwin),fileoutstart,startyear,endyear,fieldnames[0])
allcubes.append(cube1)
if len(fieldnames) >=2:
    cube2=extract_fields(filestart.get(linuxwin),fileoutstart,startyear,endyear,fieldnames[1])
    allcubes.append(cube2)
if len(fieldnames) >=3:
    cube3=extract_fields(filestart.get(linuxwin),fileoutstart,startyear,endyear,fieldnames[2])
    allcubes.append(cube3)
if len(fieldnames) >=4:
    print('you need to set up accessing more cubes')
    sys.exit(0)


fileout=fileoutstart+'data_xkrax.nc'        
iris.save(allcubes,fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
  
#sys.exit(0)
