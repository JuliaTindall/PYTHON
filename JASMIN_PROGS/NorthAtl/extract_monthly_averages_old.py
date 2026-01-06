#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

#Created on Monday June 3rd 2019

#@author: earjcti
#
# Jonathan has asked if I can extract some fields from Max's experiments so 
# Matt Jones can run these through his programs
#
# they would like monthly data averaged over 9-11ka and 5-7ka
#
# update:  30.10.2019 Would now also like 1000 year timesices
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

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
    
#####################################
def extract_fields(filestart, varnamein, exptlist, offset):


    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    filetype='a@pd'

   
    if varnamein=='Stash code = 338' or  varnamein=='Stash code = 321':    
        filetype='a@pc'
    
    # loop over months
    print(varnamein)
    
    moncube=iris.cube.CubeList([])
    for mon in range(0,len(monthnames)):
    #for mon in range(0,2):
        print(monthnames[mon])
        allcubes=iris.cube.CubeList([])
        if varnamein=='Stash code = 338':
            all16ocubes=iris.cube.CubeList([])
            all18ocubes=iris.cube.CubeList([])

        #loop over years
        for expt in range(0,len(exptlist)):
            for year in range(STARTYEAR, ENDYEAR + 1):
                
                print(filestart,exptlist[expt])
                print(filetype, EXTRA, year, monthnames[mon])
                filename=(filestart + exptlist[expt] + 
                          '/' + exptlist[expt] + filetype + EXTRA + 
                          np.str(year) + monthnames[mon] + '.nc')
                

               
                cube=iris.load_cube(filename,varnamein)          
                u = unit.Unit('months since 0800-12-01 00:00:00',
                                  calendar=unit.CALENDAR_360_DAY)
                cube.coord('t').units=u
                cube.coord('t').points=((year - STARTYEAR + 1)+
                            expt*(ENDYEAR - STARTYEAR + 1))
           
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

                    all16ocubes.append(cube16o)
                    all18ocubes.append(cube18o)


             
            
                if (varnamein=='TOTAL PRECIPITATION RATE     KG/M2/S' or
                   varnamein=='EVAPORATION FROM SEA (GBM)   KG/M2/S' or    
                   varnamein=='EVAP FROM SOIL SURF -AMOUNT KG/M2/TS' or
                   varnamein=='TRANSPIRATION RATE           KG/M2/S'):
                    cube=(cube.extract(iris.Constraint(surface=0.000000)))
                
                if (varnamein=='TEMPERATURE AT 1.5M' or 
                    varnamein=='RELATIVE HUMIDITY AT 1.5M'):
                    cube=(cube.extract(iris.Constraint(ht=-1.000000)))

                if (varnamein=='EVAP FROM CANOPY - AMOUNT   KG/M2/TS'):
                     cube=(cube.extract(iris.Constraint(level275=-1.000000)))

                if varnamein!='Stash code = 338':
                    allcubes.append(cube)


        ##################################
        # process monthly average for non d18op field

        if varnamein!='Stash code = 338':
           #make sure the metadata on all cubes are the same
           #print(allcubes)
            equalise_attributes(allcubes)
            unify_time_units(allcubes)
            for i in range(1,len(allcubes)):
                allcubes[i].coord('t').attributes=(
                    allcubes[0].coord('t').attributes)

            catcube=allcubes.concatenate_cube()
            print('concat cube is',np.shape(catcube))
       
            nc,ny,nx=np.shape(catcube)
            if nc != (ENDYEAR - STARTYEAR + 1)*(len(exptlist)):
                print('the cube has not been concatenated correctly')
                print('we should have', ENDYEAR - STARTYEAR + 1,'times')
                print('we have',nc,'times')
                sys.exit(0)
         
            # find mean across cube dimension
            tempcube = catcube.collapsed('t', iris.analysis.MEAN)
            meancube=iris.util.new_axis(tempcube, 't')
            meancube.coord('t').points=mon+1
            meancube.coord('t').bounds=None

        else: # do for d18op
           # sort out 16o cube
            equalise_attributes(all16ocubes)
            unify_time_units(all16ocubes)
            for i in range(1,len(all16ocubes)):
                all16ocubes[i].coord('t').attributes=(
                    all16ocubes[0].coord('t').attributes)

            catcube16o=all16ocubes.concatenate_cube()
       
            nc,ny,nx=np.shape(catcube16o)
            if nc != (ENDYEAR - STARTYEAR + 1)*(len(exptlist)):
                print('the 16o cube has not been concatenated correctly')
                print('we should have', ENDYEAR - STARTYEAR + 1,'times')
                print('we have',nc,'times')
                sys.exit(0)
            # sort out 18o cube
            equalise_attributes(all18ocubes)
            unify_time_units(all18ocubes)
            for i in range(1,len(all18ocubes)):
                all18ocubes[i].coord('t').attributes=(
                    all18ocubes[0].coord('t').attributes)

            catcube18o=all18ocubes.concatenate_cube()
       
            nc,ny,nx=np.shape(catcube18o)
            if nc != (ENDYEAR - STARTYEAR + 1)*(len(exptlist)):
                print('the 18o cube has not been concatenated correctly')
                print('we should have', ENDYEAR - STARTYEAR + 1,'times')
                print('we have',nc,'times')
                sys.exit(0)
         
            # find mean across cube dimension
            tempcube16o = catcube16o.collapsed('t', iris.analysis.MEAN)
            tempcube18o = catcube18o.collapsed('t', iris.analysis.MEAN)

            tempcube=((tempcube18o / tempcube16o)-2005.2E-6)/2005.2E-9
            # add average offset
            tempcube.data=tempcube.data+offset
            # replace nan with 0
            tempcube.data=np.where(
                tempcube16o.data==0.,0.0,tempcube.data) 
            tempcube.rename('d18o')
            tempcube.units='unknown'
                   
            meancube=iris.util.new_axis(tempcube, 't')
            meancube.coord('t').points=mon+1
            meancube.coord('t').bounds=None
        # end of stash code 338 check loop

       
       
        # append to cube containing all months
        moncube.append(meancube)
        

    print('end of months')       
    # unifty attributes on cubes
    equalise_attributes(moncube)
    unify_time_units(moncube)
    for i in range(1,len(moncube)):
            moncube[i].coord('t').attributes=moncube[0].coord('t').attributes
   
    
    outcube=moncube.concatenate_cube()
    


    if varnamein=='TOTAL PRECIPITATION RATE     KG/M2/S': 
        outcube.data=outcube.data * 60.*60.*24.
        outcube.long_name='TOTAL PRECIPITATION RATE    MM/DAY'
        outcube.units='mm/day'
            
    if varnamein=='EVAPORATION FROM SEA (GBM)   KG/M2/S':  
        outcube.data=outcube.data * 60.*60.*24.
        outcube.long_name='EVAPORATION FROM SEA    MM/DAY'
        outcube.units='mm/day'  
    if varnamein=='EVAP FROM SOIL SURF -AMOUNT KG/M2/TS': 
        outcube.data=outcube.data * 60.*60.*24. / 1800.
        outcube.long_name='EVAP FROM SOIL SURF   MM/DAY'
        outcube.units='mm/day'
    if varnamein=='TRANSPIRATION RATE           KG/M2/S': 
        outcube.data=outcube.data * 60.*60.*24.
        outcube.long_name='TRANSPIRATION RATE    MM/DAY'
        outcube.units='mm/day'
    if varnamein=='EVAP FROM CANOPY - AMOUNT   KG/M2/TS':
        outcube.data=outcube.data * 60.*60.*24. / 1800.
        outcube.long_name='EVAP FROM CANPOPY    MM/DAY'
        outcube.units='mm/day'
    
    return(outcube)
  
    
def extract_all_avg(filestart, exptnames, offset):
    """
    this will provide the loop that will extract the fieldnames one by one
    input:  filestart,  a list of experiments
    returns:  a cube containing monthly averaged data for each of the fields
    """
        
    cubelist = iris.cube.CubeList([])

    for field in range(0, len(FIELDNAMES)):
        cube = extract_fields(filestart, FIELDNAMES[field], exptnames, offset)
        cubelist.append(cube)
            
    
    return cubelist

##########################################################
# main program

# this is regridding where all results are in a single file
# create a dictionary with the long field names in and the field names we want
# we are also using dictionaries so that we only have to change timeperiod name
# when rerunning

LINUXWIN = 'j' # l-linux, w=windows, j jasmin
            	
FILESTART = {"l":'/nfs/hera1/earjcti/um/netcdf/xkrax_netcdf/xkraxa@pd',
             "w":'C:/Users/julia/OneDrive/WORK/DATA/TEMPORARY/',
             "j":'/home/users/jctindall/umoutput/BAS_timeslices/'}


STARTYEAR=20
if LINUXWIN == 'w':
    ENDYEAR = 31
else:
    ENDYEAR=50


EXTRA='v'
EXPT_11_9=['xlubl','xlubk','xlubj']
EXPT_7_5=['xlubf','xlubg','xlubh']

EXPTNAME = {"0ka" : "xluba",
            "1ka" : "xlubb",
            "2ka" : "xlubc",
            "3ka" : "xlubd",
            "4ka" : "xlube",
            "5ka" : "xlubf",
            "6ka" : "xlubg",
            "7ka" : "xlubh",
            "8ka" : "xlubi",
            "9ka" : "xlubj",
            "10ka" : "xlubk",
            "11ka" : "xlubl",
            "12ka" : "xlubm",
            "13ka" : "xlubn",
            "14ka" : "xlubo",
            "15ka" : "xlubp",
            "16ka" : "xlubq",
            "17ka" : "xlubr",
            "18ka" : "xlubs",
            "19ka" : "xlubt",
            "20ka" : "xlubu",
            "21ka" : "xlubv"
            }


# based on Lambeck et al 2014 table S3.  assuming 140m of sea level drop equates to
# 1permil of ocean d18o (because their LGM sea level was 140m).
D18O_OFFSET = {"0ka" : 0.0,
               "1ka" : 0.0,
               "2ka" : 0.0,
               "3ka" : 0.01,
               "4ka" : 0.01,
               "5ka" : 0.02,
               "6ka" : 0.02,
               "7ka" : 0.04,
               "8ka" : 0.09,
               "9ka" : 0.19,
               "10ka" : 0.29,
               "11ka" : 0.39,
               "12ka" : 0.46,
               "13ka" : 0.5,
               "14ka" : 0.55,
               "15ka" : 0.7,
               "16ka" : 0.77,
               "17ka" : 0.82,
               "18ka" : 0.86,
               "19ka" : 0.9,
               "20ka" : 0.96,
               "21ka" : 1.0
               }

FILEOUTSTART='modeloutput/'

FIELDNAMES=['TOTAL PRECIPITATION RATE     KG/M2/S',
            'TEMPERATURE AT 1.5M',
            'Stash code = 338',
            'EVAPORATION FROM SEA (GBM)   KG/M2/S',
            'RELATIVE HUMIDITY AT 1.5M',
            'EVAP FROM SOIL SURF -AMOUNT KG/M2/TS',
            'EVAP FROM CANOPY - AMOUNT   KG/M2/TS',
            'TRANSPIRATION RATE           KG/M2/S']
#FIELDNAMES=['Stash code = 338']

########################################

# do expt 11-9ka

#allcubes = extract_all_avg(FILESTART.get(LINUXWIN), EXPT_11_9, 0.29)


#fileout = FILEOUTSTART + 'data_11-9ka.nc'        
#iris.save(allcubes, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=2.0E20)

##########################################
# do individual years

PERIODS = ['21ka']
for time in range(0, len(PERIODS)):
    exptlist = [(EXPTNAME.get(PERIODS[time]))]
    offset = (D18O_OFFSET.get(PERIODS[time]))
    timecubes = extract_all_avg(FILESTART.get(LINUXWIN), exptlist, offset)
    print('timecubes is',timecubes)
    fileout=FILEOUTSTART + 'data_' + PERIODS[time] + '.nc'        
    iris.save(timecubes, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=2.0E20)

    

  
#sys.exit(0)
