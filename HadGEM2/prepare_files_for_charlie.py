#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 18 14:13:50 2019

#@author: earjcti1
#
# This program will extract the fields that charlie asked for
#
#
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
from iris.experimental.equalise_cubes import equalise_attributes
import sys
import warnings

warnings.filterwarnings("ignore")

def simplify_cube(cube):
    """
    gets cube and makes sure dimensions are longitude, latitude surface and
    t. 
    """    
    for coord in cube.coords():
        if coord.var_name == 'level275':
            coord.var_name = 'surface'
    
    cube.coord('surface').points = 0.0
    cube.coord('surface').units = 'm'
    cube.coord('surface').attributes = None
    
    cube.data = np.where(cube.data > 1.0E10, 0., cube.data)
    return cube

#####################################

 
def extract_fields(filestart,expt,startyear,endyear,varnamein):

    # load required cubes
    #cubes=iris.load(filename)
    #print(cubes)
    #sys.exit(0)
    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    
    
    # loop over years
    
   
    for year in range(startyear,endyear):
        stringyear=np.str(year).zfill(2)
        for mon in range(0,len(monthnames)):
            filename=(filestart + stringyear + monthnames[mon] + '_temp.nc')
            print(filename,varnamein)
            
            cube = iris.load_cube(filename, varnamein)
            allcubes.append(cube)
       
    print(allcubes)
    #make sure the metadata on all cubes are the same
    equalise_attributes(allcubes)
    catcube=allcubes.concatenate()
        
   
    return catcube

     
   

##########################################################
# main program

# this is regridding where all results are in a single file
# create a dictionary with the long field names in and the field names we want
# we are also using dictionaries so that we only have to change timeperiod name
# when rerunning
            
FIELEXTRA = {"TOTAL PRECIPITATION RATE     KG/M2/S" : "pd",
		"SURFACE TEMPERATURE AFTER TIMESTEP" : "pd",
        "POTENTIAL TEMPERATURE (OCEAN)  DEG.C" : "pf",
      	}

SHORTNAME = {
		"SURFACE TEMPERATURE AFTER TIMESTEP" : "SurfaceTemperature",
        "TOTAL PRECIPITATION RATE     KG/M2/S" : "TotalPrecipitation",
        "POTENTIAL TEMPERATURE (OCEAN)  DEG.C" : "SST",
  	}


#FIELDNAMES = [ "TOTAL PRECIPITATION RATE     KG/M2/S",
#               "SURFACE TEMPERATURE AFTER TIMESTEP",
#        "POTENTIAL TEMPERATURE (OCEAN)  DEG.C",
#       	]

FIELDNAMES = ["POTENTIAL TEMPERATURE (OCEAN)  DEG.C"]
	       
#fieldname=["V COMPNT OF WIND ON PRESSURE LEVELS"]

expt = 'xkvjg'
linux_win='l'
startyear=50
endyear=100


filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+expt+'/'

allcubes=iris.cube.CubeList([])
      

for i, varnamein in enumerate(FIELDNAMES):
    varnameout = SHORTNAME.get(varnamein)
    if varnameout == 'SST':
        filestart2 = filestart + 'temp_data/' + expt + 'o@pfn'
    if varnameout == 'TotalPrecipitation':
        filestart2 = filestart + 'precip_data/' + expt + 'a@pdn'
    if varnameout == 'SurfaceTemperature':
        filestart2 = filestart + '/netcdf/pdfiles/' + expt + 'a@pdn'
   
    print(varnameout, filestart2)
    cube = extract_fields(filestart2,expt,startyear,endyear, varnamein)

    fileout=('/nfs/hera1/earjcti/um/HadGEM_data/'
             +expt+'/'+expt+ '_' + varnameout + '_timeseries.nc')

    print(cube)
    iris.save(cube, fileout, fill_value=2.0E20)
    
#sys.exit(0)
