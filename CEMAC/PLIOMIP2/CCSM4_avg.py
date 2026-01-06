#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 29th


#@author: earjcti
"""
This program will average all of the CCSM4 models.
We will average the data file (text files) and also the mean average temperature file
"""

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
import iris.quickplot as qplt
import iris.plot as iplt
import cf_units as unit
from iris.experimental.equalise_cubes import equalise_attributes
import sys


#########################################################################
# stuff for doing text file is here  
def get_text_data(filename):
    """
    gets the text data for each filename
    returns globalmean, monthly mean, latitudinal mean
    """
    f=open(filename,"r")
    f1=f.readlines()
    f2 = [x.replace('\n', '') for x in f1]
            
    # get the means according to their position in the file
    all_mean_sd=f2[2]
    all_mon_mean_sd=f2[5:5+12]
    all_lat_mean_sd=f2[20:20+180]
           
    # extract global mean
    meanglob,sd=all_mean_sd.split(',')
    
   
    monmeans = np.zeros(12)
    latmeans = np.zeros(180)
    lats = np.zeros(180)
    # extract monthly means
    for x in all_mon_mean_sd:
        mon,mean,sd=x.split(',')
        monmeans[int(mon)-1]=float(mean)
            
            
    # extract latitude means
    for x in all_lat_mean_sd:
        lat,mean,sd=x.split(',')
        latss=int(float(lat)+89.5) # convert latitude to a subscript
               
        if mean != ' --' and mean != '--':
            latmeans[latss]=float(mean) # stores latitudinal means
        else:
            latmeans[latss]=np.nan
            
        lats[latss]=lat # stores latitudes
      
   
    return meanglob, monmeans, latmeans, lats
   
def writeout_textdata(mean_global, mean_mon, mean_lat, lats, expt):
    """
    writeout the text data to a file
    replace stdev with -999
    """
    
    textout = FILESTART + OUTMODEL + '/' + expt + '.' + FIELDNAME + '.data.txt'
    file1 =  open(textout, "w")

    # write out to a file
    file1.write('global annual mean and standard deviation\n')
    file1.write('------------------------------------------\n')
   
    # write out global temperautre
    file1.write(np.str(np.round(mean_global, 2))+', -99.99\n')

    # get monthly means and standard deviation
    file1.write('monthly means and standard deviations \n')
    file1.write('----------------------------------------')
    file1.write('month    mean    sd  \n')

    for i in range(0, 12):
        file1.write(np.str(i+1)+', '+np.str(np.round(mean_mon[i], 2))+', -99.99\n')

    # get latitudinal means and standard deviation
    file1.write('zonal means and standard deviations \n')
    file1.write('----------------------------------------\n')
    file1.write('latitude    mean    sd  \n')
    for i in range(0, len(mean_lat)):
        file1.write(np.str(lats[i])+', '+np.str(np.round(mean_lat[i], 2))+', -99.99\n')

    file1.close()
    

def main_avg_text(expt):
    """
    loop over all the models and extract mean, monthlymeans, latitudinal means
    average all the means
    write out to a file in the format of the input.  (Put standard deviation to -999.999)
    """
    
    for i, model in enumerate(MODELNAMES):
        filename = FILESTART + model + '/' + expt + '.' + FIELDNAME + '.data.txt'
        
        globmean, monmeans, latmeans, lats = get_text_data(filename)
        
        if i == 0:
            allmeans = np.zeros(NMODELS)
            allmonmeans = np.zeros((NMODELS, len(monmeans)))
            alllatmeans = np.zeros((NMODELS, len(latmeans)))
        
        allmeans[i] = globmean
        allmonmeans[i, :] = monmeans
        alllatmeans[i, :] = latmeans
        
    mean_global = np.mean(allmeans)
    mean_mon = np.mean(allmonmeans, axis=0)
    mean_lat = np.mean(alllatmeans, axis=0)
    
    writeout_textdata(mean_global, mean_mon, mean_lat, lats, expt)
 
###############################################
## stuff for doing netcdf file is here
def main_avg_netcdf(expt):
    """
    loop over all the models and extract the global average netcdf file
    average all the means
    write out to a file in the format of the input.  
    """ 
    all_cubes=iris.cube.CubeList([])     
    for i, model in enumerate(MODELNAMES):
        print(i)
        filename = FILESTART + model + '/' + expt + '.' + FIELDNAME + '.allmean.nc'
        cube = iris.load_cube(filename)
        modelcube = resort_coords(cube, i)
        modelcube.data=modelcube.data.astype('float32') 
        all_cubes.append(modelcube)
    
    
    iris.experimental.equalise_cubes.equalise_attributes(all_cubes)
  
    cat_cubes = all_cubes.concatenate_cube()
    meancube = cat_cubes.collapsed(['model_level_number'], iris.analysis.MEAN)  
    
    fileout = FILESTART + OUTMODEL + '/' + expt + '.' + FIELDNAME + '.allmean.nc'
    iris.save(meancube, fileout, netcdf_format = 'NETCDF3_CLASSIC', fill_value = 2.0E20)
  
    
def resort_coords(cube,levelno):
    """
    this will make all the dimensions of the cube match.  They will all be
    longitude, latitude, level-no (ie 1 for first model, 2 for second model...)
    
    input is the cube and the level number
    output is the cube with the new dimensions
    """
    
    for coord in cube.coords():        
        name=coord.standard_name
        if name !='latitude' and name!='longitude':
            if name==None:
                if coord.long_name==None:
                    cube.remove_coord(coord.var_name)
                else:
                    cube.remove_coord(coord.long_name)
            else:
                cube.remove_coord(name)
                
    for coord in cube.coords():   # now this will be longitude or latitude
        coord.points=coord.points.astype('float32') 
        coord.var_name = coord.standard_name
        coord.long_name = coord.standard_name
       
       
     
    newcube=iris.util.new_axis(cube)
    newcube.add_dim_coord(iris.coords.DimCoord(levelno, 
            standard_name='model_level_number', long_name='model', 
            var_name='model', 
            units=None,
            bounds=None,
            coord_system=None, circular=False),0) 
   
    # this will make sure cell_methods match and that cubes can
    # be concatenated
    newcube.cell_methods = None
    newcube.rename('tas')
    
        
    return newcube      
    
    
def main():
    """ 
    main program
    1. average the text files for each of the models and writeout
    2. average the netcdf files from each of the models and writeout
    """
    
    for i, expt in enumerate(EXPTNAMES):
        avgtext = main_avg_text(expt)
        avgnetcdf = main_avg_netcdf(expt)
    
    
   

##########################################################
# fixed constants
        
LINUX_WIN='w'
if LINUX_WIN == 'w':
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'#
else:
    FILESTART = ' '

MODELNAMES=['CCSM4-1deg', 'CCSM4-2deg','CCSM4-UoT']
NMODELS = len(MODELNAMES)
OUTMODEL = 'CCSM4-avg'

FIELDNAME='NearSurfaceTemperature'
EXPTNAMES=['EOI400','E280']
#EXPTNAMES=['EOI400']


main()