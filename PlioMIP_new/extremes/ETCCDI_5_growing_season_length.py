#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia

We are looking at ETCCDI Climate change indicies.  This program will write
indices 5 (growing season length to a file.  This is:

GSL, Growing season length: Annual (1st Jan to 31st Dec in Northern Hemisphere (NH), 1st July to 30th June in Southern Hemisphere (SH)) count between first span of at least 6 days with daily mean temperature TG>5oC and first span after July 1st (Jan 1st in SH) of 6 days with TG<5oC.

Let TGij be daily mean temperature on day i in year j. Count the number of days between the first occurrence of at least 6 consecutive days with:

TGij > 5oC.

and the first occurrence after 1st July (1st Jan. in SH) of at least 6 consecutive days with:

TGij < 5oC. 
"""
import numpy as np
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import sys


def get_HCM3_year_data_NH(filestart,year):
    """
    reads in the mean daily temperature for the year and puts it in 
    a single cube
    """
    meanTcubelist = iris.cube.CubeList([])
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    for month in months:
        filename = filestart + np.str(year).zfill(2) + month + '.nc'
        # load in data
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp'))
        cube = iris.load(filename, constraints=variable_constraint)
        maxTcube = cube[0]

        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_1'))
        cube = iris.load(filename, constraints=variable_constraint)
        minTcube = cube[0]

        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_2'))
        cube = iris.load(filename, constraints=variable_constraint)
        meanTcube = cube[0]

        # check you have got maxT, meanT and minT in correct order        
        if np.max(maxTcube.data) < np.max(meanTcube.data):
            print('cubes not in right order')
            sys.exit(0)

        if np.max(meanTcube.data) < np.max(minTcube.data):
            print('cubes not in right order2')
            sys.exit(0)

        meanTcubelist.append(meanTcube-273.15)

    equalise_attributes(meanTcubelist)
    temporarycube = meanTcubelist.concatenate_cube()
    meanTyearcube = iris.util.squeeze(temporarycube)
  
    return meanTyearcube
   
def get_HCM3_year_data_SH(filestart,year):
    """
    reads in the mean daily temperature for the year and puts it in 
    a single cube
    However for the SH a year will go from July to june
    """
    meanTcubelist = iris.cube.CubeList([])
    months = ['jl','ag','sp','ot','nv','dc','ja','fb','mr','ar','my','jn',]
    for monthno, month in enumerate(months):
        if monthno <=5:  
            filename = filestart + np.str(year).zfill(2) + month + '.nc'
        else:
            filename = filestart + np.str(year+1).zfill(2) + month + '.nc'
        print(filename)
        # load in data
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp'))
        cube = iris.load(filename, constraints=variable_constraint)
        maxTcube = cube[0]

        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_1'))
        cube = iris.load(filename, constraints=variable_constraint)
        minTcube = cube[0]

        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_2'))
        cube = iris.load(filename, constraints=variable_constraint)
        meanTcube = cube[0]

        # check you have got maxT, meanT and minT in correct order        
        if np.max(maxTcube.data) < np.max(meanTcube.data):
            print('cubes not in right order')
            sys.exit(0)

        if np.max(meanTcube.data) < np.max(minTcube.data):
            print('cubes not in right order2')
            sys.exit(0)

        meanTcubelist.append(meanTcube-273.15)

    equalise_attributes(meanTcubelist)
    temporarycube = meanTcubelist.concatenate_cube()
    meanTyearcube = iris.util.squeeze(temporarycube)
  
    return meanTyearcube
   



def calc_grow_seas(meanTcube,nhshind):
    """
    note nhshind = 1.0 for nh and -1.0 for sh

    input is a cube of yearly temperature data.  We want to find the growing
    season length.  Do this as follows:

    1. find the first span of at least 6 days with temperature > 5degC
    2. find the first span after this with 6 days with temperature < 5degC
    """

    temporary_cube = meanTcube.collapsed('t',iris.analysis.SUM)
    growing_seas_arr = np.zeros(np.shape(temporary_cube.data))

    for j, lat in enumerate(meanTcube.coord('latitude').points):
        if lat * nhshind >= 0:
            for i, lon in enumerate(meanTcube.coord('longitude').points):
                timeseries = meanTcube[:, j, i].data

                growstart = -99
                growend = -99

                for day, temp in enumerate(timeseries):
                    if growstart < 0 and temp > 5.0 and day < 355 :
                       # print(day)
                        # is this start of growing seas
                        if (timeseries[day+1] > 5.0 and 
                            timeseries[day+2] > 5.0 and
                            timeseries[day+3] > 5.0 and 
                            timeseries[day+4] > 5.0 and
                            timeseries[day+5] > 5.0):
                            growstart = day
                            growend = 360
                        
                    if (growstart>=0 and growend == 360 
                        and temp < 5.0 and 180 < day < 355):
                        # is this the end of the growing season
                        # only check if end of growing season has not already
                        # been found
                     #    print('j',day, len(timeseries))
                         if (timeseries[day+1] < 5.0 and 
                            timeseries[day+2] < 5.0 and
                            timeseries[day+3] < 5.0 and 
                            timeseries[day+4] < 5.0 and
                            timeseries[day+5] < 5.0):
                            growend = day
                growing_seas_arr[j,i] = growend - growstart

    growing_seas_cube = temporary_cube.copy(data=growing_seas_arr)

    return growing_seas_cube

def get_HadCM3_diagnostics(expt, extra):
    """
    gets the diagnostics (frost days, summer days, icing days tropical nights
    from HadCM3)
    """
    filestart = '/nfs/hera1/earjcti/um/' + expt + '/pb/' + expt + 'a@pb' + extra
  
    for year in range(99, 100):
        meanTcube  = get_HCM3_year_data_NH(filestart, year)
        NH_growing_seas_len_cube = calc_grow_seas(meanTcube, 1.0)

        meanTcube  = get_HCM3_year_data_SH(filestart, year)
        SH_growing_seas_len_cube = calc_grow_seas(meanTcube, -1.0)

        growing_seas_cube = NH_growing_seas_len_cube.copy(np.maximum(NH_growing_seas_len_cube.data, SH_growing_seas_len_cube.data))

        growing_seas_cube.long_name = 'length of growing season (days)'
        growing_seas_cube.units = None
        
        outfile = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                   expt + '/' + extra + '_' + 
                   'diag5_' + np.str(year) + '.nc')
        iris.save(growing_seas_cube, outfile, netcdf_format="NETCDF3_CLASSIC")


##########################################################
def read_data(expt,extra,startyear,endyear):
    """
    reads in the data for each year, finds the sum and returns
    """
    filestart = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/'

    growcubelist = iris.cube.CubeList([])
   
    for year in range(startyear,endyear):
        filename = (filestart + expt + '/diag5/' + extra + 
                    '_diag5_'+ np.str(year) + '.nc')

        growcubelist.append(iris.load_cube(filename,'length of growing season (days)'))
        
    equalise_attributes(growcubelist)
    allgrowcube = growcubelist.merge_cube()
    meangrowcube = allgrowcube.collapsed('t',iris.analysis.MEAN)
   
  
    return meangrowcube


def plot_anom(cube,ocn_mask, expt, cntl, cube_cntl, cube_expt):
    """
    this will plot the anomaly cube between the pliocene and the pi
    """

    fieldname = 'len growing seas  '
    if ocn_mask == 'y':
       maskcube=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
      
       cube.data.mask = (maskcube.data - 1.0) * (-1.0)
       cube_cntl.data.mask = (maskcube.data - 1.0) * (-1.0)
       cube_expt.data.mask = (maskcube.data - 1.0) * (-1.0)
      
    plt.subplot(2,2,1)
    qplt.contourf(cube_cntl, levels=np.arange(100,320,20), extend='both')
    plt.gca().coastlines()
    plt.title(fieldname + ':' +  TIME.get(cntl))

    plt.subplot(2,2,2)
    qplt.contourf(cube_expt, levels=np.arange(100,320,20), extend='both')
    plt.gca().coastlines()
    plt.title(fieldname + ':' +  TIME.get(expt))

    plt.subplot(2,2,3)
    cmapname = 'Reds'
    vals = np.arange(0,30,5)
    qplt.contourf(cube,cmap=cmapname, levels=vals,extend='both')
    plt.gca().coastlines()
    plt.title(fieldname + ': ' + TIME.get(expt) + '-' +TIME.get(cntl))

    plt.subplot(2,2,4)
    qplt.contourf((cube / cube_cntl) * 100.,cmap=cmapname, 
                  levels=vals,extend='both')
    plt.gca().coastlines()
    plt.title('percentage change: ' + TIME.get(expt) + '-' +TIME.get(cntl))


    plt.tight_layout()
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag5/' + 
               'grow_seas_len_' + expt + '-' + cntl)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

       
def plot_growing_season(expt, cntl, nyears):
    """
    plots the growing season for plio and pi and difference between them
    """
    grow_expt_cube = read_data(expt,'o',0, nyears)

    grow_cntl_cube = read_data(cntl,'o',0, nyears)

    grow_anom_cube = grow_expt_cube - grow_cntl_cube
    plot_anom(grow_anom_cube,'y', expt,cntl,
              grow_cntl_cube, grow_expt_cube)
  
  
    
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'
TIME = {'tenvj' : 'mPlio', 'xozza' : 'PI', 'xozzb' : 'Plio'}


# this is for if you actually want to get the diagnostics.
# ie write them to the file /nfs/hera1/earjcti/PLIOMIP2/..ETCCDI....diags 5  
#if MODELNAME == 'HadCM3':
#    get_HadCM3_diagnostics('xozza','o')
 

# plot info about growing season.  a) growing season_pi, growing_season_plio
# and diff between them

plot_growing_season('tenvj','xozza',99)  # tenvj plio, xozza control, xozzb control
