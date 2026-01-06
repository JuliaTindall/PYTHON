#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia

This will count the number of days when the potential wet bulb temperature 
exceeds certain thresholds
20
25
30
35
and write to a file 
There will be one record in the file for each year
"""
import numpy as np
import iris
#from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import sys


def get_HCM3_year_data(filestart,year):
    """
    reads in the wet bulb temperature for the year and puts it in 
    a single cube
    """
    wetbulb_cubelist = iris.cube.CubeList([])
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    for month in months:
        filename = filestart + np.str(year).zfill(2) + month + '.nc'
        cube=iris.load_cube(filename,'WET BULB POTENTIAL TEMPERATURE')
        wetbulb_cubelist.append(cube)
        
    iris.util.equalise_attributes(wetbulb_cubelist)
    wetbulbcube = wetbulb_cubelist.concatenate_cube(wetbulb_cubelist)
   
    return wetbulbcube - 273.15

def count_wet_bulb_days(wetbulbcube, year):
    """
    input the wet bulb cube for the year
    output the number of days where the wet bulb temperature exceeds
           thresholds
    """
    wetbulbdata = wetbulbcube.data
    # find where wetbulb pot temp > vals
    wet_gt_20 = np.where(wetbulbdata > 20, 1.0, 0)
    wet_gt_25 = np.where(wetbulbdata > 25, 1.0, 0)
    wet_gt_30 = np.where(wetbulbdata > 30, 1.0, 0)
    wet_gt_35 = np.where(wetbulbdata > 35, 1.0, 0)

    # set up cubes
    wet_gt_20_cube = wetbulbcube.copy(data = wet_gt_20)
    wet_gt_25_cube = wetbulbcube.copy(data = wet_gt_25)
    wet_gt_30_cube = wetbulbcube.copy(data = wet_gt_30)
    wet_gt_35_cube = wetbulbcube.copy(data = wet_gt_35)

    # count how many days
    
    ndays_wet_gt_20_cube = wet_gt_20_cube.collapsed('t',iris.analysis.SUM)
    ndays_wet_gt_25_cube = wet_gt_25_cube.collapsed('t',iris.analysis.SUM)
    ndays_wet_gt_30_cube = wet_gt_30_cube.collapsed('t',iris.analysis.SUM)
    ndays_wet_gt_35_cube = wet_gt_35_cube.collapsed('t',iris.analysis.SUM)

    # set up for concatenation by year
    ndays_wet_gt_20_cube.coord('p').rename('year')
    ndays_wet_gt_25_cube.coord('p').rename('year')
    ndays_wet_gt_30_cube.coord('p').rename('year')
    ndays_wet_gt_35_cube.coord('p').rename('year')

    ndays_wet_gt_20_cube.coord('year').points = [year]
    ndays_wet_gt_25_cube.coord('year').points = [year]
    ndays_wet_gt_30_cube.coord('year').points = [year]
    ndays_wet_gt_35_cube.coord('year').points = [year]

    ndays_wet_gt_20_cube.remove_coord('t')
    ndays_wet_gt_25_cube.remove_coord('t')
    ndays_wet_gt_30_cube.remove_coord('t')
    ndays_wet_gt_35_cube.remove_coord('t')
    
    return (ndays_wet_gt_20_cube, ndays_wet_gt_25_cube, ndays_wet_gt_30_cube,
            ndays_wet_gt_35_cube)

def get_HadCM3_diagnostics(expt, extra):
    """
    gets the diagnostics (frost days, summer days, icing days tropical nights
    from HadCM3)
    """
    filestart = '/nfs/hera1/earjcti/um/' + expt + '/pb/' + expt + 'a@pb' + extra
    
    allcubes_wbgt20 = iris.cube.CubeList([])
    allcubes_wbgt25 = iris.cube.CubeList([])
    allcubes_wbgt30 = iris.cube.CubeList([])
    allcubes_wbgt35 = iris.cube.CubeList([])

    for year in range(1, 79):
        print(year)
        yearuse = np.str(year).zfill(2)
        wetbulbcube  = get_HCM3_year_data(filestart, year)

      
        (wet_gt_20_cube, 
         wet_gt_25_cube, 
         wet_gt_30_cube, 
         wet_gt_35_cube) = count_wet_bulb_days(wetbulbcube, year)
        

        allcubes_wbgt20.append(wet_gt_20_cube)
        allcubes_wbgt25.append(wet_gt_25_cube)
        allcubes_wbgt30.append(wet_gt_30_cube)
        allcubes_wbgt35.append(wet_gt_35_cube)
    
    iris.util.equalise_attributes(allcubes_wbgt20)
    print(allcubes_wbgt20[0],allcubes_wbgt20[1])
    wbgt20_cube = allcubes_wbgt20.concatenate_cube()
    wbgt20_cube.long_name='ndays with wetbulb pottemp gt 20degC'
    
    wbgt25_cube = allcubes_wbgt25.concatenate_cube()
    wbgt25_cube.long_name='ndays with wetbulb pottemp gt 25degC'

    wbgt30_cube = allcubes_wbgt30.concatenate_cube()
    wbgt30_cube.long_name='ndays with wetbulb pottemp gt 30degC'
   
    wbgt35_cube = allcubes_wbgt35.concatenate_cube()
    wbgt35_cube.long_name='ndays with wetbulb pottemp gt 35degC'
   
    
    cubelist = [wbgt20_cube, wbgt25_cube, wbgt30_cube, wbgt35_cube]

    outfile = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                   expt + '/wetbulb_thresholds_' + extra + '.nc')
    iris.save(cubelist, outfile, netcdf_format="NETCDF3_CLASSIC")

  
    
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'

TIME = {'tenvj' : 'mPlio', 'xozza' : 'PI', 'xozzb' : 'Plio','tenvs':'E560',
        'xozzm' : 'E560'}
# this is for if you actually want to get the diagnostics.
# ie write them to the file /nfs/hera1/earjcti/PLIOMIP2/..ETCCDI....diags 1-4  
get_HadCM3_diagnostics('xozzb','p')


