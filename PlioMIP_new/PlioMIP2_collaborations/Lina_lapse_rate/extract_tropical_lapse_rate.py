#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#created May 2025
#
#@author: earjcti
#"""
#
#  Lina wants the tropical lapse rate
#  I am going to send
#  1. Netcdf file with Pliocene and preindustrial temperatures averaged over
#     30 years
#  2. png file plotted of avg lapse rate with pressure
#  3. txt file of what is in the png file
#
#
#

#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import iris
from iris.cube import CubeList
#import xlwt
#from xlwt import Workbook
import os
import matplotlib.cm as cm
from matplotlib.colors import Normalize
from iris.cube import CubeList
import sys


def get_avg(exptname,startyear,endyear):
    """
    gets the average data
    """
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    filestart = ('/nfs/hera1/earjcti/um/'+exptname+'/pcpd/' + exptname + 
                 'a#pc00000')
    
    cubes = CubeList([])
    for year in range(startyear,endyear):
        for month in months:
            filename = filestart + str(year) + month + '+.nc'
            cube = iris.load_cube(filename,'TEMPERATURE ON PRESSURE LEVELS')
            cubes.append(cube)

    iris.util.equalise_attributes(cubes)
    allcubes = cubes.concatenate_cube()
    meancube = allcubes.collapsed('t',iris.analysis.MEAN)
    
    return meancube


def get_tropical_lapse_rate(exptname,cntlname,startyear,endyear):
    """
    assuming we already have the mean averaged file
    """
    filename='avg_temp_on_pressure_levels'+str(startyear)+'_'+str(endyear)+'.nc'
    plio_cube=iris.load_cube(filename,
                             'TEMPERATURE ON PRESSURE LEVELS xqbwd mPWP')
    pi_cube = iris.load_cube(filename,
                             'TEMPERATURE ON PRESSURE LEVELS xqbwc PI')

    # get tropics cubes
    tropics_constraint = iris.Constraint(latitude=lambda cell: -20 < cell < 20)
    plio_tropics_cube = plio_cube.extract(tropics_constraint)
    pi_tropics_cube = pi_cube.extract(tropics_constraint)

    # get weighted average
    plio_tropics_cube.coord('latitude').guess_bounds()
    plio_tropics_cube.coord('longitude').guess_bounds()
    plio_tropics_grid_areas = iris.analysis.cartography.area_weights(
        plio_tropics_cube)
    plio_tropics_mean_cube = plio_tropics_cube.collapsed(
        ["latitude","longitude"],
        iris.analysis.MEAN,weights=plio_tropics_grid_areas)
   

    pi_tropics_cube.coord('latitude').guess_bounds()
    pi_tropics_cube.coord('longitude').guess_bounds()
    pi_tropics_grid_areas = iris.analysis.cartography.area_weights(
        pi_tropics_cube)
    pi_tropics_mean_cube = pi_tropics_cube.collapsed(
        ["latitude","longitude"],
        iris.analysis.MEAN,weights=pi_tropics_grid_areas)
    
    fig,ax=plt.subplots()
    #ax.set_yscale('log')
    ax.set_yticks(plio_tropics_mean_cube.coord('p').points)
    ax.set_ylim(1000,100)
    ax.plot(plio_tropics_mean_cube.data-273.15,
             plio_tropics_mean_cube.coord('p').points,
             label='pliocene')
    ax.plot(pi_tropics_mean_cube.data-273.15,
             pi_tropics_mean_cube.coord('p').points,
             label='pi')
    plt.xlabel('Temperature deg C')
    plt.ylabel('Pressure hPa')
    ax.legend()
    fileout = 'Temp_on_pressure_levels_'+str(startyear)+'_'+str(endyear)+'.png'

    plt.savefig(fileout)
    
    f=open('Temp_on_pressure_levels_'+str(startyear)+'_'+str(endyear)+'.txt',
           'w')
    f.write('pressure,Pliocene temperature,PI temperature\n')
    for i,press in enumerate(pi_tropics_mean_cube.coord('p').points):
        f.write(str(press) + ',' + str(plio_tropics_mean_cube.data[i]-273.15)
                + ',' + str(pi_tropics_mean_cube.data[i]-273.15) + '\n')

    f.close()
            
                
    
    

#################
# MAIN PROGRAM
################

exptname = 'xqbwd'
cntlname = 'xqbwc'
startyear = 3900
endyear=4000


#mean_expt_cube = get_avg(exptname,startyear,endyear)
#mean_expt_cube.long_name='TEMPERATURE ON PRESSURE LEVELS xqbwd mPWP'
#mean_cntl_cube = get_avg(cntlname,startyear,endyear)
#mean_cntl_cube.long_name='TEMPERATURE ON PRESSURE LEVELS xqbwc PI'

#fileout = 'avg_temp_on_pressure_levels'+str(startyear)+'_'+str(endyear)+'.nc'
#iris.save([mean_expt_cube,mean_cntl_cube],fileout)

# get tropical lapse rate
get_tropical_lapse_rate(exptname,cntlname,startyear,endyear)
