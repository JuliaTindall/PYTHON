#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Dec.2024
#@author: earjcti
#"""
#
# This will calculate an annual mean for surface evaporation (and also 
# surface runoff, subsurface runoff
#

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import iris
from iris.cube import CubeList
import sys

def get_avg_evap_mine():
    """
    gets the mean evaporation from my experiments that have been got for
    the pliomip2 database
    """        
    
    filestart = ('/nfs/hera1/earjcti/um/'+EXPT+'/evap/' + EXPT + '.evap.')

    cubes = CubeList([])
    for year in range(YEARSTART, YEARSTART+NYEARS):
      filename = filestart + str(year).zfill(3) + '.nc'
      print(filename)
      cube = iris.load_cube(filename)
      cubes.append(cube)

    allcubes = cubes.concatenate_cube()
    meancube = allcubes.collapsed('t',iris.analysis.MEAN)

    return meancube


def get_avg_runoff_mine(fieldname):
    """
    gets the mean surface and subsurface runoff from the raw files
    """        
    
    extraname = {
        "xozza" : "p",
        "xozzb" : "p",
        "xozzc" : "o",
        "xozzf" : "o",
        "tenvl" : "o",
        "tenvk" : "o",
        "tenvm" : "o"}
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']

    filestart = ('/nfs/hera1/earjcti/um/'+EXPT+'/pd/' + EXPT + 'a@pd' + 
                 extraname.get(EXPT))

    cubes = CubeList([])
    for year in range(YEARSTART, YEARSTART+NYEARS):
        print(year)
        for month in months:
            filename = filestart + str(year) + month + '.nc'
            cube = iris.load_cube(filename,fieldname)
            cubes.append(cube)

    iris.util.equalise_attributes(cubes)
    iris.util.unify_time_units(cubes)
    allcubes = cubes.concatenate_cube()
    print(allcubes)
    meancube = allcubes.collapsed('t',iris.analysis.MEAN)
    meandata=meancube.data

    print(meancube.data.mask)
    meancube_corr_data = np.ma.where(meancube.data.mask, -99999.,meancube.data)
    meancube_runoff=meancube.copy(data=meancube_corr_data)
    return meancube_runoff

#################
# MAIN PROGRAM
################





FIELD = 'evap' # evap or runoff
EXPT='tenvm'
YEARSTART=70
NYEARS=30

if FIELD == 'evap':
    #if EXPT == 'xozzb' or EXPT == 'xozzc' or EXPT == 'xozzf':
        avgcubes = get_avg_evap_mine()

if FIELD == 'runoff':
    #if EXPT == 'xozzb' or EXPT == 'xozzc' or EXPT == 'xozzf':
        avgsurf = get_avg_runoff_mine('SURFACE RUNOFF RATE          KG/M2/S')
        avgsubs = get_avg_runoff_mine('SUB-SURFACE RUNOFF RATE      KG/M2/S')
        avgcubes = [avgsurf,avgsubs]



fileout = ('/nfs/hera1/earjcti/um/'+EXPT+'/means/' + EXPT + '.' + FIELD + '.' + str(YEARSTART) + '-' + str(YEARSTART + NYEARS-1) + '.nc')

iris.save(avgcubes,fileout,fill_value=-99999.)

