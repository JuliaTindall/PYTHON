#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia


This will plot the weather at a particular location like we might see on the internet.

For each month we will get average Tmax, average Tmin and total precipitation
and plot them all

"""
import numpy as np
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import sys


def get_HCM3_year_data(filestart,year):
    """
    reads in the maximum  and minimum temperature for the year and puts it in 
    a single cube
    """
    maxTcubelist = iris.cube.CubeList([])
    minTcubelist = iris.cube.CubeList([])
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    for month in months:
        filename = filestart + np.str(year).zfill(2) + month + '.nc'
        cubes=iris.load(filename)
        temperature_cubes = iris.cube.CubeList([])
        max_temp = []
        for cube in cubes:
            if (cube.var_name == 'temp' or cube.var_name == 'temp_1' 
                or cube.var_name == 'temp_2'):
               try:
                   cube.coord('t_1').rename('t')
               except:
                   pass
               cube_reg_1 = cube.extract(iris.Constraint(latitude = LAT_REQ,
                                                       longitude = LON_REQ))
               cube_reg = cube_reg_1.collapsed('t',iris.analysis.MEAN)
               temperature_cubes.append(cube_reg)
               max_temp.append(np.max(cube_reg.data))
     
        maxindex = np.argmax(max_temp)
        minindex = np.argmin(max_temp)
        indexes = np.argsort(max_temp)
        minTcube = temperature_cubes[indexes[0]]
        meanTcube = temperature_cubes[indexes[1]]
        maxTcube = temperature_cubes[indexes[2]]
    
        # check you have got maxT, meanT and minT in correct order        
        if np.max(maxTcube.data) < np.max(meanTcube.data):
            print('cubes not in right order')
            sys.exit(0)

        if np.max(meanTcube.data) < np.max(minTcube.data):
            print('cubes not in right order2')
            sys.exit(0)

        maxTcubelist.append(maxTcube-273.15)
        minTcubelist.append(minTcube-273.15)

        Tmax = []
        Tmin = []
        for cube in maxTcubelist:
            Tmax.append(cube.data[0])
        for cube in minTcubelist:
            Tmin.append(cube.data[0])
   
    return np.asarray(Tmax),np.asarray(Tmin)
   

def get_HadCM3_diagnostics(expt, extra):
    """
    gets the diagnostics (frost days, summer days, icing days tropical nights
    from HadCM3)
    """
    filestart = '/nfs/hera1/earjcti/um/' + expt + '/pb/' + expt + 'a@pb' + extra

    count=0
    summaxT = np.zeros(12)
    summinT = np.zeros(12)
    for year in range(0, 10):
        yearuse = np.str(year).zfill(2)
        maxT, minT  = get_HCM3_year_data(filestart, year)
        
        summaxT = summaxT + maxT
        summinT = summinT + minT
        count=count+1
    
    avg_maxT = summaxT / count
    avg_minT = summinT / count
    print(avg_maxT, avg_minT)
    
    plt.plot(avg_maxT)
    plt.plot(avg_minT)
    plt.show()
    sys.exit(0)

        
##########################################################
def read_data(expt,extra,startyear,endyear):
    """
    reads in the data for each year, finds the sum and returns
    """
    filestart = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/'

    frostcubelist = iris.cube.CubeList([])
    summercubelist = iris.cube.CubeList([])
    icingcubelist = iris.cube.CubeList([])
    tropcubelist = iris.cube.CubeList([])
 
    for year in range(startyear,endyear):
        filename = (filestart + expt + '/diag1-4/' + extra + 
                    '_diag1-4_'+ np.str(year) + '.nc')

        frostcubelist.append(iris.load_cube(filename,'number of frost days'))
        summercubelist.append(iris.load_cube(filename,'number of summer days'))
        icingcubelist.append(iris.load_cube(filename,'number of icing days'))
        tropcubelist.append(iris.load_cube(filename,
                                           'number of tropical nights'))

    equalise_attributes(frostcubelist)
    allfrostcube = frostcubelist.merge_cube()
    meanfrostcube = allfrostcube.collapsed('t',iris.analysis.MEAN)
   
    equalise_attributes(summercubelist)
    allsummercube = summercubelist.merge_cube()
    meansummercube = allsummercube.collapsed('t',iris.analysis.MEAN)
    meansummercube.units = None
  
    equalise_attributes(icingcubelist)
    allicingcube = icingcubelist.merge_cube()
    meanicingcube = allicingcube.collapsed('t',iris.analysis.MEAN)
   
    equalise_attributes(tropcubelist)
    alltropcube = tropcubelist.merge_cube()
    meantropcube = alltropcube.collapsed('t',iris.analysis.MEAN)

    return (meanfrostcube, meansummercube, meanicingcube,meantropcube)
  
def plot_anom(cube,field,ocn_mask, expt, cntl, cube_cntl, cube_expt):
    """
    this will plot the anomaly cube between the pliocene and the pi
    """
    fieldname = {'frost' : 'number of frost days',
                 'summer': 'number of summer days',
                 'icing': 'number of icing days',
                 'tropical_nights' : 'number of tropical nights'}

    if ocn_mask == 'y':
       maskcube=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
      
       cube.data.mask = (maskcube.data - 1.0) * (-1.0)
       cube_cntl.data.mask = (maskcube.data - 1.0) * (-1.0)
       cube_expt.data.mask = (maskcube.data - 1.0) * (-1.0)
      
    plt.subplot(2,2,1)
    qplt.contourf(cube_cntl, levels=np.arange(0,300,20), extend='max')
    plt.gca().coastlines()
    plt.title(fieldname.get(field) + ':' +  TIME.get(cntl))

    plt.subplot(2,2,2)
    qplt.contourf(cube_expt, levels=np.arange(0,300,20), extend='max')
    plt.gca().coastlines()
    plt.title(fieldname.get(field) + ':' +  TIME.get(expt))

    plt.subplot(2,2,3)
    if field == 'frost' or field == 'icing':
        cmapname = 'Blues_r'
        vals = np.arange(-40,5,5)
    else:
        cmapname = 'Reds'
        vals = np.arange(0,65,5)
    qplt.contourf(cube,cmap=cmapname, levels=vals,extend='both')
    plt.gca().coastlines()
    plt.title(fieldname.get(field) + ': ' + TIME.get(expt) + '-' +TIME.get(cntl))
    plt.subplot(2,2,4)
    qplt.contourf((cube / cube_cntl) * 100.,cmap=cmapname, 
                  levels=vals,extend='both')
    plt.gca().coastlines()
    plt.title('percentage change: ' + TIME.get(expt) + '-' +TIME.get(cntl))


    plt.tight_layout()
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag1-4/' + 
               field + '_' + expt + '-' + cntl)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()
    
    
def plot_n_extremes(expt, exptextra, cntl, cntlextra, nyears):
    """
    plots the number of days that are: 1. frost days, 2. summer days,
    3. icing days, 4. tropical nights and how it changes between plio and cntl
    """
    (frost_expt_cube, summerday_expt_cube, 
    icing_expt_cube, tropnight_expt_cube) = read_data(expt,exptextra,0, nyears)

    (frost_cntl_cube, summerday_cntl_cube, 
    icing_cntl_cube, tropnight_cntl_cube) = read_data(cntl,cntlextra,0, nyears)

    frost_anom_cube = frost_expt_cube - frost_cntl_cube
    plot_anom(frost_anom_cube,'frost','y', expt,cntl,
              frost_cntl_cube, frost_expt_cube)
  
    summer_anom_cube = summerday_expt_cube - summerday_cntl_cube
    plot_anom(summer_anom_cube,'summer','y', expt,cntl,
              summerday_cntl_cube, summerday_expt_cube)
  
    icing_anom_cube = icing_expt_cube - icing_cntl_cube
    plot_anom(icing_anom_cube,'icing','y', expt,cntl,
              icing_cntl_cube, icing_expt_cube)
  
    trop_anom_cube = tropnight_expt_cube - tropnight_cntl_cube
    plot_anom(trop_anom_cube,'tropical_nights','y', expt,cntl,
              tropnight_cntl_cube, tropnight_expt_cube)
  
    
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'

TIME = {'tenvj' : 'mPlio', 'xozza' : 'PI', 'xozzb' : 'Plio','tenvs':'E560',
        'xozzm' : 'E560'}

LAT_REQ = 52.5
LON_REQ = 0.0

# this is for if you actually want to get the diagnostics.

Tmax_mon, Tmin_mon = get_HadCM3_diagnostics('xozzm','w')

# plot map of number of days that are 'extreme' according to the ETCCDI 1-4
# criteria
# plot anomaly from plio and pi

('xozzm','w','xozza','o',99)  # tenvj plio, xozza control, xozzb control
 

