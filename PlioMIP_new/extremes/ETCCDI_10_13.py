#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia

We are looking at ETCCDI Climate change indicies.  This program will deal with indices 10-13 these are

10.  percentage of days when TN < 10th percentile for preindustrial (TN=daily minimum temperature)
11. percentage of days when TX < 10th percentile for preindustrial (TX=daily maximum temperature)
12. percentage of days when TN > 90th percentile
13. percentage of days when TX > 90th percentile

"""
import numpy as np
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
from iris.coords import DimCoord
import matplotlib.pyplot as plt
import sys
from iris.experimental.equalise_cubes import equalise_attributes
from os.path import exists



#######################################################################
#  STEP 1 WRITE ALL THE PERCENTILES TO A FILE   
def get_months_days(day):
    """
    To find percentiles we find all the temperatres for the 100 years which 
    are within 2 days. 
    So for 15th May, we would include temperatures 
    from 13th 14th 15th 16th 17th May
    """

    # find central month
    monthreq = np.int(np.floor(day / 30))
    dayreq = day - monthreq * 30

    # find day before
    daym1 = dayreq -1 
    monthm1 = monthreq
    if daym1 < 0:
        monthm1 = monthreq - 1
        if monthm1 < 0: monthm1 = monthm1 + 12
        daym1 = daym1 + 30


    # find two days before
    daym2 = dayreq -2 
    monthm2 = monthreq
    if daym2 < 0:
        monthm2 = monthreq - 1
        if monthm2 < 0: monthm2 = monthm2 + 12
        daym2 = daym2 + 30
    

    # find day after
    dayp1 = dayreq +1 
    monthp1 = monthreq
    if dayp1 >= 30:
        monthp1 = monthreq + 1
        if monthp1 >= 12: monthp1 = monthp1 - 12
        dayp1 = dayp1 - 30

    # find two days after
    dayp2 = dayreq +2 
    monthp2 = monthreq
    if dayp2 >= 30:
        monthp2 = monthreq + 1
        if monthp2 >= 12: monthp2 = monthp2 - 12
        dayp2 = dayp2 - 30

    monthsreq = [monthm2, monthm1, monthreq, monthp1, monthp2]
    daysreq = [daym2, daym1, dayreq, dayp1, dayp2]
    return [monthsreq, daysreq]

def get_temperatures(monthreq,daysreq):
    """
    get the temperatures that are appropriate for finding the percentiles
    for the day of interest
    """
    monthnames = {0:'ja',1:'fb',2:'mr',3:'ar',4:'my',5:'jn',6:'jl',
                  7:'ag',8:'sp',9:'ot',10:'nv',11:'dc'}
    if MIN_MAX == 'min':
         # this is minimum temperature
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_1'))
    if MIN_MAX == 'max':
         # this is maximum temperature
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp'))
    if MIN_MAX == 'mean':
         # this is maximum temperature
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_2'))

    temperature_cubes = iris.cube.CubeList([])
  
    for i, month in enumerate(monthreq):
        day = daysreq[i]
        for year in range(0,100):
            filename = (FILESTART + np.str(year).zfill(2) 
                            + monthnames.get(month) + '.nc')
            if exists(filename):

                cubelist = iris.load(filename, constraints=variable_constraint)
                cube = cubelist[0]
                monthprev = month
        
                tempinit = cube[day, :, :, :]
                temperature = iris.util.squeeze(tempinit)
                temperature_cubes.append(temperature)
            
    equalise_attributes(temperature_cubes)
    print(temperature_cubes)
    temperatures = temperature_cubes.merge_cube()
    print(daysreq,'merged')
   
    return temperatures
    
def percentiles(nperc, temperatures_cube):
    """
    get the percentiles
    in nperc ; the percentiles we want
       temperatures cube;  the temperatures for that day for which we
                           need to find the percentiles, 
                           DIMENSIONS: t latitude, longitude
    """
    data = temperatures_cube.data
    # cube shape for putting the new data
    cube_shape = temperatures_cube.collapsed('t', iris.analysis.MEAN)  
    percentile_array = np.zeros((len(nperc), 73, 96))
    for lat in range(0, 73):
        for lon in range(0, 96):
            # find the sorted temperature
            temperature = np.sort(data[:, lat, lon])
            ntemp = len(temperature)
            for i, perc in enumerate(nperc):
                index = np.int(np.rint(perc * ntemp / 100))
                percentile_array[i, lat, lon] = temperature[index]
                
    # put the percentiles into the cube  
    percentile_cubelist = iris.cube.CubeList([])
    for i in range(0, len(nperc)):
        cube = cube_shape.copy(data=percentile_array[i, :, :])
        if MIN_MAX == 'min':
            cube.long_name = np.str(nperc[i]) + 'th percentile of minimum temperatures'
        if MIN_MAX == 'max':
            cube.long_name = np.str(nperc[i]) + 'th percentile of maximum temperatures'
        if MIN_MAX == 'mean':
            cube.long_name = np.str(nperc[i]) + 'th percentile of maximum temperatures'
        cube.var_name = np.str(nperc[i]) + 'th_percentile'
        percentile_cubelist.append(cube)

    return percentile_cubelist
    

def get_HadCM3_percentiles(expt, extra):
    """
    gets the diagnostics percentiles for writing to a file
    """
    filestart = '/nfs/hera1/earjcti/um/' + expt + '/pb/' + expt + 'a@pb' + extra
  
    # THIS PART OF THE PROGRAM WILL CALCULATE THE PERCENTILES FOR EACH DAY   
    for day in range(0, 360):
        # find days that we need to get temperature.  This is the day required
        # two days before and two days afterwards
        monthsreq, daysreq = get_months_days(day)

        print(monthsreq,daysreq)
        # get all the minimum/ maximum temperatures that fall on one of the 
        # correct days
        temperatures_cube = get_temperatures(monthsreq,daysreq)
        print(temperatures_cube)
        # we should have 500 temperatures so we can find the percentiles
        percentiles_cubelist = percentiles([5, 10, 90, 95],temperatures_cube) 
      
        # write the percentiles to a file for that day.  
        print(np.str(day))
        filename = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                    EXPT + '/diag_10_13/daily_percentiles/' + 
                    MIN_MAX + '_temperature_' + np.str(day) + '.nc')

        iris.save(percentiles_cubelist, filename, 
                  netcdf_format="NETCDF3_CLASSIC")
###########################################################################
#  STEP 2 
#  READ ALL THE PERCENTILES FROM THE FILE (OBTAINED IN STEP 1) AND
#  SEE HOW MANY DAYS ARE MORE EXTREME  
def get_month_temperatures(month):      
    """
    get the temperatures over 100 years for this month
    """
    monthnames = {0:'ja',1:'fb',2:'mr',3:'ar',4:'my',5:'jn',6:'jl',
                  7:'ag',8:'sp',9:'ot',10:'nv',11:'dc'}
    # orig code that didn't quite work
   # if MIN_MAX == 'min':
   #      # this is minimum temperature
   #     variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == '#temp_1'))
   # if MIN_MAX == 'max':
   #      # this is maximum temperature
   #     variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp'))
   # if MIN_MAX == 'mean':
   #      # this is maximum temperature
   #     variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == '#temp_2'))

    temperature_cubes = iris.cube.CubeList([])
  
    for year in range(0,NYEARS):
        filename = (FILESTART + np.str(year).zfill(2) 
                            + monthnames.get(month) + '.nc')
        #if exists(filename):
        #cubelist = iris.load(filename, constraints=variable_constraint)
        #cube = cubelist[0]
        # put year in ht field for concatenation

        # new code
        cubes=iris.load(filename)
        cubelist = iris.cube.CubeList([])
        max_temp = []
        for cube in cubes:
            if (cube.var_name == 'temp' or cube.var_name == 'temp_1' 
                or cube.var_name == 'temp_2'):
               try:
                   cube.coord('t_1').rename('t')
               except:
                   pass
               cubelist.append(cube)
               max_temp.append(np.max(cube.data))
        maxindex = np.argmax(max_temp)
        minindex = np.argmin(max_temp)
        indexes = np.argsort(max_temp)
        if MIN_MAX == 'min':
            cubereq = cubelist[indexes[0]]
        if MIN_MAX == 'max':
            cubereq = cubelist[indexes[2]]
        if MIN_MAX == 'mean':
            cubereq = cubelist[indexes[1]]
       
        cubereq.coord('ht').points = year
        cubereq.coord('t').points = np.arange(0,30,1)
        temperature_cubes.append(cubereq)
            

    equalise_attributes(temperature_cubes)
    print(temperature_cubes)
    temperatures = temperature_cubes.concatenate_cube()
    return temperatures

def percentiles_from_file(day):
    """
    reads in the percentiles that we have obtained in step 1
    and returns
    """
    maxminind = {'min': 'minimum', 'max':'maximum'}
   
    filename = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                   CNTL + '/diag_10_13/daily_percentiles/' + 
                   MIN_MAX + '_temperature_' + np.str(day) + '.nc')
    perc5_cube = iris.load_cube(filename, '5th percentile of ' + 
                               maxminind.get(MIN_MAX) + ' temperatures')
    perc10_cube = iris.load_cube(filename, '10th percentile of ' + 
                               maxminind.get(MIN_MAX) + ' temperatures')
    perc90_cube = iris.load_cube(filename, '90th percentile of ' + 
                               maxminind.get(MIN_MAX) + ' temperatures')
    perc95_cube = iris.load_cube(filename, '95th percentile of ' + 
                               maxminind.get(MIN_MAX) + ' temperatures')

    return perc5_cube, perc10_cube, perc90_cube, perc95_cube

def nextr(temperatures_cube, percentile_cube, lt_gt, percentile, day):
    """
    here we find the number of temperatures that is more extreme than the
    percentile
    """      
    if lt_gt == 'less than':
        extreme_data = np.where(temperatures_cube.data <= percentile_cube.data, 
                                1.0, 0.0)
    if lt_gt == 'more than':
        extreme_data = np.where(temperatures_cube.data >= percentile_cube.data, 
                                1.0, 0.0)
    
    extreme_allyears = ((np.sum(extreme_data, axis=0) * 100)  /
                        np.size(extreme_data,axis=0))
    extreme_cube = percentile_cube.copy(data=extreme_allyears)
    name = ('pcent of years ' + lt_gt + ' ' + CNTL + 
            ' ' + np.str(percentile) + 'th percentile')
    extreme_cube.long_name = name
    extreme_cube.cell_methods = None
    extreme_cube.remove_coord('t')
    extreme_cube.attributes = None
    extreme_cube.var_name = 'pcent_extreme'
    time = DimCoord(day,standard_name = 'time', long_name = 't', 
                              var_name = None, units = 'day')
    print(extreme_cube)
    extreme_cube.add_aux_coord(time)
    ext3d_cube = iris.util.new_axis(extreme_cube, time)
   
    return ext3d_cube
       
def find_nextremes():
    """
    this will find the number of days which are more extreme than the
    percentiles (5th, 10th, 90th, 95th) in the control  (normally preindustrial)
    """  
    extreme5_cubelist = iris.cube.CubeList([])
    extreme10_cubelist =iris.cube.CubeList([])
    extreme90_cubelist = iris.cube.CubeList([])
    extreme95_cubelist = iris.cube.CubeList([])

    for month in range(0,12):
        month_temp_cube = get_month_temperatures(month)
        for day in range(0, 30):
            # get percentile cubes for day
            dayperc = (month *30) + day
            print('dayperc is',dayperc)
            (perc5_cube, perc10_cube, 
             perc90_cube, perc95_cube) = percentiles_from_file(dayperc)
      
             # get temperatures for this day in the year in expt
            Tcube = month_temp_cube[day, :, :, :]
           
            # find how many temperatures are more extreme than percentiles
            extreme5_cube = nextr(Tcube, perc5_cube, 'less than', 5, day)
            extreme10_cube = nextr(Tcube, perc10_cube, 'less than', 10, day)
            extreme90_cube = nextr(Tcube, perc90_cube, 'more than', 90, day)
            extreme95_cube = nextr(Tcube, perc95_cube, 'more than', 95, day)

            # put days in 0-360 format
            extreme5_cube.coord('time').points=[dayperc]
            extreme10_cube.coord('time').points=[dayperc]
            extreme90_cube.coord('time').points=[dayperc]
            extreme95_cube.coord('time').points=[dayperc]

            extreme5_cubelist.append(extreme5_cube)
            extreme10_cubelist.append(extreme10_cube)
            extreme90_cubelist.append(extreme90_cube)
            extreme95_cubelist.append(extreme95_cube)

    equalise_attributes(extreme5_cubelist)
    equalise_attributes(extreme10_cubelist)
    equalise_attributes(extreme90_cubelist)
    equalise_attributes(extreme95_cubelist)
   
  
    extreme5_year_cube = extreme5_cubelist.concatenate_cube()
    extreme10_year_cube = extreme10_cubelist.concatenate_cube()
    extreme95_year_cube = extreme95_cubelist.concatenate_cube()
    extreme90_year_cube = extreme90_cubelist.concatenate_cube()
   
    filename = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + EXPT + '/diag_10_13/percentage_more_extreme_than_T'+ MIN_MAX + '_' + CNTL + '_5th_percentile.nc'
    iris.save(extreme5_year_cube, filename ,netcdf_format="NETCDF3_CLASSIC")
  
    filename = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + EXPT + '/diag_10_13/percentage_more_extreme_than_T'+ MIN_MAX + '_'+  CNTL + '_10th_percentile.nc'
    iris.save(extreme10_year_cube, filename ,netcdf_format="NETCDF3_CLASSIC")
   
    filename = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + EXPT + '/diag_10_13/percentage_more_extreme_than_T'+ MIN_MAX + '_' + CNTL + '_90th_percentile.nc'
    iris.save(extreme90_year_cube, filename ,netcdf_format="NETCDF3_CLASSIC")
   
    filename = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + EXPT + '/diag_10_13/percentage_more_extreme_than_T'+ MIN_MAX + '_' + CNTL + '_95th_percentile.nc'
    iris.save(extreme95_year_cube, filename ,netcdf_format="NETCDF3_CLASSIC")
   
#########################################################
#STEP 3 plot information about extremes
def get_percentile(percentile):
    """
    gets the control temperature of percentile
    """
    maxminind = {'min': 'minimum', 'max':'maximum'}
   
    cubelist = iris.cube.CubeList([])
    for day in range(0, 360):
        filename = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                    CNTL + '/diag_10_13/daily_percentiles/' + 
                    MIN_MAX + '_temperature_' + np.str(day) + '.nc')
        perc_cube = iris.load_cube(filename, np.str(percentile) + 
                                    'th percentile of ' + 
                                    maxminind.get(MIN_MAX) + ' temperatures')
        perc_cube.coord('t').points=day
        perc_cube.coord('t').attributes=None
        perc_cube.coord('t').bounds=None
        cubelist.append(perc_cube)

    equalise_attributes(cubelist)
    percentiles_cube = cubelist.merge_cube()
    mean_percentile = percentiles_cube.collapsed('t',iris.analysis.MEAN)
    mean_percentile.convert_units('celsius')
   
    return mean_percentile

    
def plot_extremes_exceeded_annual(percentile, ocn_mask):
    """
    reads in the files showing when the extremes have been exceeded and find
    an annual averages them and plots (land only)
    """
    gtlt = {10 : ' < ', 5 : ' < ' ,90 : ' > ', 95  : ' > '}
    filename = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + EXPT 
                + '/diag_10_13/percentage_more_extreme_than_T' + MIN_MAX 
                + '_' + CNTL + '_' + np.str(percentile) 
                + 'th_percentile.nc')
    cube = iris.load_cube(filename)
   
    if ocn_mask == 'y':
        if EXPT == 'tenvj':
            filelsm = '/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrparm.mask.nc'
            cube = iris.load(filelsm)
            print(cube)
            exptlsmcube=iris.load_cube(filelsm,'LAND MASK (LOGICAL: LAND=TRUE)')
            cube_ann_avg.data.mask = (exptlsmcube.data - 1.0) * (-1.0)
    
    if percentile==10:
        valmin=0
        valmax=16
        valrange=2
    if percentile==90 or percentile==95:
        valmin=0
        valmax=55
        valrange=5

    cube_ann_avg.units = None
   # qplt.pcolormesh(cube_ann_avg,vmin=minval, vmax=maxval)
    qplt.contourf(cube_ann_avg,levels=np.arange(valmin,valmax,valrange),
                  extend='max')
    plt.gca().coastlines()
    plt.title('% days: T' + MIN_MAX + gtlt.get(percentile) + 
              TIME.get(CNTL) + ' ' + np.str(percentile)  + 
              ' percentile',fontsize=8)
   
def plot_temp(temp_cube, ocn_mask, title, vals):
    """
    plots a cube of temperatures 
    """   
  
    if ocn_mask == 'y':
        if EXPT == 'tenvj':
            filelsm = '/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrparm.mask.nc'
            cube = iris.load(filelsm)
            print(cube)
            exptlsmcube=iris.load_cube(filelsm,'LAND MASK (LOGICAL: LAND=TRUE)')
            temp_cube.data.mask = (exptlsmcube.data - 1.0) * (-1.0)
        
    qplt.contourf(temp_cube, levels=vals, extend='both')
    plt.gca().coastlines()
    plt.title(title,fontsize=8)
   
    
  
def plot_extremes(perc_low, perc_high):
    """
    we will do a 4 column plot.
    1. cntl 90th percentile
    2. average difference between highest and lowest percentiles
       (to get an idea of natural variability) 
percencentage of days in pliocene which exceed perc_high for pi
    3. difference between perc_high and perc low (to get an idea of 
       natural variability    
    4. percentage of days in pliocene which are lower than perc_low for pi
    """
    temp_percentile_low = get_percentile(perc_low)
    temp_percentile_high = get_percentile(perc_high)
    temp_perc_diff = temp_percentile_high - temp_percentile_low
    
    plt.subplot(2,2,1)
    title = ('ann_avg ' + np.str(perc_low) + 'percentile of T' + 
             MIN_MAX +  ' for' + TIME.get(CNTL))
    vals = np.arange(-45,45,15)
    plot_temp(temp_percentile_low,'y',title, vals )
    
    plt.subplot(2,2,2)
    title = (TIME.get(CNTL) + ' T' + MIN_MAX + ': ' + np.str(perc_high) + 'th -'
             + np.str(perc_low) + 'th percentile')
    vals = np.arange(0,21,3)
    plot_temp(temp_perc_diff,'y',title, vals)

    plt.subplot(2,2,3)
    plot_extremes_exceeded_annual(perc_low, 'y')
   
    plt.subplot(2,2,4)
    plot_extremes_exceeded_annual(perc_high, 'y')
    
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag10-13/' +
               EXPT + '_' + CNTL + np.str(perc_low) + '_' + np.str(perc_high) + 
               '_percentiles_annual_' + 'T' + MIN_MAX)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()
   
   
   
    
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'
MIN_MAX = 'max'
NYEARS = 100
TIME = {'xozza': 'PI', 'tenvj' : 'plio','xozzm':'E560'}
  

# this is for obtaining the percentiles used in climate change indices   
# STEP1
#EXPT= 'tenvj'
EXPT = 'xozzm'  # this is the experiments we are checking
CNTL = 'xozza'  # we are seeing how many of the days  in the experiment
                    # have temperatures more extreme than those in the control 
FILESTART = '/nfs/hera1/earjcti/um/'+EXPT+'/pb/'+EXPT+'a@pbw'
#get_HadCM3_percentiles(EXPT,'w')
#get_HadCM3_percentiles(CNTL,'o')

# STEP 2
# count the number of days with the temperature greater than all the percentiles
find_nextremes()

# STEP3
# plot information from the extremes.  Will need to use the files derived in step 2
plot_extremes(5,95)  # number in brackets is the percentile
