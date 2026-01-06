
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia

We are looking at ETCCDI Climate change indicies.  This program will deal with indices 14-15 these are warm spell and cold spell duration index

14 WSDI, Warm spell duration index: Annual count of days with at least 6 consecutive days when TX > 90th percentile

Let TXij be the daily maximum temperature on day i in period j and let TXin90 be the calendar day 90th percentile centred on a 5-day window for the base period 1961-1990. Then the number of days per period is summed where, in intervals of at least 6 consecutive days:

TXij > TXin90

15 CSDI, Cold spell duration index: Annual count of days with at least 6 consecutive days when TN < 10th percentile

Let TNij be the daily maximum temperature on day i in period j and let TNin10 be the calendar day 10th percentile centred on a 5-day window for the base period 1961-1990. Then the number of days per period is summed where, in intervals of at least 6 consecutive days:


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


def percentiles_from_file():
    """
    reads in the percentiles that we have obtained in step 1
    and returns 12 cubes one for each month
    """
    if MIN_MAX == 'min':
        fieldname = '10th percentile of minimum temperatures'
    if MIN_MAX == 'max':
        fieldname = '90th percentile of maximum temperatures'
    
   
    perc_cube_allmonths = iris.cube.CubeList([])

    for month in range(0,12):
        cubelist = iris.cube.CubeList([])
        daystart = month * 30
        dayend = daystart + 30

        for day in range(daystart,dayend):
            filename = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                        CNTLNAME + '/diag_10_13/daily_percentiles/' + 
                        MIN_MAX + '_temperature_' + np.str(day) + '.nc')
            perc_cube = iris.load_cube(filename, fieldname)
            # stuff to make sure they can concatenate
            perc_cube.cell_methods = None
            perc_cube.remove_coord('t')
            perc_cube.attributes = None
            perc_cube.var_name = 'pcent_extreme'
            time = DimCoord(day,standard_name = 'time', long_name = 't', 
                              var_name = None, units = 'day')
            perc_cube.add_aux_coord(time)
            perc_cube_3d = iris.util.new_axis(perc_cube, time)

            cubelist.append(perc_cube_3d)

    
        percentile_month_cube = cubelist.concatenate_cube()
        perc_cube_allmonths.append(percentile_month_cube)

    return perc_cube_allmonths

       
    
def calculate_extreme_spell():
    """
    this function will calculate the number of days which belong to an
    extreme spell for each gridbox in the 100 year period
    """

    # read in percentiles file.  we need one of these per month
    # the cubelist will contain one cube per month
    perc_cubelist = percentiles_from_file()
    data_cube = perc_cubelist[0][0,:,:] # shape reqd for putting warm spell data
    warm = {'min' : 'cold', 'max' : 'warm'}
    warm_cold = {'min' : 'coldspell', 'max' : 'warmspell'}
   
    # get_temperature data and find where it is extreme
    if MIN_MAX == 'min':
        # this is minimum temperature
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_1'))
    if MIN_MAX == 'max':
        # this is maximum temperature
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp'))
    if MIN_MAX == 'mean':
         # this is maximum temperature
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_2'))

    # find where the data is extreme
    extreme_data = np.zeros((NYEARS * 12,30, 73, 96))
    filestart = ('/nfs/hera1/earjcti/um/' + EXPTNAME + '/pb/' + EXPTNAME
                 + 'a@pb' + EXTRA )
    
    for year in range(0,0 + NYEARS):
        print(year)
        cubelist_ann = iris.cube.CubeList([])
        for month in range(0,12):
            monix = (year * 12) + month
            percentile_cube = perc_cubelist[month]
            filename = (filestart + np.str(year).zfill(2) 
                        + MONTHNAMES.get(month) + '.nc')
            if exists(filename):
                print(filename)
            else:
                print('replacing',filename)
                filename = (filestart + np.str(year-1).zfill(2) 
                        + MONTHNAMES.get(month) + '.nc')
                print('with',filename)

            
            cubes = iris.load(filename, constraints=variable_constraint)
            cube = iris.util.squeeze(cubes[0])
            print(np.shape(cube.data))
            print(np.shape(percentile_cube.data))
         
            if MIN_MAX == 'min':
                extreme_data[monix,:,:] = np.where(cube.data 
                                                   <= percentile_cube.data, 
                                                   1.0, 0.0)
            if MIN_MAX == 'max':
                extreme_data[monix,:,:] = np.where(cube.data 
                                                   >= percentile_cube.data, 
                                                   1.0, 0.0)
         
    #use the extreme data array to find out how many wam/cold spells there
    #are 
    spell_array = np.zeros((73,96))
    for lat in range(0,73):
        print('lat is',lat)
        for lon in range(0,96):
            poss_warm_spell = 0.0
            for day in range(0,360*NYEARS):
                monix=np.int(np.floor(day/30))
                dayix = day - monix * 30
                # check extremes
                if extreme_data[monix,dayix,lat,lon] == 1.0:
                    poss_warm_spell = poss_warm_spell+1.0
                # if not an extreme reset start of warm spell
                if (extreme_data[monix,dayix,lat,lon]==0.0 
                    and poss_warm_spell > 0.0):
                    if poss_warm_spell >=6: 
                        spell_array[lat,lon] = (spell_array[lat,lon] + 
                                                poss_warm_spell)
                    poss_warm_spell=0.0
          
    spell_array = (spell_array *100.) / (360. * NYEARS)
    spell_cube = data_cube.copy(data=spell_array)
    spell_cube.units=None
    spell_cube.long_name = 'percentage of ' + warm.get(MIN_MAX) + ' spell days'
    outfile = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                   EXPTNAME + '/diag14_15/' + EXTRA + '_' + CNTLNAME + '_' + 
                   warm_cold.get(MIN_MAX) + 
                   '_diag14-15.nc') 
    iris.save(spell_cube, outfile, netcdf_format="NETCDF3_CLASSIC")
    qplt.contourf(spell_cube)
    plt.show()
#########################################################
def plot_nextreme_spell(ocn_mask):
    """
    plots the number of days in a warm spll or a cold spell
    """
    warm_cold = {'min' : 'coldspell', 'max' : 'warmspell'}
    period = {'tenvj':'mPWP', 'xozza':'PI', 'xozzb':'mPWP'}
   
    gtlt = {10 : ' < ', 5 : ' < ' ,90 : ' > ', 95  : ' > '}
    filename = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                   EXPTNAME + '/diag14_15/' + EXTRA + '_' + CNTLNAME + '_' + 
                   warm_cold.get(MIN_MAX) + 
                   '_diag14-15.nc')

    cube = iris.load_cube(filename)
    print(cube)
 
    if ocn_mask == 'y':
        if EXPTNAME == 'tenvj':
            filelsm = '/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrparm.mask.nc'
            exptlsmcube=iris.load_cube(filelsm,'LAND MASK (LOGICAL: LAND=TRUE)')
            cube.data.mask = (exptlsmcube.data - 1.0) * (-1.0)
        if EXPTNAME == 'xozza':
            filelsm = '/nfs/hera1/earjcti/ancil/preind2/qrparm.mask.nc'
            exptlsmcube=iris.load_cube(filelsm,'LAND MASK (LOGICAL: LAND=TRUE)')
            cube.data.mask = (exptlsmcube.data - 1.0) * (-1.0)
    
   
  
    if MIN_MAX == 'max' and EXPTNAME == 'tenvj' and CNTLNAME == 'xozza':
        valmin=0.
        valmax=55.
        valdiff=5.
    else:
        valmin=0.
        valmax=6.
        valdiff=1.


    cube.units = '%'
   # qplt.pcolormesh(cube_ann_avg,vmin=minval, vmax=maxval)
    qplt.contourf(cube,levels=np.arange(valmin,valmax,valdiff),
                  extend='max')
    plt.gca().coastlines()
    plt.title(('percentage of ' + period.get(EXPTNAME) + ' days that are ' 
              + warm_cold.get(MIN_MAX) + ' days in '+ period.get(CNTLNAME)),
              fontsize=10)
    filestart = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/' +
                 'diag14_15/' + EXTRA + '_' + EXPTNAME + '_' + CNTLNAME + 
                 '_percentage_' + 
                   warm_cold.get(MIN_MAX))
    plt.savefig(filestart + '.eps')
    plt.savefig(filestart + '.png')
   
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'
MIN_MAX = 'max' # max is warm spell duration index
                # min is cold spell duration index

NYEARS = 100
MONTHNAMES = {0:'ja',1:'fb',2:'mr',3:'ar',4:'my',5:'jn',6:'jl',
                  7:'ag',8:'sp',9:'ot',10:'nv',11:'dc'}
  
EXPTNAME = 'xozzb'
CNTLNAME = 'xozza'
EXTRA='o'



nextreme_spell = calculate_extreme_spell()

plot_nextreme_spell('y')
