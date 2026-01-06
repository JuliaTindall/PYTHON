#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on July 21 2020

#
# This program will produce a regridded 1X1degree timeseries of a given field.  
# We will remove the annual cycle in order to look for interannual variability 
# etc.  


import numpy as np
from netCDF4 import Dataset
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import iris.analysis.cartography
import iris.coord_categorisation
import cf_units as unit
from iris.experimental.equalise_cubes import equalise_attributes
import sys
#import os

###################################################
# get all data from files as a single cube
##################################################

def get_hadcm3_cube():
    """
    get's the datacube from hadcm3 
    """
    allcubes = iris.cube.CubeList([])
    monthnames = ['ja', 'fb', 'mr', 'ar','my','jn',
                  'jl','ag','sp','ot','nv','dc']
    
    extraalt = {'t' : 'u',
                'n' : 'o' }
    
    for year in range(STARTYEAR, ENDYEAR):
        if year < 100:
            yearuse = str(year).zfill(2)
            filestart2 = FILESTART + EXTRA + yearuse
        elif year < 200:
            yearuse = str(year -100).zfill(2)
            filestart2 = FILESTART + extraalt.get(EXTRA) + yearuse
        for month in monthnames:
            filename = filestart2 + month + '.nc'
            cube = iris.load_cube(filename, depth = 5.0)
            allcubes.append(cube)

#         if model == 'HadCM3':
 #           cubetemp.coord('t').rename('time')
 #       cubetemp.coord('time').points = (np.arange(0, 12)+((i-startyear)*12))*30.

 #       cubetemp.coord('time').units = u

 #       allcubes.append(cubetemp)


    equalise_attributes(allcubes)
    cube_temp = allcubes.concatenate_cube()
    print(cube_temp)

 #   if model == 'MRI2.3':
 #       cube_temp.coord('pressure level').rename('surface')

 #   if model == 'HadCM3' and fielduse == 'SST':
 #       cube_temp.coord('unspecified').rename('surface')

 #   if model == 'HadCM3' and fielduse == 'NearSurfaceTemperature':
 #       cube_temp.coord('ht').rename('surface')


  #  cube_temp.coord('surface').points = 0.
  #  cube = cube_temp.extract(iris.Constraint(surface=0.))

    return cube


######################################################
def cube_avg(cube):
    """
    Extract monthly averaged data from a cube

    Parameters:
    cube (iris cube): A cube with montly data that we average

    Returns:
    meanmonthcube (iris cube): the input cube averaged over each of the 12 calendar months

    """

    meanmonthcube = cube.aggregated_by('month', iris.analysis.MEAN)

    # NORESM and CESM1.2 does not start at month = 1,
    #it starts at month = 2. but should be 1
    # we are doing dome roundabout way of reordering the data
    #if modelname in ('NorESM1-F', 'NorESM-L', 'CESM1.2', 'CCSM4', 'CESM2'):
    if modelname in ('NorESM1-F', 'NorESM-L'):
        allcubes = iris.cube.CubeList([])
        for mon in range(2, 13):
            slice = meanmonthcube.extract(iris.Constraint(month=mon))
            # attempt to reorder time coordinate
            slice.coord('time').points = mon-1
            slice.coord('month').points = mon-1
            slice.coord('time').bounds = None
            allcubes.append(slice)
        # do december (month 1)
        slice = meanmonthcube.extract(iris.Constraint(month=1))
        slice.coord('time').points = 12
        slice.coord('month').points = 12
        slice.coord('time').bounds = None
        allcubes.append(slice)
        #process
        meanmonthcube = (allcubes.merge_cube())

    meanmonthcube.long_name = fieldnameout

    iris.util.promote_aux_coord_to_dim_coord(meanmonthcube, 'month')

    return meanmonthcube


def remove_ann_cycle(cube, mon_avg_cube):
    """
    removes the annual cycle stored in mon_avg_cube from cube
    """

    cubedata = cube.data
    for i, monthno in enumerate(cube.coord('month').points):
        mon_avg_data = (mon_avg_cube.extract(iris.Constraint(month=monthno))).data
        cubedata[i, :, :] = cubedata[i, :, :] - mon_avg_data


    timeseries_cube = cube.copy(data=cubedata)

    return timeseries_cube


##############################################
def main():
    """
    get average seasonal cycle for the cubes
    """

    cube = get_hadcm3_cube()
    sys.exit(0)

    # now regrid the cube onto a 1X1 grid (we will first try regridding the raw data)
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'

    # do not need to regrid CCSM4_UoTdata or a field that was originally on a tripolar grid

    if ((modelname == 'CCSM4-UoT')
            or (modelname == 'IPSLCM5A' and fieldnamein == 'tos')
            or (modelname == 'IPSLCM5A2' and fieldnamein == 'tos')):
        regridded_cube = cube
    else:
        cubegrid = iris.load_cube('one_lev_one_deg.nc')
        regridded_cube = cube.regrid(cubegrid, iris.analysis.Linear())


    refdate = 'days since 0800-01-01 00:00:00'

    # for cosmos
    if modelname == 'COSMOS':
        # cosmos data is in a strange time coordinate line yyyymmdd
        # we need to convert it to days since reference time
        origpoints = regridded_cube.coord('time').points
        npoints = len(origpoints)
        yeararr = np.zeros(npoints)
        montharr = np.zeros(npoints)
        dayarr = np.zeros(npoints)
        daydecimal = np.zeros(npoints)
        dayssinceref = np.zeros(npoints)
        for i in range(0, npoints):
            origstr = str(origpoints[i])
            yeararr[i] = origstr[:][0:4]
            montharr[i] = origstr[:][4:6]
            dayarr[i] = origstr[:][6:8]
            daydecimal[i] = origstr[:][8:]
            dayssinceref[i] = dayssinceref[i-1]+dayarr[i]+daydecimal[i]-daydecimal[i-1]
        # subtract 1 from days since reference date (as reference date will be 1st Jan)
        dayssinceref = dayssinceref-1


        regridded_cube.coord('time').points = dayssinceref
        #  end of COSMOS loop

    # for EC-Earth3.1
    if modelname == 'EC-Earth3.1':
    # convert from hours to days
        origpoints = regridded_cube.coord('time').points
        newpoints = origpoints/24.
        regridded_cube.coord('time').points = newpoints
        refdate = 'days since 2390-01-01 00:00:00'


    # regrid to mm/day from kg/m2/s if required
    if modelname in ('EC-Earth3.1', 'EC-Earth3.3', 'IPSLCM5A',
                     'IPSLCM5A2', 'IPSLCM6A', 'CCSM4-Utr', 'GISS2.1G'):
        if fieldnamein == 'pr':
            regridded_cube.data = regridded_cube.data * 60. *60. *24.
            cube.data = cube.data* 60. *60. *24.
            regridded_cube.name = 'Total precipitation'
            regridded_cube.long_name = 'Total precipitation'
            regridded_cube.units = 'mm/day'



    if modelname in ('NorESM1-F', 'NorESM-L', 'CESM1.2', 'CESM2', 'CCSM4'):

       # if precipitation is in m/s convert to mm/day
        if fieldnamein == 'pr':
            regridded_cube.data = regridded_cube.data * 60. * 60. * 24. * 1000.
            cube.data = cube.data * 60. * 60. * 24. * 1000.
            regridded_cube.name = 'Total precipitation'
            regridded_cube.long_name = 'Total precipitation'
            regridded_cube.units = 'mm/day'

    if modelname in ('CCSM4-UoT', 'NorESM1-F', 'NorESM-L', 'IPSLCM6A',
                     'EC-Earth3.1', 'EC-Earth3.3', 'IPSLCM5A', 'IPSLCM5A2',
                     'HadCM3', 'GISS2.1G'):
         # convert to celcius
        if fieldnamein in ('tas', 'tos'):
            regridded_cube.convert_units('Celsius')
            cube.convert_units('Celsius')


    if modelname in ('COSMOS', 'MIROC4m', 'IPSLCM6A', 'EC-Earth3.1'):
        regridded_cube.coord('time').units = refdate


    # add auxillary coordinates month and year
    iris.coord_categorisation.add_month_number(regridded_cube, 'time', name='month')
    iris.coord_categorisation.add_year(regridded_cube, 'time', name='year')


     # correct the start month if required
    regridded_cube = correct_start_month(regridded_cube)
    
    # calculate averages
    mean_mon_cube = cube_avg(regridded_cube)



    # remove annual cyble
    anom_cube = remove_ann_cycle(regridded_cube, mean_mon_cube)

    # to check we have removed the average properly get the monthly
    # average of the anomaly cube it should be zero

    new_mean_mon_cube = cube_avg(anom_cube)
    qplt.contourf(new_mean_mon_cube[2, :, :], levels=np.arange(-0.01, 0.011, 0.001), extend='both')
    plt.show()



    # write the cubes out to a file

    outfile = outstart+'timeseries_no_ann_cycle.nc'
    iris.save(anom_cube, outfile, netcdf_format='NETCDF3_CLASSIC', fill_value=2.0E20)






##########################################################
# main program

exptname = {
        "e280" : "tenvo",
        "e400" : "tenvq",
        "e560":"tenvs",
        "eoi400" : "tenvj",
        "eoi350" : "tenvk",
        "eoi450" : "tenvl",
        "eoi280" : "tenvm"
      
}

extraname = {
        "e280" : "t",
        "e400" : "t",
        "eoi400" : "n",
        "eoi350" : "o",
        "eoi450" : "o",
        "eoi280" : "o",
        "e560" : "v"}

LINUX_WIN='l'
STARTYEAR=0
ENDYEAR=200
TIMEPERIOD='e280'
expt=exptname.get(TIMEPERIOD)
EXTRA=extraname.get(TIMEPERIOD)

FILESTART = '/nfs/hera1/earjcti/um/' + expt + '/' + expt + 'o@pf'
main()
