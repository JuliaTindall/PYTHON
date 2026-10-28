#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 3rd 2019
#modified October 2026


#
# This program will regrid some of the data that is needed for PLIOMIP2.
# We will put 100 year average fields onto a 1deg X 1deg standard grid
# it can be used where experiments have been uploaded with 100 years in
# one file

# ammended to do PlioMIP3 fields



import numpy as np
import iris
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import iris.analysis.cartography
import iris.coord_categorisation
from iris.cube import CubeList
import cf_units as unit
import cftime
import sys
#import os

###################################################
# get all data from files as a single cube
##################################################



def get_miroc_tos():
    """
    get miroc data
    """
    cube = iris.load_cube(filename)

    # MIROC didn't have the units set up correctly for 'tos'
    cube.coord('latitude').units = 'degrees'
    cube.coord('longitude').units = 'degrees'
    return(cube)


    
def get_ccsm4_2deg():
    """
    get ccsm4_2deg utrecht data
    """
    
    allcube = iris.load(filename)
    ncubes = len(allcube)
    for i in range(0, ncubes):
        print(allcube[i].var_name,fielduse)
        if allcube[i].var_name == fielduse:
            cube = allcube[i]
    
    # put units as celcius if required
    if fielduse == 'tas':
        print(cube.units)
        cube.units='Celsius'
    
    cube2 = iris.util.new_axis(cube, 'time')
    print('julia check')
    
    return cube


def get_cesm13_singlecube(filename_, fielduse_,expt_):
    """
    get a single cube from a single file
    this is needed because total precipitation is the sum of two cubes
    """
    fileends = []
    if expt_ == 'LP':
        fileends.append('.089001-093912_regrid.nc')
        fileends.append('.094001-098912_regrid.nc')
    if expt_ == 'PI':
        fileends.append('.050101-055012_regrid.nc')
        fileends.append('.055101-060012_regrid.nc')

    cubes = CubeList([])
    for end in fileends:
        fileuse = filename_ + end
        allcube = iris.load(fileuse)
        ncubes = len(allcube)
        for i in range(0, ncubes):
            #print(allcube[i].var_name,fielduse_)
            if allcube[i].var_name == fielduse_:
                filecube = allcube[i]
                cubes.append(filecube)

    iris.util.equalise_attributes(cubes)
    singlecube = cubes.concatenate_cube()
            
    return singlecube
    
def get_cesm13(expt):
    """
    get cesm1.3 data
    """
    
    if type(fielduse) is list: # perhaps a list containing large scale and convective precip
                               # that need to be added together for total precipitaiton
       cube1 = get_cesm13_singlecube(filename[0], fielduse[0],expt)
       cube2 = get_cesm13_singlecube(filename[1], fielduse[1],expt)

       
       if fieldnameout == 'TotalPrecipitation':
           cube = cube1 + cube2
       if fieldnameout == 'SST':
           # mask temperature (cube1) by ice fraction (cube2)
           cube1.convert_units('Celsius')
           cubedata = np.where((cube2.data > 0.01), 
                               -1.8, cube1.data)
           cube = cube1.copy(data=cubedata)
    else:
       cube = get_cesm12_singlecube(filename, fielduse)
    
    # put units as celcius if required
    if fielduse == 'TREFHT':
        print(cube.units)
        cube.convert_units('Celsius')

    if fieldnameout == 'SST':
        # we need to mask out the land
        if exptnamein == 'Eoi400':
            filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' + 
                       'b40.B1850.f09_g16.PMIP4-pliomip2.' + 
                       'cam.h0.LANDFRAC.0851.0950.nc')
        if exptnamein == 'E280':
            filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' + 
                       'b.e12.B1850.f09_g16.preind.' + 
                       'cam.h0.LANDFRAC.0701.0800.nc')
        lsmcube = get_cesm12_singlecube(filelsm, 'LANDFRAC')
   
        nt, nx, ny = np.shape(cube.data)
        mymask_2d = lsmcube[0].data
        mymask_2d = np.where(mymask_2d > 0.01,  1.0, 0.0)
        mymask = np.vstack([mymask_2d] * nt)

        cube.data  =  np.ma.array(cube.data,  mask = mymask)
    
    
    return cube

######################################################
def correct_start_month(cube):
    """
    parameters: cube
    returns: the same cube
    
    if month doesn't start on january we will have to change some of the
    years from the end to match those at the start to give 100 full years

    
    """
    print('julia3')

    print(cube.coord('time').points)
    print(cube.coord('time').units)
    print(cube.coord('year').points)
    print(cube.coord('month').points)

    # also check there are 1200 points
    if len(cube.coord('time').points) != 1200.:
           print('wrong number of points',len(cube.coord('time').points))
           sys.exit(0)
   
    startyear = (cube.coord('year').points[0])
    endyear = (cube.coord('year').points[-1])
    # count the number of months that have the same year as the first index
    nstart = 0
    nend = 0
    for i in range(0, 12):
        if cube.coord('year').points[i]  == startyear:
            nstart = nstart+1
    for i in range(-13, 0):
        if cube.coord('year').points[i]  == endyear:
            nend = nend+1
    if nend!= 12 or nstart!= 12:
        # oops we don't have a full year at the start and the end
        # check to see whether we can aggregate into one year
        if nend+nstart  == 12:
            for i in range(-1 * nend, 0):
                # move the last few to the start of the cube
                cube.coord('year').points[i] = startyear


        else:

            print('you have a partial year somewhere')
            print('correct input data to provide full years')
            print(nend, nstart)
            sys.exit(0)

    return cube

######################################################    
def cube_avg(cube):
    """
    Extract averaged data and standard deviation values from a cube

    Parameters:
    cube (iris cube): A cube with montly data that we average

    Returns:
    meanmonthcube (iris cube): the input cube averaged over each of the 12 calendar months
    sdmonthcube (iris cube): the standard deviation on the values in meanmonthcube
    meanyearcube (iris cube): annual average values for each calendar year
    meancube (iris cube) : the long term mean value at each gridpoint
    sdcube (iris cube): the interannual standard deviation on the mean at each gridpoint
    """

    meanmonthcube = cube.aggregated_by('month', iris.analysis.MEAN)
    
    # NORESM and CESM1.2 does not start at month = 1,
    #it starts at month = 2. but should be 1
    # we are doing dome roundabout way of reordering the data
    if (modelname  == 'NorESM1-F' or modelname  == 'NorESM-L'
       or modelname == 'iCESM1.3' or modelname == 'CCSM4-1deg'):
        allcubes = iris.cube.CubeList([])
        for mon in range(2, 13):
            slice  =  meanmonthcube.extract(iris.Constraint(month = mon))
            # attempt to reorder time coordinate
            slice.coord('time').points = mon-1
            slice.coord('month').points = mon-1
            slice.coord('time').bounds = None
            allcubes.append(slice)
        # do december (month 1)
        slice  =  meanmonthcube.extract(iris.Constraint(month = 1))
        slice.coord('time').points = 12
        slice.coord('month').points = 12
        slice.coord('time').bounds = None
        allcubes.append(slice)
        #process
        meanmonthcube = (allcubes.merge_cube())
        print(meanmonthcube.coord('time').points)
        print(meanmonthcube.coord('month').points)

    meanmonthcube.long_name = fieldnameout

    print('months',meanmonthcube.coord('month').points)
    iris.util.promote_aux_coord_to_dim_coord(meanmonthcube, 'month')
    sdmonthcube = cube.aggregated_by('month', iris.analysis.STD_DEV)
    sdmonthcube.long_name = fieldnameout

    # mean and standard deviation for each year
    meanyearcube = cube.aggregated_by('year', iris.analysis.MEAN)
    meanyearcube.long_name=fieldnameout

    # mean and interannual standard deviation over dataset
    meancube = meanmonthcube.collapsed('time', iris.analysis.MEAN)
    meancube.long_name=fieldnameout
    sdcube = meanyearcube.collapsed('time', iris.analysis.STD_DEV)
    sdcube.long_name=fieldnameout

    return [meanmonthcube, sdmonthcube, meanyearcube, sdcube, meancube]


def mon_avg(cube):
    """
    get the monthly average data from a cube

    Parameters:
    cube (iris cube) contains monthly data (ie t=12)

    Returns:
    meanmon (np array) with mean global values for each of the 12 months
    """

    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    grid_areas2 = iris.analysis.cartography.area_weights(cube)
    tempcube = cube.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas2)
    meanmon = tempcube.data

    return meanmon

def get_monthly_sd(cube):
    """
    get monthly values of the standard deviation on the
    globally averaged mean from a cube

    Parameters:
        cube (iris cube) : a cube with monthly data from the whole dataset
                           (tdim=nyears*nmonths)
    Returns:
        stdevmon (numpy array) : an array with the monthly standard deviation
                                 on the monthly global mean

    """


    stdevmon = np.zeros(12)
    for mon in range(1, 13):
        # mon_slice is the slice from this month.  It should have a time
        # dimension of t where t is the number of years
        mon_slice = cube.extract(iris.Constraint(month=mon))
        mon_slice.coord('latitude').guess_bounds()
        mon_slice.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(mon_slice)
        mon_avg = mon_slice.collapsed(['latitude', 'longitude'],
                                      iris.analysis.MEAN, weights=grid_areas)
        stdevcube = mon_avg.collapsed(['time'], iris.analysis.STD_DEV)
        stdevmon[mon-1] = stdevcube.data

    return stdevmon



def area_means(allmeancube, yearmeancube):
    """
    get the globally averaged mean and standard deviation
    and latitudinal means and standard deviation

    Parameters:
    allmeancube (iris cube) lat/lon cube of temporal mean tdim=1
    yearmeancube (iris cube)  lat/lon/time cube of yearly averages tdim=nyears

    Returns
    meanann (float) global long term mean
    stdevann (float)  global long term standard deviation of the mean
    mean lat (np.array(nlats) latitudinally averaged long term mean
    stdevlat (np.array(nlats)standard deviation of the latitudinal mean

    """

    allmeancube.coord('latitude').guess_bounds()
    allmeancube.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(allmeancube)
    tempcube = mean_data.collapsed(['latitude', 'longitude'],
                                   iris.analysis.MEAN, weights=grid_areas)
    meanann = tempcube.data


    # get mean for each latitude
    tempcube = allmeancube.collapsed(['longitude'], iris.analysis.MEAN)
    meanlat = tempcube.data
    meanlat = np.squeeze(meanlat)


    # get standard deviation
    # 1. mean for each year

    yearmeancube.coord('latitude').guess_bounds()
    yearmeancube.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(yearmeancube)
    tempcube = yearmeancube.collapsed(['latitude', 'longitude'],
                                      iris.analysis.MEAN, weights=grid_areas)
    stdevcube = tempcube.collapsed(['time'], iris.analysis.STD_DEV)
    stdevann = stdevcube.data


    # get standard deviation for each latitude
    tempcube = yearmeancube.collapsed(['longitude'], iris.analysis.MEAN)
    stdevcube = tempcube.collapsed(['time'], iris.analysis.STD_DEV)
    #stdevann=stdevcube.data
    stdevlat = stdevcube.data
    stdevlat = np.squeeze(stdevlat)

    return [meanann, stdevann, meanlat, stdevlat]

##############################################
def regrid_data(fieldnamein, exptnamein):
    """
    regrid the data
    """


    print('moodelname is', modelname)
    print('filename is', filename)
    print('fielduse is', fielduse)
    print('exptnamein is', exptnamein)
    print('exptnameout is',exptnameout)
    
    # outfile
    if linux_win  == 'l':
        outstart = ('/uolstore/Research/a/hera1/earjcti/regridded100_P3/'+modelname+'/'+exptnameout+'.'+
        fieldnameout+'.')
        lsmstart = '/uolstore/Research/a/hera1/earjcti/regridded100_P3/'
    else:
        outstart = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
              +modelname+'\\'+exptnameout+'.'+fieldnameout+'.')
        lsmstart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'


   
    #####################################
    # get all data in a single cube


    if modelname == 'iCESM1.3':
        cube = get_cesm13(exptnamein)
    else:
        cube = iris.load_cube(filename)

    # if too many points extract last 1200 times (assuming time is last dim)
    cube = cube[-1200:]
   
    # sort out inconsistencies
    if modelname == 'ACCESS':
        cube.coord('latitude').coord_system = None
        cube.coord('longitude').coord_system = None
    if modelname == 'HadCM3':
        cube.coord('t').rename('time')
        cube =iris.util.squeeze(cube)
        
    ndim = cube.ndim
   


    # now regrid the cube onto a 1X1 grid (we will first try regridding the raw data)
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'

    # do not need to regrid UofTdata or a field that was originally on a tripolar grid
    print('julia1')
    if ((modelname   == 'UofT-CCSM4')
        or (modelname  == 'IPSLCM5A' and fieldnamein  == 'tos')
        or (modelname  == 'IPSLCM5A2' and fieldnamein  == 'tos')):
        regridded_cube = cube
    else:
        cubegrid = iris.load_cube('one_lev_one_deg.nc')
        regridded_cube = cube.regrid(cubegrid, iris.analysis.Linear())


  
    #refdate = 'days since 0800-01-01 00:00:00'

    # awi models date is in a strange time coordinate line yyyymmdd
    # we need to convert it to days since reference time
  
    if (modelname  == 'COSMOS' or modelname == 'AWI-CM2'
        or modelname =='AWI-ESM2') and (fieldnamein !='tos'):
        # cosmos data is in a strange time coordinate line yyyymmdd
        # we need to convert it to days since reference time
        origpoints = regridded_cube.coord('time').points
        dates = []
        for point in origpoints:
            # in case point is something like 18500115.0
            datestr = str(int(point))

            year = int(datestr[0:4])
            month = int(datestr[4:6])
            day = int(datestr[6:8])

            dates.append(cftime.DatetimeGregorian(year, month, day))

        new_units = unit.Unit('days since 0001-01-01 00:00:00',
                                  calendar='gregorian')
            
        newpoints = new_units.date2num(dates)

        regridded_cube.coord('time').points = newpoints
        regridded_cube.coord('time').units = new_units
        


    # for EC-Earth3.1
    if modelname  == 'EC-Earth3.1':
    # convert from hours to days
        origpoints = regridded_cube.coord('time').points
        newpoints = origpoints/24.
        regridded_cube.coord('time').points = newpoints
        refdate = 'days since 2390-01-01 00:00:00'


    # regrid to mm/day from kg/m2/s if required
    if (modelname  == 'EC-Earth3.1' or modelname  == 'IPSLCM5A'
             or modelname  == 'IPSLCM5A2' or modelname  == 'IPSLCM6A'
             or modelname == 'CCSM4-2deg' or modelname =='GISS'
             or modelname == 'ACCESS'):
        if fieldnamein  == 'pr':
            regridded_cube.data = regridded_cube.data * 60. *60. *24.
            cube.data = cube.data* 60. *60. *24.
            regridded_cube.name = 'Total precipitation'
            regridded_cube.long_name = 'Total precipitation'
            regridded_cube.units = 'mm/day'



    if (modelname  == 'UofT' or modelname  == 'NorESM1-F' 
        or modelname  == 'NorESM-L' or modelname == 'CESM1.2'
        or modelname == 'CCSM4-1deg'):
       # if precipitation is in m/s convert to mm/day
        if fieldnamein  == 'pr':
            regridded_cube.data = regridded_cube.data * 60. *60. *24. *1000.
            cube.data = cube.data* 60. *60. *24. *1000.
            regridded_cube.name = 'Total precipitation'
            regridded_cube.long_name = 'Total precipitation'
            regridded_cube.units = 'mm/day'

    if (modelname  == 'UofT' or modelname  == 'NorESM1-F' or modelname  == 'NorESM-L'
        or modelname  == 'IPSLCM6A' or modelname  == 'EC-Earth3.1'
        or modelname  == 'IPSLCM5A' or modelname  == 'IPSLCM5A2'
        or modelname  == 'HadCM3' or modelname == 'GISS'
        or modelname == 'ACCESS'):
         # convert to celcius
        if fieldnamein  == 'tas' or fieldnamein  == 'tos':
            regridded_cube.convert_units('Celsius')
            cube.convert_units('Celsius')


    if modelname  == 'UofT':
        # we need to add the missing time coordinate
        points = (np.arange(0, 1200)*30)+15. # go for middle of month
        u  =  unit.Unit('days since 0800-01-01 00:00:00',
               calendar = unit.CALENDAR_360_DAY) # put as 360 day calendar because of the way
                                               # the data was sent.

        regridded_cube.add_dim_coord(iris.coords.DimCoord(points,
                standard_name = 'time',  long_name = 'time',
                var_name = 'time',
                units = u,
                bounds = None,
                coord_system = None,  circular = False), 0)
    elif (modelname  == 'IPSLCM6A' or modelname  == 'EC-Earth3.1'):
        regridded_cube.coord('time').units = refdate



       # end of Uof T loop

    print('julia2')
    #print(regridded_cube.coord('time'))

    # add auxillary coordinates month and year
    iris.coord_categorisation.add_month_number(regridded_cube,  'time',  name = 'month')
    iris.coord_categorisation.add_year(regridded_cube,  'time',  name = 'year')

    print(regridded_cube.coord('month').points[0:100])
    # correct the start month if required
    regridded_cube=correct_start_month(regridded_cube)
    
    # calculate averages
    mean_mon_data, sd_mon_data, mean_year_data, sd_data, mean_data = cube_avg(regridded_cube)


   
    print('j4 all mon', np.shape(mean_data))
    print('j5 all year', np.shape(sd_data))

    # extract monthly data from cubes

    jan_slice  =  regridded_cube.extract(iris.Constraint(month = 1))
    feb_slice  =  regridded_cube.extract(iris.Constraint(month = 2))
    mar_slice  =  regridded_cube.extract(iris.Constraint(month = 3))
    apr_slice  =  regridded_cube.extract(iris.Constraint(month = 4))
    may_slice  =  regridded_cube.extract(iris.Constraint(month = 5))
    jun_slice  =  regridded_cube.extract(iris.Constraint(month = 6))
    jul_slice  =  regridded_cube.extract(iris.Constraint(month = 7))
    aug_slice  =  regridded_cube.extract(iris.Constraint(month = 8))
    sep_slice  =  regridded_cube.extract(iris.Constraint(month = 9))
    oct_slice  =  regridded_cube.extract(iris.Constraint(month = 10))
    nov_slice  =  regridded_cube.extract(iris.Constraint(month = 11))
    dec_slice  =  regridded_cube.extract(iris.Constraint(month = 12))




    # write the cubes out to a file

    outfile = outstart+'mean_month.nc'
    iris.save(mean_mon_data, outfile, netcdf_format = 'NETCDF3_CLASSIC', fill_value = 2.0E20)

    outfile = outstart+'sd_month.nc'
    iris.save(sd_mon_data, outfile, netcdf_format = 'NETCDF3_CLASSIC', fill_value = 2.0E20)

    print('mean cube',mean_data)
    
    outfile = outstart+'allmean.nc'
    iris.save(mean_data, outfile, netcdf_format = 'NETCDF3_CLASSIC', fill_value = 2.0E20)

    outfile = outstart+'allstdev.nc'
    iris.save(sd_data, outfile, netcdf_format = 'NETCDF3_CLASSIC', fill_value = 2.0E20)

    ##########################################################################
    # get the global mean and standard deviation and write them all out to a file
    #
    textout = outstart+'data.txt'

    file1 =  open(textout, "w")

    # get mean field for cube

    mean_data.coord('latitude').guess_bounds()
    mean_data.coord('longitude').guess_bounds()
    grid_areas  =  iris.analysis.cartography.area_weights(mean_data)
    tempcube = mean_data.collapsed(['latitude', 'longitude'],
                                iris.analysis.MEAN, weights = grid_areas)
    meanann = tempcube.data


    # get mean for each latitude
    tempcube = mean_data.collapsed(['longitude'], iris.analysis.MEAN)
    meanlat = tempcube.data
    meanlat = np.squeeze(meanlat)


    # get standard deviation
    # 1. mean for each year

    mean_year_data.coord('latitude').guess_bounds()
    mean_year_data.coord('longitude').guess_bounds()
    grid_areas  =  iris.analysis.cartography.area_weights(mean_year_data)
    tempcube = mean_year_data.collapsed(['latitude', 'longitude'],
                                iris.analysis.MEAN, weights = grid_areas)
    stdevcube = tempcube.collapsed(['time'], iris.analysis.STD_DEV)
    stdevann = stdevcube.data

    #plt.plot(tempcube.data)
    #plt.plot([0, 100], [meanann, meanann])
    #plt.plot([0, 100], [meanann+stdevann+stdevann, meanann+stdevann+stdevann])
    #plt.plot([0, 100], [meanann-stdevann-stdevann, meanann-stdevann-stdevann])
    #plt.title('data and data+/-2sd')
    #plt.show()
    #plt.close()

    # get standard deviation for each latitude
    tempcube = mean_year_data.collapsed(['longitude'],  iris.analysis.MEAN)
    stdevcube = tempcube.collapsed(['time'], iris.analysis.STD_DEV)
    #stdevann = stdevcube.data
    stdevlat = stdevcube.data
    stdevlat = np.squeeze(stdevlat)






    # write out to a file
    file1.write('global annual mean and standard deviation\n')
    file1.write('------------------------------------------\n')
    if ndim>= 4:
        file1.write(str(np.round(meanann[0], 2))+', '+str(np.round(stdevann[0], 3))+'\n')
    else:
        file1.write(str(np.round(meanann, 2))+', '+str(np.round(stdevann, 3))+'\n')

    # get monthly means and standard deviation
    file1.write('monthly means and standard deviations \n')
    file1.write('----------------------------------------')
    file1.write('month    mean    sd  \n')

    mean_mon_data.coord('latitude').guess_bounds()
    mean_mon_data.coord('longitude').guess_bounds()
    grid_areas2  =  iris.analysis.cartography.area_weights(mean_mon_data)
    tempcube = mean_mon_data.collapsed(['latitude', 'longitude'],
                                iris.analysis.MEAN, weights = grid_areas2)
    meanmon = tempcube.data

    # get monthly average using grid areas from year average
    # to calculate standard deviation
    #print(np.shape(jan_slice))
    #print(np.shape(grid_areas))
    #sys.exit(0)

    jan_avg = jan_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    feb_avg = feb_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    mar_avg = mar_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    apr_avg = apr_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    may_avg = may_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    jun_avg = jun_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    jul_avg = jul_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    aug_avg = aug_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    sep_avg = sep_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    oct_avg = oct_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    nov_avg = nov_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)
    dec_avg = dec_slice.collapsed(['latitude', 'longitude'], iris.analysis.MEAN, weights = grid_areas)

    stdevmon = np.zeros(12)

    stdevcube = jan_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[0] = stdevcube.data
    stdevcube = feb_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[1] = stdevcube.data
    stdevcube = mar_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[2] = stdevcube.data
    stdevcube = apr_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[3] = stdevcube.data
    stdevcube = may_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[4] = stdevcube.data
    stdevcube = jun_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[5] = stdevcube.data
    stdevcube = jul_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[6] = stdevcube.data
    stdevcube = aug_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[7] = stdevcube.data
    stdevcube = sep_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[8] = stdevcube.data
    stdevcube = oct_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[9] = stdevcube.data
    stdevcube = nov_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[10] = stdevcube.data
    stdevcube = dec_avg.collapsed(['time'], iris.analysis.STD_DEV)
    stdevmon[11] = stdevcube.data

    for i in range(0, 12):
        if ndim>= 4:
            file1.write(str(i+1)+', '+str(np.round(meanmon[i].data[0], 2))+', '+str(np.round(stdevmon[i], 3))+'\n')
        else:
            file1.write(str(i+1)+', '+str(np.round(meanmon[i], 2))+', '+str(np.round(stdevmon[i], 3))+'\n')

    # get latitudinal means and standard deviation
    file1.write('zonal means and standard deviations \n')
    file1.write('----------------------------------------\n')
    file1.write('latitude    mean    sd  \n')
    for i in range(0, len(meanlat)):
        file1.write(str(mean_data.coord('latitude').points[i])+', '+str(np.round(meanlat[i], 2))+', '+str(np.round(stdevlat[i], 3))+'\n')

    file1.close()


    ########################################################
    # check that we have averaged properly.  To do this we are
    # going to plot the annual cycle of the global mean field
    # for the regridded dataset and also for each year

   
    #global mean
    #plt.subplot(2, 2, 1) # global mean from each year
    #subcube = subcube_mean_mon.copy(data = )  # set up structure of subcube

    #if cube.coord('latitude').has_bounds():
    #    cube.coord('latitude').bounds
    #else:
    #    cube.coord('latitude').guess_bounds()

    #if cube.coord('longitude').has_bounds():
    #    cube.coord('longitude').bounds
    #else:
    #    cube.coord('longitude').guess_bounds()

    #grid_areas  =  iris.analysis.cartography.area_weights(cube)
    #newcube = cube.collapsed(['latitude', 'longitude'],  iris.analysis.MEAN, weights = grid_areas)
    #nt = len(newcube.data)
    #nyears = int(nt/12)
    #print('nyears is', nyears)

    #for i in range(0, nyears):
    #    tstart = i*12
    #    tend = (i+1)*12
    #    plotdata = newcube.data[tstart:tend]
    #    plt.plot(plotdata, color = 'r')



    # global mean from average

    #grid_areas  =  iris.analysis.cartography.area_weights(mean_mon_data)
    #temporal_mean  =  mean_mon_data.collapsed(['latitude', 'longitude'],  iris.analysis.MEAN, weights = grid_areas)

    #plt.plot(temporal_mean.data, color = 'b', label = 'avg')
    #plt.title('globavg '+fieldnamein)
    #plt.legend()




    # check at 30N
    #plt.subplot(2, 2, 2)

    #bounds = cube.coord('latitude').bounds
    #nbounds = cube.coord('latitude').nbounds
    #nbounds, dummy = np.shape(bounds)
    #for i in range(0, nbounds):
    #    if (bounds[i, 0]>= 32. >= bounds[i, 1] or bounds[i, 0]<= 32. <bounds[i, 1]):
    #        index = i
    #print(cube)
    #if ndim>= 4:
    #    subcube = cube[:, :, index, :]
    #else:
    #    subcube = cube[:, index, :]
    #cube_avg_30N = subcube.collapsed(['longitude'], iris.analysis.MEAN)

    #for i in range(0, nyears):
    #    tstart = i*12
    #    tend = (i+1)*12
    #    plotdata = cube_avg_30N.data[tstart:tend]
    #    plt.plot(plotdata, color = 'r')

    #mean at 30N
    #slice_30N =  mean_mon_data.extract(iris.Constraint(latitude = 32))
    #mean_30N = slice_30N.collapsed(['longitude'],  iris.analysis.MEAN)


    #plt.plot(mean_30N.data, color = 'b', label = 'avg')
    #plt.title('average at 30N by month')
    #plt.legend()
    #plt.show()
    #plt.close()


#############################################################################
def getnames(modelname, filestart, fieldnamein, exptnamein):

# this program will get the names of the files and the field for each
# of the model

    # set up model specific dictionaries
    MIROC_FIELDS  = {"pr" : "Precipitation",
        "tas" : "Near-Surface Air Temperature",
        "sic" : "SeaIceAreaFraction",
        "tos" : "tos"
        }


    ECearth_FIELDS  = {"pr" : "totp",
        "tas" : "tas",
        "tos" : "sst",
        "sic" : "SeaIceAreaFraction"
        }

    IPSLCM5A_FIELDS  = {"pr" : "TotalPrecip_pr",
        "tas" : "NearSurfaceTemp_tas",
        "sic" : "SeaIceAreaFraction",
        "tos": "SeasurfaceTemp_sst"
        }

    NorESM_FIELDS = {"pr" : "PRECT",
        "tas" : "TREFHT",
        "sic" : "SeaIceAreaFraction",
        "tos" : "sst"
        }
    
    CCSM42_FIELDS = {"pr" : "TotalPrecip",
                      "tas" : "NearSurfaceAirTemp",
                      "sic" : "SeaIceAreaFraction",
                      "tos" : "SeaSurfaceTemp"
                      }
    ACCESS_FIELDS = {"tas": "TEMPERATURE AT 1.5M",
                     "tos": "Conservative temperature",
                     "pr": "TOTAL PRECIPITATION RATE     KG/M2/S"}
    HadCM3_FIELDS = {"tas": "TEMPERATURE AT 1.5M",
                     "pr": "TOTAL PRECIPITATION RATE"}
    AWI_FIELDS = {"tas": "2m temperature",
                     "pr": "large scale precipitation"}
    AWI_FILE_FIELDS = {"pr" : "totalprecip"}
    CESM12_FIELDS = {"pr" : "TotalPrecip",
                     "tas" : "TREFHT",
                     "sic" : "SeaIceAreaFraction",
                     "tos" : "TS"
                     }

    CESM12_EXTRA =  {"Eoi400": "b.e12.B1850.f09_g16.PMIP4-pliomip2.cam.h0.",
                     "E280": "b.e12.B1850.f09_g16.preind.cam.h0.",
                    }
    
    CCSM4_EXTRA =  {"Eoi400": "b40.B1850.f09_g16.PMIP4-pliomip2.cam.h0.",
                     "E280": "b40.B1850.f09_g16.preind.cam.h0.",
                    }

    ECearth_EXPT = {"Eoi400": "mPlio",
              "E280":"PI"
              }
    HadCM3_EXPT = {"LP": "xqbwd","PI": "xqbwc", "EP" : "xqbwg",
                   "LP_min_LSM_change":"xqbwr"}
    
    CESM12_EXPT = {"Eoi400": "PlioMIP2",
              "E280":"PI"
              }

    IPSLCM5A_EXPT = {"Eoi400": "Eoi400",
              "E280":"PI"
              }

    CESM12_TIME = {"E280" : ".0701.0800",
                   "Eoi400" : ".0806.0905"
                   }
    
    CCSM4_TIME = {"Eoi400" : ".0851.0950",
                  "E280" : ".0081.0180"
                  }

    IPSLCM5A_TIME = {"Eoi400": "3581_3680",
              "E280":"3600_3699"
              }

    IPSLCM5A21_TIME = {"Eoi400": "3381_3480",
              "E280":"6110_6209",
              }
    IPSLCM6A_TIME = {"Eoi400": "midPliocene-eoi400_r1i1p1f1_gr_185001-204912",
              "E280":"piControl_r1i1p1f1_gr_285001-304912",
              }
    IPSLCM6A_TIME_ALT = {"Eoi400": "midPliocene-eoi400_r1i1p1f1_gn_185001-204912",
              "E280":"piControl_r1i1p1f1_gn_285001-304912",
              }
    GISS_TIME1 = {"Eoi400": "midPliocene-eoi400_r1i1p1f1_gn_305101-310012",
              "E280":"piControl_r1i1p1f1_gn_490101-495012",
              "E560": "abrupt-2xCO2_r1i1p1f1_gn_190101-195012"
              }
    GISS_TIME2 = {"Eoi400": "midPliocene-eoi400_r1i1p1f1_gn_310101-315012",
              "E280":"piControl_r1i1p1f1_gn_495101-500012",
              "E560": "abrupt-2xCO2_r1i1p1f1_gn_195101-200012"
              }
    atm_ocn_ind = {"tas": "Amon",
                 "pr": "Amon",
                 "tos":"Omon"}
   
    # get names for each model
    if modelname == 'ACCESS':
        filename = (filestart + modelname + '/' + exptnamein + '/' +
                    fieldnamein)
        #print(fieldname)
        if fieldnamein == 'tos':
            filename = filename + '_Omon_ACCESS-ESM1-5_'
        else:
            filename = filename +  '_Amon_ACCESS-ESM1-5_'
        if exptnamein == 'LP':
            #filename = filename + 'LP_09020116-10011216.nc'
            filename = filename + 'LP_09021231-10011231.nc'
        if exptnamein == 'PI':
            filename = filename + 'pre-industrial_01010116-02001216.nc'
        fielduse = ACCESS_FIELDS.get(fieldnamein)

    if modelname == 'AWI-CM2':
        
        filename = (filestart + 'AWI/AWI-CM2/' + exptnamein +
                    '/AWI-CM2_' + exptnamein + '_' +
                    AWI_FILE_FIELDS.get(fieldnamein,fieldnamein))
        if fieldnamein == 'tos':
            filename = filename + '_monthly_timeseries_r360x181_'
        else:            
            filename = filename + '_monthly_timeseries_native_'

        filename = filename + '580001-589912.nc'
        fielduse = AWI_FIELDS.get(fieldnamein)
                    
    if modelname == 'AWI-ESM2':
        filename = (filestart + 'AWI/AWI-ESM2/' + exptnamein +
                    '/AWI-ESM2_' + exptnamein + '_' +
                    AWI_FILE_FIELDS.get(fieldnamein,fieldnamein))
        if fieldnamein == 'tos':
            filename = filename + '_monthly_timeseries_r360x180_'
        else:            
            filename = filename + '_monthly_timeseries_native_'

        if exptnamein == 'LP':
            filename = filename + '603201-654912.nc'
        if exptnamein == 'PI':
            filename = filename + '645001-694912.nc'
       
        fielduse = AWI_FIELDS.get(fieldnamein)
                    
    if modelname == 'COSMOS':
        filename = (filestart + 'AWI/COSMOS/' + exptnamein +
                    '/COSMOS_' + exptnamein + '_' +
                    AWI_FILE_FIELDS.get(fieldnamein,fieldnamein))
        if fieldnamein == 'tos':
            filename = filename + '_monthly_timeseries_r360x180_'
            if exptnamein == 'LP':
                filename = filename + '2800-2899.nc'
            if exptnamein == 'PI':
                filename = filename + '2800-2899.nc'

        else:            
            filename = filename + '_monthly_timeseries_native_'
            if exptnamein == 'LP':
                filename = filename + '280001-289912.nc'
            if exptnamein == 'PI':
                filename = filename + '280001-289912.nc'
      

       
        fielduse = AWI_FIELDS.get(fieldnamein)
                    
   
    if modelname   ==  'MIROC4m':
        filename = (filestart+'MIROC/MIROC4m/' + exptnamein + '/mo-timeseries/'
                    + fieldnamein + '_MIROC4m_pliomip3-' + exptnamein +
                    '_mo-timeseries.nc')
        fielduse = MIROC_FIELDS.get(fieldnamein)
        
    if modelname  == 'HadCM3':
        fielduse = HadCM3_FIELDS.get(fieldnamein)
        filename = (filestart+'HadCM3/' + exptnamein + 
                    '/monthly_timeseries/HadCM3_'+exptnamein+'_'
                    +HadCM3_EXPT.get(exptnamein) + '_' + fieldnamein + '_' +
                    'monthly_timeseries.nc')
        
    if modelname  == 'MRI-CGCM2.3':
        exptuse = exptname_l.get(exptnamein)
        fielduse = MIROC_FIELDS.get(fieldnamein)
        filename = (filestart+modelname+'/'+fielduse+'/'
                      +exptuse+'.'+fielduse+'.')
    if modelname  == 'EC-Earth3.1':
        exptuse = exptname_l.get(exptnamein)
        fielduse = ECearth_FIELDS.get(fieldnamein)
        print(fielduse)
        filename = (filestart+'EC-Earth3.1/'
                  +ECearth_EXPT.get(exptnamein)
                  +'.EC-Earth3.1.surface.nc')
    if modelname  == 'IPSLCM5A' or modelname  == 'IPSLCM5A2':
        exptuse = exptname_l.get(exptnamein)
        if modelname  == 'IPSLCM5A':
            timeuse = IPSLCM5A_TIME.get(exptnamein)
        if modelname  == 'IPSLCM5A2':
            timeuse = IPSLCM5A21_TIME.get(exptnamein)
        fielduse = IPSLCM5A_FIELDS.get(fieldnamein)
        if fieldnamein  == 'tos':
            filename = (filestart+modelname+'/'
                  +IPSLCM5A_EXPT.get(exptnamein)+'.'
                  +fielduse+'_'+timeuse+'_monthly_TS_rectilinear.nc')
        else:
            filename = (filestart+modelname+'/'
                  +IPSLCM5A_EXPT.get(exptnamein)+'.'
                  +fielduse+'_'+timeuse+'_monthly_TS.nc')

    if modelname  == 'NorESM1-F' or modelname  == 'NorESM-L':
        fielduse = NorESM_FIELDS.get(fieldnamein)
        filename = (filestart+modelname+'/'+modelname+'_'+
                 exptnamein+'_'+fielduse+'.nc')
    if modelname  == 'IPSLCM6A':
        fielduse = MIROC_FIELDS.get(fieldnamein)
        if fieldnamein  == 'tos':
            filename = (filestart+modelname+'/'+fielduse+
                  '_Omon_IPSL-CM6A-LR_'+IPSLCM6A_TIME_ALT.get(exptnamein)+'_rectilinear.nc')
        else:
            print(filestart, modelname, fielduse, IPSLCM6A_TIME.get(exptnamein), atm_ocn_ind.get(fieldnamein))
            filename = (filestart+modelname+'/'+fielduse+
                  '_'+atm_ocn_ind.get(fieldnamein)+'_IPSL-CM6A-LR_'+IPSLCM6A_TIME.get(exptnamein)+'.nc')
    if modelname  == 'GISS':
        fielduse = fieldnamein
        exptuse = exptname_l.get(exptnamein)
        filename = []
        filename.append(filestart+modelname+'/'+exptuse+'/'+fielduse+
                  '_'+atm_ocn_ind.get(fieldnamein)+
                  '_GISS-E2-1-G_'+GISS_TIME1.get(exptnamein)+'.nc')
        filename.append(filestart+modelname+'/'+exptuse+'/'+fielduse+
                  '_'+atm_ocn_ind.get(fieldnamein)+
                  '_GISS-E2-1-G_'+GISS_TIME2.get(exptnamein)+'.nc')


    if modelname == 'CCSM4-2deg':
        filename=(filestart + 'CESM1.0.5/' + exptnamein + '/' +
                  exptnamein + '_' + CCSM42_FIELDS.get(fieldnamein) +
                  '.nc')
        fielduse = fieldnamein
        
    if modelname == 'iCESM1.3':
        fileinit = filestart + 'CESM/iCESM1.3/' + exptnamein
        if fieldnamein == 'tos':
            fileinit = fileinit + '/OCN'
        else:
            fileinit = fileinit + '/ATM/b.e13.B1850C5'
        
        if exptnamein == 'LP':
            fileinit = fileinit + 'CN.ne30_g16.icesm13_ihesp.'
            fileinit = fileinit + 'midPliocene.001.cam.h0.'
        if exptnamein == 'PI':
            fileinit = fileinit + '.ne30_g16.icesm13_ihesp.'
            fileinit = fileinit + 'PI.002.cam.h0.'
        if fieldnamein == 'pr':
            # this has been passed in two files. 
            # for convective and large scale precipitaiton put filename in a list
            filename1 = (fileinit + 'PRECC')
            filename2 = (fileinit + 'PRECL')
            filename = [filename1, filename2]
            fielduse = ['PRECC', 'PRECL']
        if fieldnamein == 'tos':
            # if we are passing SST we also need to pass ice fraciton
            filename1 = (filestart + 'NCAR/' + 
                  CESM12_EXTRA.get(exptnamein) + 
                  'TS' +
                  CESM12_TIME.get(exptnamein) + '.nc')
            filename2 = (filestart + 'NCAR/' + 
                  CESM12_EXTRA.get(exptnamein) + 
                  'ICEFRAC' +
                  CESM12_TIME.get(exptnamein) + '.nc')
            filename = [filename1, filename2]
            fielduse = ['TS', 'ICEFRAC']
        if fieldnamein =='tas':
            filename=(filestart + 'NCAR/' + 
                  CESM12_EXTRA.get(exptnamein) + 
                  CESM12_FIELDS.get(fieldnamein) +
                  CESM12_TIME.get(exptnamein) + '.nc')
            print(filename)
            print(fieldnamein)
            fielduse = CESM12_FIELDS.get(fieldnamein,fieldnamein)

            
            
      
    retdata = [fielduse, filename]
    return(retdata)
 

##########################################################
# main program

filename  =  ' '
linux_win  =  'l'
modelname  = "HadCM3"  # 'ACCESS' 'AWI-CM2' 'AWI-ESM2' 'COSMOS' 'MIROC4m'
                        # 'iCESM1.3'
exptname  =  {
        "E280" : "E280",
        "Eoi400" : "EOI400",
        "E400":"E400",
        "E560": "E560"}

exptname_l  =  {
        "E280" : "e280",
        "Eoi400" : "eoi400",
        "E400":"e400",
        "E560": "e560"}

fieldname  =  {
        "pr" : "TotalPrecipitation",
        "tas" : "NearSurfaceTemperature",
        "sic" : "SeaIceConcentration",
        "tos": "SST"
        }


# this is regridding where all results are in a single file
#fieldnamein = ['pr']
#exptnamein = ['Eoi400']

#fieldnamein = ['pr']
fieldnamein = ['tos'] # ocean tempeature or sst
#exptnamein = ['Eoi400']

fieldnamein = ['pr','tas']
#exptnamein = ['PI','LP']
exptnamein = ['LP_min_LSM_change']
if linux_win  == 'l':
    filestart = '/uolstore/Research/b/b0164/Data/PlioMIP3/'
else:
    filestart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'




for expt in range(0, len(exptnamein)):
    for field in range(0, len(fieldnamein)):

      
        # call program to get model dependent names
        # fielduse,  and  filename
        retdata = getnames(modelname, filestart, fieldnamein[field], exptnamein[expt])

        print(retdata)
        fielduse = retdata[0]
        filename = retdata[1]

        fieldnameout = fieldname.get(fieldnamein[field])
        exptnameout = exptname.get(exptnamein[expt],exptnamein[expt])

        print('filename is', filename)


        regrid_data(fieldnamein[field], exptnamein[expt])

#sys.exit(0)
