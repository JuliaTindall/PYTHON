#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 3rd 2019


#
# This program will regrid some of the data that is needed for PLIOMIP2.
# We will put 100 year average fields onto a 1deg X 1deg standard grid
# it can be used where experiments have been uploaded with 100 years in
# one file
#
# it can currently do MIROC4 and COSMOS


import numpy as np
import iris
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import iris.analysis.cartography
import iris.coord_categorisation
import cf_units as unit
from iris.experimental.equalise_cubes import equalise_attributes
import sys
#import os

###################################################
# get all data from files as a single cube
##################################################
def get_ecearth_cube(exptname,lsmstart):
    """
    get's the datacube from ecearth
    """
    allcube = iris.load(filename)
    ncubes = len(allcube)
    for i in range(0, ncubes):
        if allcube[i].var_name == fielduse:
            cube = allcube[i]
    if fielduse  == "sst":
        if exptname  == 'Eoi400':
            lsmfile = lsmstart+modelname+'/lsm_plio_new.nc'
        if exptname  == 'E280':
            lsmfile = lsmstart+modelname+'/lsm.nc'
      
        lsmcube = iris.load(lsmfile)
        #  qiong has said to take >0.5 as land and >0.5 as sea

        nt, nx, ny = np.shape(cube.data)
        mymask_2d = lsmcube[0].data
        mymask_2d = np.where(mymask_2d > 0.5,  1.0, 0.0)
        mymask = np.vstack([mymask_2d] * nt)


        cube.data  =  np.ma.array(cube.data,  mask = mymask)
    return(cube)

def get_hadcm3_cube(model):
    """
    get's the datacube from hadcm3 or mri-cgcm2.3
    """
    allcubes = iris.cube.CubeList([])
    startyear = 0
    endyear = 100
    if model  == 'MRI-CGCM2.3':
        startyear = startyear+1
        endyear = endyear+1

    for i in range(startyear, endyear):
        yearuse = str(i).zfill(3)
        filenameuse = (filename+yearuse+'.nc')
        cubetemp = iris.load_cube(filenameuse)

        u  =  unit.Unit('days since 0800-01-01 00:00:00',
                                  calendar = unit.CALENDAR_360_DAY)
        if model  == 'HadCM3':
            cubetemp.coord('t').rename('time')
        cubetemp.coord('time').points = (np.arange(0, 12)+((i-startyear)*12))*30.

        cubetemp.coord('time').units = u

        allcubes.append(cubetemp)


    equalise_attributes(allcubes)
    cube_temp = allcubes.concatenate_cube()

    if model  == 'MRI-CGCM2.3':
        cube_temp.coord('pressure level').rename('surface')

    if model  == 'HadCM3' and fielduse  == 'SST':
        cube_temp.coord('unspecified').rename('surface')

    if model  == 'HadCM3' and fielduse  == 'NearSurfaceTemperature':
        cube_temp.coord('ht').rename('surface')


    cube_temp.coord('surface').points = 0.
    cube  =  cube_temp.extract(iris.Constraint(surface = 0.))


    return(cube)

def get_ipslcm5_atm(exptname, fieldnamein):
    """
      get data from the atmospheric files in ipslcm5

      there is a bit of an error in the file calendar so we will
    """
    # copy the data to a new file but without the error
    with Dataset(filename) as src,  Dataset("temporary.nc",  "w", format = 'NETCDF3_CLASSIC') as dst:
        # copy attributes
        for name in src.ncattrs():
            dst.setncattr(name,  src.getncattr(name))
        # copy dimensions
        for name,  dimension in src.dimensions.iteritems():

            if name !=  'tbnds':   # don't copy across time counter bounds
                dst.createDimension(name,  (len(dimension)))

        # copy all file data
        for name,  variable in src.variables.iteritems():
            print('name is', name, variable)
            if name != 'time_counter_bnds' and name!= 'time_centered':
                x  =  dst.createVariable(name,  variable.datatype,
                                       variable.dimensions)
                if name  == 'time_counter':
                    # convert from seconds to days and start at middle of month
                    dst.variables[name][:]  =  (src.variables[name][:] / (60.*60.*24))-(src.variables[name][0] / (60.*60.*24))+15.
                else:
                    dst.variables[name][:]  =  src.variables[name][:]
                # copy attributes for this variable
                for ncattr in src.variables[name].ncattrs():
                    attribute = src.variables[name].getncattr(ncattr)
                    print(ncattr, exptname)
                    if ncattr  == 'calendar' and exptname  == 'Eoi400':
                        dst.variables[name].setncattr(ncattr, '360_day')
                    else:
                        if (ncattr  == 'units' and name  == 'time_counter'):
                    # change units from seconds to days
                            dst.variables[name].setncattr(ncattr, attribute.replace('seconds', 'days'))
                        else:
                            dst.variables[name].setncattr(ncattr, attribute)

        fieldreq = fieldnamein
        if fieldnamein  == 'pr':
            fieldreq = 'Precip Totale liq+sol'
        if fieldnamein  == 'tas':
            fieldreq = 'Temperature 2m'


        cube = iris.load_cube('temporary.nc', fieldreq)

        if fieldnamein  == 'ts' or fieldnamein  == 'tas':
            cube.convert_units('Celsius')

        if exptname  == 'Eoi400':
            u  =  unit.Unit('days since 0800-01-01 00:00:00',
                              calendar = unit.CALENDAR_360_DAY)
        else:
            u  =  unit.Unit('days since 0800-01-01 00:00:00',
                          calendar = unit.CALENDAR_365_DAY)
        cube.coord('time').units = u

        return(cube)

def get_ipslcm6():
    """
    get ipslcm6 data
    here 200 years have been supplied.  We only want the last 100 years
    """
    cubeall = iris.load_cube(filename)

    cubelist = iris.cube.CubeList([])
    for i,  t_slice in enumerate(cubeall.slices(['latitude', 'longitude'])):
        if i >= 1200:
            t_slice.coord('time').bounds = None
            t_slice2 = iris.util.new_axis(t_slice, 'time')
            cubelist.append(t_slice2)

    cube = cubelist.concatenate_cube()
    return(cube)

def get_miroc_tos():
    """
    get miroc data
    """
    cube = iris.load_cube(filename)

    # MIROC didn't have the units set up correctly for 'tos'
    cube.coord('latitude').units = 'degrees'
    cube.coord('longitude').units = 'degrees'
    return(cube)

def get_giss():
    """
    get giss data
    here are multiple files containing the data
    """
    allcubes = iris.cube.CubeList([])
    for file in range(0, len(filename)):
        cubetemp = iris.load_cube(filename[file])
        allcubes.append(cubetemp)

    equalise_attributes(allcubes)

    cube = allcubes.concatenate_cube()

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


def get_cesm12_singlecube(filename_, fielduse_):
    """
    get a single cube from a single file
    this is needed because total precipitation is the sum of two cubes
    """
    allcube = iris.load(filename_)
    ncubes = len(allcube)
    for i in range(0, ncubes):
        print(allcube[i].var_name,fielduse_)
        if allcube[i].var_name == fielduse_:
            singlecube = allcube[i]
            
    return singlecube
    
def get_cesm12(exptnamein):
    """
    get cesm1.2 data
    """
    
    if type(fielduse) is list: # perhaps a list containing large scale and convective precip
                               # that need to be added together for total precipitaiton
       cube1 = get_cesm12_singlecube(filename[0], fielduse[0])
       cube2 = get_cesm12_singlecube(filename[1], fielduse[1])
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
       or modelname == 'CESM1.2' or modelname == 'CCSM4-1deg'):
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



    # outfile
    if linux_win  == 'l':
        outstart = ('/nfs/hera1/earjcti/regridded/'+modelname+'/'+exptnameout+'.'+
        fieldnameout+'.')
        lsmstart = '/nfs/hera1/earjcti/regridded/'
    else:
        outstart = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
              +modelname+'\\'+exptnameout+'.'+fieldnameout+'.')
        lsmstart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'


   
    #####################################
    # get all data in a single cube
    if modelname  == 'EC-Earth3.1': # all fields in one file
        cube = get_ecearth_cube(exptnamein,lsmstart)
    elif (modelname  == 'HadCM3' or modelname  == 'MRI-CGCM2.3'):
        cube = get_hadcm3_cube(modelname)
    elif ((modelname  == 'IPSLCM5A' or modelname  == 'IPSLCM5A2') and
        (fieldnamein != 'tos')):
        cube = get_ipslcm5_atm(exptnamein, fieldnamein)
    elif (modelname  == 'IPSLCM6A'):
        cube = get_ipslcm6()
    elif (modelname  == 'MIROC4m' and fieldnamein  == 'tos'):
        cube = get_miroc_tos()
    elif (modelname  == 'GISS'):
        cube = get_giss()
    elif (modelname  == 'CCSM4-2deg'):
        cube = get_ccsm4_2deg()
    elif (modelname  == 'CESM1.2' or modelname == 'CCSM4-1deg'):
        cube = get_cesm12(exptnamein)
    else:
        cube = iris.load_cube(filename)


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
        cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
        regridded_cube = cube.regrid(cubegrid, iris.analysis.Linear())



    refdate = 'days since 0800-01-01 00:00:00'

    # for cosmos
    if modelname  == 'COSMOS':
        # cosmos data is in a strange time coordinate line yyyymmdd
        # we need to convert it to days since reference time
        origpoints = regridded_cube.coord('time').points
        npoints = len(origpoints)
        print(npoints)

        print(origpoints)
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
    if modelname  == 'EC-Earth3.1':
    # convert from hours to days
        origpoints = regridded_cube.coord('time').points
        newpoints = origpoints/24.
        regridded_cube.coord('time').points = newpoints
        refdate = 'days since 2390-01-01 00:00:00'


    # regrid to mm/day from kg/m2/s if required
    if (modelname  == 'EC-Earth3.1' or modelname  == 'IPSLCM5A'
             or modelname  == 'IPSLCM5A2' or modelname  == 'IPSLCM6A'
             or modelname == 'CCSM4-2deg' or modelname =='GISS'):
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
        or modelname  == 'HadCM3' or modelname == 'GISS'):
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
    elif (modelname  == 'COSMOS' or modelname  == 'MIROC4m' or
          modelname  == 'IPSLCM6A' or modelname  == 'EC-Earth3.1'):
        regridded_cube.coord('time').units = refdate



       # end of Uof T loop

    print('julia2')
    #print(regridded_cube.coord('time'))

    # add auxillary coordinates month and year
    iris.coord_categorisation.add_month_number(regridded_cube,  'time',  name = 'month')
    iris.coord_categorisation.add_year(regridded_cube,  'time',  name = 'year')
    
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

    plt.plot(tempcube.data)
    plt.plot([0, 100], [meanann, meanann])
    plt.plot([0, 100], [meanann+stdevann+stdevann, meanann+stdevann+stdevann])
    plt.plot([0, 100], [meanann-stdevann-stdevann, meanann-stdevann-stdevann])
    plt.title('data and data+/-2sd')

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
        file1.write(np.str(np.round(meanann[0], 2))+', '+np.str(np.round(stdevann[0], 3))+'\n')
    else:
        file1.write(np.str(np.round(meanann, 2))+', '+np.str(np.round(stdevann, 3))+'\n')

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
            file1.write(np.str(i+1)+', '+np.str(np.round(meanmon[i].data[0], 2))+', '+np.str(np.round(stdevmon[i], 3))+'\n')
        else:
            file1.write(np.str(i+1)+', '+np.str(np.round(meanmon[i], 2))+', '+np.str(np.round(stdevmon[i], 3))+'\n')

    # get latitudinal means and standard deviation
    file1.write('zonal means and standard deviations \n')
    file1.write('----------------------------------------\n')
    file1.write('latitude    mean    sd  \n')
    for i in range(0, len(meanlat)):
        file1.write(np.str(mean_data.coord('latitude').points[i])+', '+np.str(np.round(meanlat[i], 2))+', '+np.str(np.round(stdevlat[i], 3))+'\n')

    file1.close()


    ########################################################
    # check that we have averaged properly.  To do this we are
    # going to plot the annual cycle of the global mean field
    # for the regridded dataset and also for each year



    #global mean
    plt.subplot(2, 2, 1) # global mean from each year
    #subcube = subcube_mean_mon.copy(data = )  # set up structure of subcube

    if cube.coord('latitude').has_bounds():
        cube.coord('latitude').bounds
    else:
        cube.coord('latitude').guess_bounds()

    if cube.coord('longitude').has_bounds():
        cube.coord('longitude').bounds
    else:
        cube.coord('longitude').guess_bounds()

    grid_areas  =  iris.analysis.cartography.area_weights(cube)
    newcube = cube.collapsed(['latitude', 'longitude'],  iris.analysis.MEAN, weights = grid_areas)
    nt = len(newcube.data)
    nyears = np.int(nt/12)
    print('nyears is', nyears)

    for i in range(0, nyears):
        tstart = i*12
        tend = (i+1)*12
        plotdata = newcube.data[tstart:tend]
        plt.plot(plotdata, color = 'r')



    # global mean from average

    grid_areas  =  iris.analysis.cartography.area_weights(mean_mon_data)
    temporal_mean  =  mean_mon_data.collapsed(['latitude', 'longitude'],  iris.analysis.MEAN, weights = grid_areas)

    plt.plot(temporal_mean.data, color = 'b', label = 'avg')
    plt.title('globavg '+fieldnamein)
    plt.legend()




    # check at 30N
    plt.subplot(2, 2, 2)

    bounds = cube.coord('latitude').bounds
    nbounds = cube.coord('latitude').nbounds
    nbounds, dummy = np.shape(bounds)
    for i in range(0, nbounds):
        if (bounds[i, 0]>= 32. >= bounds[i, 1] or bounds[i, 0]<= 32. <bounds[i, 1]):
            index = i
    #print(cube)
    if ndim>= 4:
        subcube = cube[:, :, index, :]
    else:
        subcube = cube[:, index, :]
    cube_avg_30N = subcube.collapsed(['longitude'], iris.analysis.MEAN)

    for i in range(0, nyears):
        tstart = i*12
        tend = (i+1)*12
        plotdata = cube_avg_30N.data[tstart:tend]
        plt.plot(plotdata, color = 'r')

    #mean at 30N
    slice_30N =  mean_mon_data.extract(iris.Constraint(latitude = 32))
    mean_30N = slice_30N.collapsed(['longitude'],  iris.analysis.MEAN)


    plt.plot(mean_30N.data, color = 'b', label = 'avg')
    plt.title('average at 30N by month')
    plt.legend()
    plt.show()
    plt.close()


#############################################################################
def getnames(modelname, filestart, fieldnamein, exptnamein):

# this program will get the names of the files and the field for each
# of the model

    # set up model specific dictionaries
    MIROC_FIELDS  = {"pr" : "pr",
        "tas" : "tas",
        "sic" : "SeaIceAreaFraction",
        "tos" : "tos"
        }

    COSMOS_FIELDS  = {"pr" : "TotalPrecip",
        "tas" : "NearSurfaceAirTemp",
        "sic" : "SeaIceAreaFraction",
        "tos" : "SeaSurfaceTemp"
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
    cosmos_version = {"tas": "",
                 "pr": "",
                 "tos":"_remapbil"}

    # get names for each model
    if modelname   ==  'MIROC4m':
        filename = filestart+modelname+'/'
        fielduse = MIROC_FIELDS.get(fieldnamein)
        filename = (filename+fielduse+
                      '/MIROC4m_'+exptnamein+'_'+atm_ocn_ind.get(fieldnamein)+'_'+fielduse+'.nc')
    if modelname   ==  'COSMOS':
        if linux_win  == 'l':
            filename = filestart+'/AWI/COSMOS/'
            filename = filename+exptnamein+'/'
        else:
            filename = filestart+'/COSMOS/'
        fielduse = COSMOS_FIELDS.get(fieldnamein)
        filename = (filename+exptnamein+'.'+fielduse+
                      '_CMIP6_name_'+fieldnamein+
                      '_2650-2749_monthly_mean_time_series'+
                      cosmos_version.get(fieldnamein)+'.nc')
    if modelname   ==  'UofT':
        if linux_win  == 'l':
            filename = filestart+modelname+'/'
            filename = filename+'UofT-CCSM4/'+exptnamein+'/'+atm_ocn_ind.get(fieldnamein)+'/'
        else:
            filename = filestart+'UofT-CCSM4\\'+exptnamein+'\\'
        fielduse = MIROC_FIELDS.get(fieldnamein)
        filename = (filename+fielduse+
                      '_'+atm_ocn_ind.get(fieldnamein)+'_'+exptnamein+
                      '_'+modelname+'-CCSM4_gr.nc')
    if modelname  == 'HadCM3':
        exptuse = exptname_l.get(exptnamein)
        fielduse = fieldname.get(fieldnamein)
        filename = (filestart+'LEEDS/HadCM3/'+exptuse+'/'+fielduse+'/'
                      +exptuse+'.'+fielduse+'.')
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
        
    if modelname == 'CESM1.2':
        if fieldnamein == 'pr':
            # this has been passed in two files. 
            # for convective and large scale precipitaiton put filename in a list
            filename1 = (filestart + 'NCAR/' + 
                  CESM12_EXTRA.get(exptnamein) + 
                  'PRECC' +
                  CESM12_TIME.get(exptnamein) + '.nc')
            filename2 = (filestart + 'NCAR/' + 
                  CESM12_EXTRA.get(exptnamein) + 
                  'PRECL' +
                  CESM12_TIME.get(exptnamein) + '.nc')
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
            fielduse = CESM12_FIELDS.get(fieldnamein)
            
    if modelname == 'CCSM4-1deg':
        if fieldnamein == 'pr':
            # this has been passed in two files. 
            # for convective and large scale precipitaiton put filename in a list
            filename1 = (filestart + 'NCAR/' + 
                  CCSM4_EXTRA.get(exptnamein) + 
                  'PRECC' +
                  CCSM4_TIME.get(exptnamein) + '.nc')
            filename2 = (filestart + 'NCAR/' + 
                  CCSM4_EXTRA.get(exptnamein) + 
                  'PRECL' +
                  CCSM4_TIME.get(exptnamein) + '.nc')
            filename = [filename1, filename2]
            fielduse = ['PRECC', 'PRECL']
        if fieldnamein == 'tos':
            # if we are using SST we also need to pass ice fraction
            filename1 = (filestart + 'NCAR/' + 
                  CCSM4_EXTRA.get(exptnamein) + 
                  'TS' +
                  CCSM4_TIME.get(exptnamein) + '.nc')
            filename2 = (filestart + 'NCAR/' + 
                  CCSM4_EXTRA.get(exptnamein) + 
                  'ICEFRAC' +
                  CCSM4_TIME.get(exptnamein) + '.nc')
            filename = [filename1, filename2]
            fielduse = ['TS', 'ICEFRAC']
        if fieldnamein == 'tas':
            filename=(filestart + 'NCAR/' + 
                  CCSM4_EXTRA.get(exptnamein) + 
                  CESM12_FIELDS.get(fieldnamein) +
                  CCSM4_TIME.get(exptnamein) + '.nc')
            fielduse = CESM12_FIELDS.get(fieldnamein)
            
      
    retdata = [fielduse, filename]
    return(retdata)


##########################################################
# main program

filename  =  ' '
linux_win  =  'l'
modelname  = "IPSLCM5A2" # MIROC4m  COSMOS UofT EC-Earth3.1
                   # HadCM3 MRI-CGCM2.3
                   # IPSLCM5A,  IPSLCM5A2
                   # NorESM1-F NorESM-L
                   # IPSLCM6A GISS
                   # new to this version CCSM4-2deg, CESM1.2
                   # CCSM4-1deg

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
#fieldnamein = ['tos'] # ocean tempeature or sst
#exptnamein = ['Eoi400']

fieldnamein = ['tos']
exptnamein = ['Eoi400', 'E280']
#exptnamein = ['E560']
if linux_win  == 'l':
    filestart = '/nfs/hera1/pliomip2/data/'
else:
    filestart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'




for expt in range(0, len(exptnamein)):
    for field in range(0, len(fieldnamein)):

        if ((modelname  == 'IPSLCM5A' or modelname  == 'IPSLCM5A2'
             or modelname == 'CCSM4-2deg')
            and (fieldnamein[field]  == 'tos')):
            filestart = '/nfs/hera1/earjcti/PLIOMIP2/'
        if (modelname  == 'IPSLCM6A' or modelname  == 'GISS'):
            filestart = '/nfs/hera1/earjcti/PLIOMIP2/'


        # call program to get model dependent names
        # fielduse,  and  filename
        retdata = getnames(modelname, filestart, fieldnamein[field], exptnamein[expt])

        fielduse = retdata[0]
        filename = retdata[1]

        fieldnameout = fieldname.get(fieldnamein[field])
        exptnameout = exptname.get(exptnamein[expt])




        print('filename is', filename)




        regrid_data(fieldnamein[field], exptnamein[expt])

#sys.exit(0)
