#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 3rd 2019


#
# This program will regrid some of the data that is needed for PLIOMIP2.
# We will put 100 year average fields onto a 1deg X 1deg standard grid
# it can be used where experiments have been uploaded with 100 years in
# one file
#
# note a similar program which calculates the means without regridding is
# noregrid_ocn_50yr_avg.py


import numpy as np
import iris
from iris.cube import CubeList
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import iris.analysis.cartography
import iris.coord_categorisation
import cf_units as unit
#from iris.experimental.equalise_cubes import equalise_attributes
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
    print(filename)
    ncubes = len(allcube)
    for i in range(0, ncubes):
        if allcube[i].var_name == fielduse:
            cube = allcube[i]
    if fielduse  == "sst":
        if exptname  == 'Eoi400':
            #lsmfile = lsmstart+modelname+'/lsm_plio_new.nc'
            lsmfile = lsmstart+modelname+'/EC-Earth3.3_mPlio_LSM.nc'
        if exptname  == 'E280':
            #lsmfile = lsmstart+modelname+'/lsm.nc'
            lsmfile = lsmstart+modelname+'/EC-Earth3.3_PI_LSM.nc'
      
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
    allcubes = CubeList([])
    startyear = 1
    endyear = 100
    if model  == 'MRI2.3':
        startyear = startyear+1
        endyear = endyear+1

    for i in range(startyear, endyear):
        yearuse = str(i).zfill(3)
        filenameuse = (filename+yearuse+'.nc')
        print(filenameuse)
        cubetemp = iris.load_cube(filenameuse)

        u  =  unit.Unit('days since 0800-01-01 00:00:00',
                                  calendar = unit.CALENDAR_360_DAY)
        if model  == 'HadCM3':
            cubetemp.coord('t').rename('time')
        cubetemp.coord('time').points = (np.arange(0, 12)+((i-startyear)*12))*30.

        cubetemp.coord('time').units = u

        allcubes.append(cubetemp)


    iris.util.equalise_attributes(allcubes)
    cube_temp = allcubes.concatenate_cube()

    cube = iris.util.squeeze(cube_temp)

    #if model  == 'MRI2.3':
    #    cube_temp.coord('pressure level').rename('surface')

    #if model  == 'HadCM3' and fielduse  == 'SST':
    #    cube_temp.coord('unspecified').rename('surface')

    #if model  == 'HadCM3' and fielduse  == 'NearSurfaceTemperature':
    #    cube_temp.coord('ht').rename('surface')


    #cube_temp.coord('surface').points = 0.
    #cube  =  cube_temp.extract(iris.Constraint(surface = 0.))


    return(cube)

def get_HadGEM3_atm(fielduse,fieldnamein):
    """
    """
    cube = iris.load_cube(filename)
    cube.coord('t').rename('time')
    print('julia',fielduse)
    if fielduse == 'Temperature T':
        cube.convert_units('Celsius')
    print(cube.data)
    cube.var_name = fieldnamein
    cube.long_name = fieldname.get(fieldnamein)
    cube.standard_name = None

    #print(cube.name, fielduse, fieldname.get(fieldnamein))
  
    return cube

def get_HadGEM3_tos(exptin, fielduse, fieldnamein):
  
    """
    here there is one file per month containing the data
    """
    months = ['01','02','03','04','05','06','07','08','09','10','11','12']
    filemid = 'o_1m_'
    fileend = '_grid-T.nc'
   
    if exptin == 'Eoi400':
    # eoi400
        startyear = 2334
        endyear = 2434
        extra = 'v963'

    if exptin == 'E280':
    #e280
        startyear=1950
        endyear = 2050
        extra='q637'
       
    #endyear=2050

    allcubes = CubeList([])
   
    for year in range(startyear, endyear):
        # eoi400
        if year >= 2394: 
            extra = 'x150'
        for i, mon in enumerate(months):
            datestart = np.str(year) + mon + '01-'
            if i == 11:
                daterange = datestart + np.str(year+1) + months[0] + '01'
            else:
                daterange = datestart + np.str(year) + months[i+1] + '01'
            file = filename + extra + filemid + daterange + fileend
            print(file, fieldname)
            
            cubetemp = iris.load(file)
            cubetemp = iris.load_cube(file, fielduse)
            u = unit.Unit('days since 0800-01-01 00:00:00',
                  calendar=unit.CALENDAR_360_DAY) # put as 360 day calendar
            cubetemp.coord('time').attributes = None
            cubetemp.coord('time').points = ((i+((year-startyear)*12))*30.)+15.
            cubetemp.coord('time').units = u
            allcubes.append(cubetemp)
       
    iris.util.equalise_attributes(allcubes)
    cube = allcubes.concatenate_cube()
    print(cube.coord('time').points)

    cube.var_name = fieldnamein
    cube.long_name = fieldname.get(fieldnamein)
    cube.standard_name = None


    return cube





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
        print(src.dimensions)
        #for name,  dimension in src.dimensions.iteritems():
        for name,  dimension in src.dimensions.items():

            if name !=  'tbnds':   # don't copy across time counter bounds
                dst.createDimension(name,  (len(dimension)))

        # copy all file data
        for name,  variable in src.variables.items():
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
                    print(ncattr, attribute,exptname)
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

    cubelist = CubeList([])
    for i,  t_slice in enumerate(cubeall.slices(['latitude', 'longitude'])):
        if i >= 1200:
            t_slice.coord('time').bounds = None
            t_slice2 = iris.util.new_axis(t_slice, 'time')
            cubelist.append(t_slice2)

    cube = cubelist.concatenate_cube()
    return(cube)

def get_noresm_400(fieldname):
    """
    get ipslcm6 data
    here 200 years have been supplied.  We only want the last 100 years
    """
    cubeall = iris.load_cube(filename, fieldname)

    cubelist = CubeList([])
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
    allcubes = CubeList([])
    for file in range(0, len(filename)):
        cubetemp = iris.load_cube(filename[file])
        allcubes.append(cubetemp)

    iris.util.equalise_attributes(allcubes)

    cube = allcubes.concatenate_cube()

    return(cube)
    
    
def get_ccsm4_2deg():
    """
    get ccsm4_2deg utrecht data
    """
    
    allcube = iris.load(filename)
    ncubes = len(allcube)
    for i in range(0, ncubes):
        print(fielduse,allcube[i].var_name)
        if allcube[i].var_name == fielduse:
            cube = allcube[i]
    
    # put units as celcius if required
    if fielduse == 'tas':
        print(cube.units)
        cube.units='Celsius'
    
    print(cube)
    sys.exit(0)
    cube2 = iris.util.new_axis(cube, 'time')
    print('julia check')
    
    return cube

def get_ccsm4_uot(fieldnamein):
    """
    get Uof T cube (need to add a dimension)
    if precip convert to mm/day
    """
    
    cube = iris.load_cube(filename)
    points = (np.arange(0, 1200)*30)+15. # go for middle of month
    u  =  unit.Unit('days since 0800-01-01 00:00:00',
               calendar = unit.CALENDAR_360_DAY) # put as 360 day calendar because of the way
                                               # the data was sent.
    
    cube.add_dim_coord(iris.coords.DimCoord(points,
                standard_name = 'time',  long_name = 'time',
                var_name = 'time',
                units = u,
                bounds = None,
                coord_system = None,  circular = False), 0)

    if fieldnamein == 'sic':
        cube.data = cube.data / 100.
   # if fieldnamein  == 'pr':
   #    cube.data = cube.data * 60. *60. *24. *1000.
   #    cube.name = 'Total precipitation'
   #    cube.long_name = 'Total precipitation'
   #    cube.units = 'mm/day'
    
    
    return cube

def get_cesm12_singlecube(filename_, fielduse_):
    """
    get a single cube from a single file
    this is needed because total precipitation is the sum of two cubes
    """
    allcube = iris.load(filename_)
    print(allcube)
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
        cube.convert_units('Celsius')

    if fieldnameout == 'SST':
        # we need to mask out the land
        if exptnamein == 'Eoi400' :
            if modelname == 'CESM2':
               filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' + 
                          'b.e21.B1850.f09_g17.' + 
                          'PMIP4-midPliocene-eoi400.001.'+
                          'cam.h0.LANDFRAC.1101.1200.nc')
            else:
                filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' + 
                           'b40.B1850.f09_g16.PMIP4-pliomip2.' + 
                           'LANDFRAC.1001.1100.nc')
        if exptnamein == 'E280' or exptnamein == 'E400':
          #  if modelname == 'CESM2':
          #      filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' + 
          #                 'b.e12.B1850.f09_g17.' +
          #                 'CMIP6-piControl.001.cam.h0.'+
          #                 'LANDFRAC.1100.1200.nc')
          #  else:
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


def reduce_years(cube100yr, fieldnamein):
    """
    this will supply a cube with 100 years of data
    we reduce this to 50 years
    """
    
    cubelist = CubeList([])
    for i,  t_slice in enumerate(cube100yr.slices(['latitude', 'longitude'])):
        if i >= 600:
            t_slice.coord('time').bounds = None
            t_slice2 = iris.util.new_axis(t_slice, 'time')
            cubelist.append(t_slice2)

    print(fieldnamein)
    if fieldnamein == 'clt':
        cube50yr = cube100yr
    else:
        cube50yr = cubelist.concatenate_cube()
    
    return cube50yr

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
    if ((modelname  == 'NorESM1-F' and exptnamein[expt] != 'E400')
        or modelname  == 'NorESM-L'
        or modelname == 'CESM1.2' or modelname == 'CCSM4'
        or modelname == 'CESM2' or modelname == 'CCSM4-Utr'):
        print('in here do we really want to reorder')
        print(modelname, exptnamein[expt])
        sys.exit(0)
      
        allcubes = CubeList([])
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

    if avg100yr == 'y':
        regridded = 'regridded100/'
    else:
        regridded = 'regridded/'


    # outfile
    if linux_win  == 'l':
        print(regridded, modelname, exptnameout, fieldnameout)
        outstart = ('/nfs/hera1/earjcti/'+ regridded +modelname+'/'+exptnameout+'.'+
        fieldnameout+'.')
        lsmstart = '/nfs/hera1/earjcti/' + regridded
    else:
        outstart = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\' + regridded
              +modelname+'\\'+exptnameout+'.'+fieldnameout+'.')
        lsmstart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'


   
    #####################################
    # get all data in a single cube
    if (modelname  == 'EC-Earth3.1' or
       modelname == 'EC-Earth3.3'): # all fields in one file
        cube100 = get_ecearth_cube(exptnamein,lsmstart)
    elif (modelname  == 'HadCM3' or modelname  == 'MRI2.3'):
        cube100 = get_hadcm3_cube(modelname)
    elif ((modelname  == 'IPSLCM5A' or modelname  == 'IPSLCM5A2') and
        (fieldnamein != 'tos')):
        cube100 = get_ipslcm5_atm(exptnamein, fieldnamein)
    elif (modelname  == 'IPSLCM6A'):
        cube100 = get_ipslcm6()
    elif (modelname  == 'MIROC4m' and fieldnamein  == 'tos'):
        cube100 = get_miroc_tos()
    elif (modelname  == 'HadGEM3' and fieldnamein  == 'tos'):
        cube100 = get_HadGEM3_tos(exptnamein, fielduse, fieldnamein)
    elif (modelname  == 'HadGEM3' and fieldnamein  != 'tos'):
        cube100 = get_HadGEM3_atm(fielduse, fieldnamein)
    elif (modelname  == 'GISS2.1G'):
        cube100 = get_giss()
    elif (modelname  == 'CCSM4-Utr'):
        cube100 = get_ccsm4_2deg()
    elif (modelname  == 'CESM1.2' 
          or modelname == 'CCSM4'
          or modelname == 'CESM2'):
        print('before',filename)
        cube100 = get_cesm12(exptnamein)
    elif (modelname == 'CCSM4-UoT'):
        cube100 = get_ccsm4_uot(fieldnamein)
    elif (modelname == 'NorESM1-F' and exptnamein == 'E400'):
        cube100 = get_noresm_400(fielduse)
    else:
        cube100 = iris.load_cube(filename)

    
    if modelname == 'COSMOS' and fieldnamein == 'sic':
        # convert from % to a fraction
        cube100.data = cube100.data / 100.
     

    ###########################################
    # reduce number of years to 50

    if avg100yr == 'y':
        cube = cube100
    else:
        print(cube100)
        cube = reduce_years(cube100, fielduse)
    cube.data = cube.data.astype('float32')
    ndim = cube.ndim


    # now regrid the cube onto a 1X1 grid (we will first try regridding the raw data)
    # we have stored the grid we want in a file 'one_lev_one_deg.nc'

    # do not need to regrid CCSM4_UoTdata or a field that was originally on a tripolar grid
   
    if ((modelname   == 'CCSM4-UoT')
        or (modelname  == 'IPSLCM5A' and fieldnamein  == 'tos')
        or (modelname  == 'IPSLCM5A2' and fieldnamein  == 'tos')):
        regridded_cube = cube
    else:
        cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
        print(cube)
        print(cubegrid)
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
    if (modelname  == 'EC-Earth3.1'):
    # convert from hours to days
        origpoints = regridded_cube.coord('time').points
        newpoints = origpoints/24.
        regridded_cube.coord('time').points = newpoints
        refdate = 'days since 2390-01-01 00:00:00'


    # regrid to mm/day from kg/m2/s if required
    if (modelname  == 'EC-Earth3.1' or modelname == 'EC-Earth3.3'
             or modelname  == 'IPSLCM5A' or modelname == 'HadGEM3'
             or modelname  == 'IPSLCM5A2' or modelname  == 'IPSLCM6A'
             or modelname == 'CCSM4-Utr' or modelname =='GISS2.1G'):
        if fieldnamein  == 'pr':
            regridded_cube.data = regridded_cube.data * 60. *60. *24.
            cube.data = cube.data* 60. *60. *24.
            regridded_cube.name = 'Total precipitation'
            regridded_cube.long_name = 'Total precipitation'
            regridded_cube.units = 'mm/day'



    if (modelname  == 'NorESM1-F' 
        or modelname  == 'NorESM-L' 
        or modelname == 'CESM1.2'
        or modelname == 'CESM2'
        or modelname == 'CCSM4'):
        print('regridded_cube.units',regridded_cube.units)
        print('j1',regridded_cube.data[:,0])
       
       # if precipitation is in m/s convert to mm/day
        if fieldnamein  == 'pr':
            regridded_cube.data = regridded_cube.data * 60. *60. *24. *1000.
            print('j2',regridded_cube.data[:,0])
            cube.data = cube.data* 60. *60. *24. *1000.
            print('j3',regridded_cube.data[:,0])
            regridded_cube.name = 'Total precipitation'
            regridded_cube.long_name = 'Total precipitation'
            regridded_cube.units = 'mm/day'

    if (modelname  == 'CCSM4-UoT' or modelname  == 'NorESM1-F' or modelname  == 'NorESM-L'
        or modelname  == 'IPSLCM6A' or modelname  == 'EC-Earth3.1'
        or modelname == 'EC-Earth3.3'
        or modelname  == 'IPSLCM5A' or modelname  == 'IPSLCM5A2'
        or modelname  == 'HadCM3' or modelname == 'GISS2.1G'):
         # convert to celcius
        if fieldnamein  == 'tas' or fieldnamein  == 'tos':
            regridded_cube.convert_units('Celsius')
            cube.convert_units('Celsius')


 
        
    if (modelname  == 'COSMOS' or modelname  == 'MIROC4m' or
        modelname  == 'IPSLCM6A' or 
        modelname  == 'EC-Earth3.1'):
          regridded_cube.coord('time').units = refdate


       
    print(regridded_cube.coord('time'))
    print('refdate is',refdate)
  

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

    if mean_data.coord('latitude').has_bounds():
        mean_data.coord('latitude').bounds
    else:
        mean_data.coord('latitude').guess_bounds()

    if mean_data.coord('longitude').has_bounds():
        mean_data.coord('longitude').bounds
    else:
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

    if mean_year_data.coord('latitude').has_bounds():
        mean_year_data.coord('latitude').bounds
    else:
        mean_year_data.coord('latitude').guess_bounds()

    if mean_year_data.coord('longitude').has_bounds():
        mean_year_data.coord('longitude').bounds
    else:
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


    if mean_mon_data.coord('latitude').has_bounds():
        mean_mon_data.coord('latitude').bounds
    else:
        mean_mon_data.coord('latitude').guess_bounds()

    if mean_mon_data.coord('longitude').has_bounds():
        mean_mon_data.coord('longitude').bounds
    else:
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
    #plt.show()
    plt.close()


#############################################################################
def getnames(modelname, filestart, fieldnamein, exptnamein):

# this program will get the names of the files and the field for each
# of the model

    # set up model specific dictionaries
    MIROC_FIELDS  = {"pr" : "pr",
        "tas" : "tas",
    #    "sic" : "sic",
         "sic" : "siconc",
        "tos" : "tos",
        "clt" : "clt"
        }

    COSMOS_FIELDS  = {"pr" : "TotalPrecip",
        "tas" : "NearSurfaceAirTemp",
        "sic" : "SeaIceAreaFraction",
        "tos" : "SeaSurfaceTemp"
        }

    ECearth_FIELDS  = {"pr" : "totp",
        "tas" : "tas",
        "tos" : "sst",
        "sic" : "ci"
        }

    IPSLCM5A_FIELDS  = {"pr" : "TotalPrecip_pr",
        "tas" : "NearSurfaceTemp_tas",
        "sic" : "Sicfraction_sic",
        "tos": "SeasurfaceTemp_sst"
        }

    NorESM_FIELDS = {"pr" : "PRECT",
        "tas" : "TREFHT",
        "sic" : "aice",
        "tos" : "sst",
        "tasE400" : "Reference height temperature",
        "prE400" : "Total (convective and large-scale) precipitation rate (liq + ice)",
        "tosE400" : "Ocean surface temperature"
        }
    
    CCSM42_FIELDS = {"pr" : "TotalPrecip",
                      "tas" : "NearSurfaceAirTemp",
                      "sic" : "siconc",
                      "tos" : "SeaSurfaceTemp",
                     "clt" : "clt_Amon_CESM1.0.5_"
                      }
    
    CESM12_FIELDS = {"pr" : "TotalPrecip",
                     "tas" : "TREFHT",
                     "sic" : "ICEFRAC",
                     "tos" : "TS"
                     }

    HadGEM3_FILEFIELD = {"tas" : "airtemp",
                         "pr" : "precip",
                         "clt" : "totalcloud"}

    HadGEM3_LONGFIELD = {"tas" : "Temperature T", 
                         "pr" : "Total precipitation rate",
                         "clt" : "Total cloud"}

    CESM12_EXTRA =  {"Eoi400": "b.e12.B1850.f09_g16.PMIP4-pliomip2.cam.h0.",
                     "E280": "b.e12.B1850.f09_g16.preind.cam.h0.",
                    }
    
    CESM2_EXTRA =  {"Eoi400": "b.e21.B1850.f09_g17.PMIP4-midPliocene-eoi400.001.cam.h0.",
                    "E400": "b.e21.B1850.f09_g17.CMIP6-piControl.400.cam.h0.",
                     "E280": "b.e21.B1850.f09_g17.CMIP6-piControl.001.cam.h0.",
                    }
    
    CCSM4_EXTRA =  {"Eoi400": "b40.B1850.f09_g16.PMIP4-pliomip2.",
                     "E280": "b40.B1850.f09_g16.preind.cam.h0.",
                    }

    ECearth_EXPT = {"Eoi400": "mPlio",
              "E280":"PI"
              }

    HadGEM3_EXPT = {"Eoi400" : "pliocene",
                    "E280" : "pi"}
    
    CESM12_EXPT = {"Eoi400": "PlioMIP2",
              "E280":"PI"
              }

    IPSLCM5A_EXPT = {"Eoi400": "Eoi400",
              "E280":"PI"
              }

    CESM12_TIME = {"E280" : ".0701.0800",
                   "Eoi400" : ".1101.1200"
                   }
    
    CESM2_TIME = {"E280" : ".110001-120012",
                  "E400" : ".0801.0900",
                  "Eoi400" : ".1101.1200"
                   }
    
    CCSM4_TIME = {"Eoi400" : ".1001.1100",
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
    CCSM4_UofT_TIME = {"Eoi400": "midPliocene-eoi400_r1i1p1f1_gr1_160101-170012",
              "E280":"piControl_r1i1p1f1_gr1_150101-160012",
              "E560": "abrupt-2xCO2_r1i1p1f1_gn_195101-200012"
              }
    GISS_TIME2 = {"Eoi400": "midPliocene-eoi400_r1i1p1f1_gn_310101-315012",
              "E280":"piControl_r1i1p1f1_gn_495101-500012",
              "E560": "abrupt-2xCO2_r1i1p1f1_gn_195101-200012"
              }
    atm_ocn_ind = {"tas": "Amon",
                 "pr": "Amon",
                 "tos":"Omon",
                 "sic":"Omon",
                 "clt":"Amon"}
    cosmos_version = {"tas": "",
                 "pr": "",
                 "tos":"_remapbil",
                 "sic":"_remapbil"}

    # get names for each model
    if modelname   ==  'MIROC4m':
        filename = filestart+modelname+'/'
        fielduse = MIROC_FIELDS.get(fieldnamein)
        filename = (filename+fielduse+
                      '/MIROC4m_'+exptnamein+'_'+atm_ocn_ind.get(fieldnamein)+'_'+fielduse+'.nc')
        print(filename,fielduse)

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

        print(filename)
        print(fielduse)
        #sys.exit(0)
    if modelname   ==  'CCSM4-UoT':
        fielduse = MIROC_FIELDS.get(fieldnamein)
        
        if linux_win  == 'l':
            if fieldnamein == 'sic':
                filename = (filestart + 'UofT/UofT-CCSM4/for_julia/'
                            + exptnamein + '/SImon/' + fielduse + 
                            '_SImon_' + exptnamein + '_UofT-CCSM4_gr.nc')

            else:
                filename = (filestart + 'UofT/UofT-CCSM4/' + exptnamein + 
                        '/Amon/1x1_grid/' + fielduse + '_Amon_UofT-CCSM4_'
                        + CCSM4_UofT_TIME.get(exptnamein) + '.nc')
        else:
            filename = filestart+'UofT-CCSM4\\'+exptnamein+'\\'
        
        print(fielduse)
        print(filename)
     
    if modelname  == 'HadCM3':
        exptuse = exptname_l.get(exptnamein)
        fielduse = fieldname.get(fieldnamein)
        filename = (filestart+'LEEDS/HadCM3/'+exptuse+'/'+fielduse+'/'
                      +exptuse+'.'+fielduse+'.')
    if modelname  == 'MRI2.3':
        print(fieldnamein)
        exptuse = exptname_l.get(exptnamein)
        fielduse = MIROC_FIELDS.get(fieldnamein)
        filename = (filestart+'MRI-CGCM2.3/'+fielduse+'/'
                      +exptuse+'.'+fielduse+'.')
        print(fielduse)
        print(filename)
#        sys.exit(0)
    if modelname  == 'EC-Earth3.1' or modelname == 'EC-Earth3.3':
        fileend = '_surface.nc'
        if fieldnamein == 'tos' or fieldnamein == 'sic':
            fileend = '_ci-sst.nc'
        exptuse = exptname_l.get(exptnamein)
        fielduse = ECearth_FIELDS.get(fieldnamein)
        filename = (filestart + modelname + '/'
                    + modelname 
                    + '_' 
                    + ECearth_EXPT.get(exptnamein) 
                    + fileend)
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
        if exptnamein == 'E400':
            filename = (filestart + modelname + '/CO2_400/' + 
                        'NorESM1-F_E400_TREFHT_PRECT_month.nc')
            fielduse = NorESM_FIELDS.get(fieldnamein + 'E400')
            if fieldnamein == 'tos':
                filename = (filestart + modelname + '/CO2_400/' + 
                        'NorESM1-F_E400_SST_month.nc')
        print(fielduse,filename)
        #sys.exit(0)
          
            

    if modelname  == 'IPSLCM6A':
        fielduse = MIROC_FIELDS.get(fieldnamein)
        if fieldnamein  == 'tos':
            filename = (filestart+modelname+'/'+fielduse+
                  '_Omon_IPSL-CM6A-LR_'+IPSLCM6A_TIME_ALT.get(exptnamein)+'_rectilinear.nc')
        else:
            print(filestart, modelname, fielduse, IPSLCM6A_TIME.get(exptnamein), atm_ocn_ind.get(fieldnamein))
            filename = (filestart+modelname+'/'+fielduse+
                  '_'+atm_ocn_ind.get(fieldnamein)+'_IPSL-CM6A-LR_'+IPSLCM6A_TIME.get(exptnamein)+'.nc')
    if modelname  == 'GISS2.1G':
        fielduse = fieldnamein
        exptuse = exptname_l.get(exptnamein)
        filename = []
        filename.append(filestart+modelname+'/'+exptuse+'/'+fielduse+
                  '_'+atm_ocn_ind.get(fieldnamein)+
                  '_GISS-E2-1-G_'+GISS_TIME1.get(exptnamein)+'.nc')
        filename.append(filestart+modelname+'/'+exptuse+'/'+fielduse+
                  '_'+atm_ocn_ind.get(fieldnamein)+
                  '_GISS-E2-1-G_'+GISS_TIME2.get(exptnamein)+'.nc')


    if modelname == 'CCSM4-Utr':
        filename=(filestart + 'Utrecht/CESM1.0.5/' + exptnamein + '/' +
                  CCSM42_FIELDS.get(fieldnamein) + '_Omon_CESM1.0.5_' + 
                  exptnamein + '_r1i1p1f1_gn_275001-285012.nc')
        if fieldnamein == 'sic':
            fielduse = 'siconc'
        else:
            fielduse = fieldnamein
        print(filename)
        print(fielduse)
        #sys.exit(0)
        
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
        if fieldnamein =='sic':
            filename = (filestart + 'NCAR/' + 
                  CESM12_EXTRA.get(exptnamein) + 
                  'ICEFRAC' +
                  CESM12_TIME.get(exptnamein) + '.nc')
            fielduse = 'ICEFRAC'
            
    if modelname == 'CESM2':
        if fieldnamein == 'pr':
            # this has been passed in two files. 
            # for convective and large scale precipitaiton put filename in a list
            if exptnamein == 'Eoi400' or exptnamein == 'E400':
                filename1 = (filestart + 'NCAR/' + 
                  CESM2_EXTRA.get(exptnamein) + 
                  'PRECC' +
                  CESM2_TIME.get(exptnamein) + '.nc')
                filename2 = (filestart + 'NCAR/' + 
                  CESM2_EXTRA.get(exptnamein) + 
                  'PRECL' +
                  CESM2_TIME.get(exptnamein) + '.nc')
                filename = [filename1, filename2]
                fielduse = ['PRECC', 'PRECL']
            if exptnamein == 'E280':
                filename = (filestart + 'NCAR/b.e21.B1850.f09_g17.' + 
                            'CMIP6-piControl.001.cam.h0.PRECT.110001-120012.nc')
                fielduse = 'PRECT'
        if fieldnamein == 'tos':
            # if we are passing SST we also need to pass ice fraciton
            filename1 = (filestart + 'NCAR/' + 
                  CESM2_EXTRA.get(exptnamein) + 
                  'TS' +
                  CESM2_TIME.get(exptnamein) + '.nc')
            filename2 = (filestart + 'NCAR/' + 
                  CESM2_EXTRA.get(exptnamein) + 
                  'ICEFRAC' +
                  CESM2_TIME.get(exptnamein) + '.nc')
            filename = [filename1, filename2]
            fielduse = ['TS', 'ICEFRAC']
        if fieldnamein =='tas' or fieldnamein == 'sic':
            print(exptnamein)
            filename=(filestart + 'NCAR/' + 
                  CESM2_EXTRA.get(exptnamein) + 
                  CESM12_FIELDS.get(fieldnamein) +
                  CESM2_TIME.get(exptnamein) + '.nc')
            fielduse = CESM12_FIELDS.get(fieldnamein)
        if fieldnamein =='totcloud':
            filestart='/nfs/hera1/earjcti/PLIOMIP2/CESM2/clt_Amon_CESM2_'
            fielduse = 'clt'
            if exptnamein == 'Eoi400':
                filename = (filestart + 'midPliocene-eoi400_r1i1p1f1_'+
                            'gn_015101-020012.nc')
            if exptnamein == 'E280':
                filename = (filestart +'piControl_r1i1p1f1_gn_090001-099912.nc')
            
            
    if modelname == 'CCSM4':
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
        if fieldnamein == 'sic':
            filename = (filestart + 'NCAR/' + 
                         CCSM4_EXTRA.get(exptnamein) + 
                         'ICEFRAC' + CCSM4_TIME.get(exptnamein) + '.nc')
            fielduse = 'ICEFRAC'
   
                    
    if modelname == 'HadGEM3':
        filename = []
        filestart = '/nfs/hera1/pliomip2/data/HadGEM3_new/timeseries/' 
        if fieldnamein == 'tos':
            fielduse = 'sea_surface_temperature'
            filename = (filestart + exptnamein + '/ocean/sst_sal_temp' 
                        + '/new_nemo_b')
        else:
            fielduse = HadGEM3_LONGFIELD.get(fieldnamein)
            filename = (filestart + exptnamein + '/atmos/times_hadgem3_' + 
                        HadGEM3_EXPT.get(exptnamein) + '_' +
                        HadGEM3_FILEFIELD.get(fieldnamein) + '_final.nc')
       
     
    print(fielduse, filename)
    retdata = [fielduse, filename]
    return(retdata)


##########################################################
# main program

filename  =  ' '
linux_win  =  'l'
modelname  = "NorESM1-F" # MIROC4m  COSMOS CCSM4_UoT 
                   # HadCM3 MRI-CGCM2.3
                   # IPSLCM5A,  IPSLCM5A2
                   # NorESM1-F NorESM-L
                   # IPSLCM6A GISS2.1G
                   # CCSM4-Utr, CESM1.2
                   # CCSM4
                   # EC-Earth3.3 CESM2 (b.e21)
                   # new to this version
                   # HadGEM3
                  

exptname  =  {
        "E280" : "E280",
        "Eoi280" : "EOI280",
        "Eoi350" : "EOI350",
        "Eoi400" : "EOI400",
        "Eoi450" : "EOI450",
        "Eoi560" : "EOI560",
        "E400":"E400",
        "E560": "E560"}

exptname_l  =  {
        "E280" : "e280",
        "Eoi280" : "eoi280",
        "Eoi350" : "eoi350",
        "Eoi400" : "eoi400",
        "Eoi450" : "eoi450",
        "Eoi560" : "eoi560",
        "E400":"e400",
        "E560": "e560"}

fieldname  =  {
        "pr" : "TotalPrecipitation",
        "tas" : "NearSurfaceTemperature",
        "sic" : "SeaIceConc",
        "tos": "SST",
        "clt" : "totcloud"
        }


# this is regridding where all results are in a single file
#fieldnamein = ['pr','tas','tos']
#exptnamein = ['Eoi450']
avg100yr = 'n'

#fieldnamein = ['tas']
fieldnamein = ['tos'] # ocean tempeature or sst
#exptnamein = ['Eoi400']

#fieldnamein = ['sic']
#exptnamein = ['E280','Eoi400']
exptnamein = ['E400']
if linux_win  == 'l':
    #filestart = '/nfs/b0164/Data/'
    filestart = '/nfs/hera1/pliomip2/data/'
else:
    filestart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'




for expt in range(0, len(exptnamein)):
    for field in range(0, len(fieldnamein)):

        if ((modelname  == 'IPSLCM5A' or modelname  == 'IPSLCM5A2')
            and (fieldnamein[field]  == 'tos')):
            filestart = '/nfs/hera1/earjcti/PLIOMIP2/'
        if (modelname  == 'IPSLCM6A' or modelname  == 'GISS2.1G'):
            filestart = '/nfs/hera1/earjcti/PLIOMIP2/'


        # call program to get model dependent names
        # fielduse,  and  filename
        retdata = getnames(modelname, filestart, fieldnamein[field], exptnamein[expt])

        fielduse = retdata[0]
        filename = retdata[1]
        
 
        fieldnameout = fieldname.get(fieldnamein[field])
        exptnameout = exptname.get(exptnamein[expt])
        print(fielduse,filename)
#        sys.exit(0)
        regrid_data(fieldnamein[field], exptnamein[expt])

#sys.exit(0)
