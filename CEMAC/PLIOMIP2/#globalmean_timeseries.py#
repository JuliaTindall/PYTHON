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
from scipy import stats
import sys
#import os

###################################################
# get all data from files as a single cube
##################################################
def get_ecearth_cube(exptname, lsmstart):
    """
    get's the datacube from ecearth
    """
    allcube = iris.load(filename)
    print(filename)
    ncubes = len(allcube)
    for i in range(0, ncubes):
        if allcube[i].var_name == fielduse:
            cube = allcube[i]
    if fielduse == "sst":
        if exptname == 'Eoi400':
            #lsmfile = lsmstart+modelname+'/lsm_plio_new.nc'
            lsmfile = lsmstart+modelname+'/EC-Earth3.3_mPlio_LSM.nc'
        if exptname == 'E280':
            #lsmfile = lsmstart+modelname+'/lsm.nc'
            lsmfile = lsmstart+modelname+'/EC-Earth3.3_PI_LSM.nc'

        lsmcube = iris.load(lsmfile)
        #  qiong has said to take >0.5 as land and >0.5 as sea

        nt, nx, ny = np.shape(cube.data)
        mymask_2d = lsmcube[0].data
        mymask_2d = np.where(mymask_2d > 0.5, 1.0, 0.0)
        mymask = np.vstack([mymask_2d] * nt)


        cube.data = np.ma.array(cube.data, mask=mymask)
    return cube

def get_hadcm3_cube(model):
    """
    get's the datacube from hadcm3 or mri-cgcm2.3
    """
    allcubes = iris.cube.CubeList([])
    startyear = 0
    endyear = 100
    if model == 'MRI2.3':
        startyear = startyear+1
        endyear = endyear+1

    for i in range(startyear, endyear):
        yearuse = str(i).zfill(3)
        filenameuse = (filename+yearuse+'.nc')
        cubetemp = iris.load_cube(filenameuse)

        u = unit.Unit('days since 0800-01-01 00:00:00',
                      calendar=unit.CALENDAR_360_DAY)
        if model == 'HadCM3':
            cubetemp.coord('t').rename('time')
        cubetemp.coord('time').points = (np.arange(0, 12)+((i-startyear)*12))*30.

        cubetemp.coord('time').units = u

        allcubes.append(cubetemp)


    equalise_attributes(allcubes)
    cube_temp = allcubes.concatenate_cube()

    if model == 'MRI2.3':
        cube_temp.coord('pressure level').rename('surface')

    if model == 'HadCM3' and fielduse == 'SST':
        cube_temp.coord('unspecified').rename('surface')

    if model == 'HadCM3' and fielduse == 'NearSurfaceTemperature':
        cube_temp.coord('ht').rename('surface')


    cube_temp.coord('surface').points = 0.
    cube = cube_temp.extract(iris.Constraint(surface=0.))

    return cube

def get_noresm_ocn(exptnamein, fieldnamein):
    """
    get noresm ocean things.  

    """
    
    print('need to do this')
    print('it is in monthly data and on a tripolar grid')
    sys.exit(0)

    return cube

def get_ipslcm5_atm(exptname, fieldnamein):
    """
      get data from the atmospheric files in ipslcm5

      there is a bit of an error in the file calendar so we will
    """
    # copy the data to a new file but without the error
    with Dataset(filename) as src, Dataset("temporary.nc", "w", format='NETCDF3_CLASSIC') as dst:
        # copy attributes
        for name in src.ncattrs():
            dst.setncattr(name, src.getncattr(name))
        # copy dimensions
        print(src.dimensions)
        #for name,  dimension in src.dimensions.iteritems():
        for name, dimension in src.dimensions.items():

            if name != 'tbnds':   # don't copy across time counter bounds
                dst.createDimension(name, (len(dimension)))

        # copy all file data
        for name, variable in src.variables.items():
            print('name is', name, variable)
            print('datatype',variable.datatype)
            print('dimensions',variable.dimensions)
            if name not in ('time_counter_bnds', 'time_centered'):
                x = dst.createVariable(name, variable.datatype,
                                       variable.dimensions)
                if name == 'time_counter':
                    # convert from seconds to days and start at middle of month
                    dst.variables[name][:] = ((src.variables[name][:] / (60.*60.*24))
                                              -(src.variables[name][0] / (60.*60.*24)) + 15.)
                else:
                    dst.variables[name][:] = src.variables[name][:]
                # copy attributes for this variable
                for ncattr in src.variables[name].ncattrs():
                    attribute = src.variables[name].getncattr(ncattr)
                    #print(ncattr, attribute, exptname)
                    if ncattr == 'calendar' and exptname == 'Eoi400':
                        dst.variables[name].setncattr(ncattr, '360_day')
                    else:
                        if (ncattr == 'units' and name == 'time_counter'):
                    # change units from seconds to days
                            dst.variables[name].setncattr(ncattr, attribute.replace('seconds', 'days'))
                        else:
                            if ncattr != "_FillValue":
                                dst.variables[name].setncattr(ncattr, attribute)

        fieldreq = fieldnamein
        if fieldnamein == 'pr':
            fieldreq = 'Precip Totale liq+sol'
        if fieldnamein == 'tas':
            fieldreq = 'Temperature 2m'


        cube = iris.load_cube('temporary.nc', fieldreq)

        if fieldnamein in ('ts', 'tas'):
            cube.convert_units('Celsius')

        if exptname == 'Eoi400':
            u = unit.Unit('days since 0800-01-01 00:00:00',
                          calendar=unit.CALENDAR_360_DAY)
        else:
            u = unit.Unit('days since 0800-01-01 00:00:00',
                          calendar=unit.CALENDAR_365_DAY)
        cube.coord('time').units = u

        return cube

def get_ipslcm6():
    """
    get ipslcm6 data
    here 200 years have been supplied.  We only want the last 100 years
    """
    cubeall = iris.load_cube(filename)

    cubelist = iris.cube.CubeList([])
    for i, t_slice in enumerate(cubeall.slices(['latitude', 'longitude'])):
        if i >= 1200:
            t_slice.coord('time').bounds = None
            t_slice2 = iris.util.new_axis(t_slice, 'time')
            cubelist.append(t_slice2)

    cube = cubelist.concatenate_cube()
    return cube

def get_miroc_tos():
    """
    get miroc data
    """
    cube = iris.load_cube(filename)

    # MIROC didn't have the units set up correctly for 'tos'
    cube.coord('latitude').units = 'degrees'
    cube.coord('longitude').units = 'degrees'
    return cube

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

    return cube


def get_ccsm4_2deg():
    """
    get ccsm4_2deg utrecht data
    """

    allcube = iris.load(filename)
    ncubes = len(allcube)
    for i in range(0, ncubes):
        if allcube[i].var_name == fielduse:
            cube = allcube[i]

    # put units as celcius if required
    if fielduse == 'tas':
        cube.units = 'Celsius'

    cube2 = iris.util.new_axis(cube, 'time')

    return cube

def get_ccsm4_uot(fieldnamein):
    """
    get Uof T cube (need to add a dimension)
    if precip convert to mm/day
    """

    cube = iris.load_cube(filename)
    points = (np.arange(0, 1200)*30)+15. # go for middle of month
    u = unit.Unit('days since 0800-01-01 00:00:00',
                  calendar=unit.CALENDAR_360_DAY) # put as 360 day calendar because of the way
                                               # the data was sent.

    cube.add_dim_coord(iris.coords.DimCoord(points,
                                            standard_name='time', long_name='time',
                                            var_name='time',
                                            units=u,
                                            bounds=None,
                                            coord_system=None, circular=False), 0)

    if fieldnamein == 'pr':
        cube.data = cube.data * 60. *60. *24. *1000.
        cube.name = 'Total precipitation'
        cube.long_name = 'Total precipitation'
        cube.units = 'mm/day'


    return cube

def get_cesm12_singlecube(filename_, fielduse_):
    """
    get a single cube from a single file
    this is needed because total precipitation is the sum of two cubes
    """
    allcube = iris.load(filename_)
    ncubes = len(allcube)
    for i in range(0, ncubes):
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
        if exptnamein == 'Eoi400':
            if modelname == 'CESM2':
                filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' +
                           'b.e21.B1850.f09_g17.' +
                           'PMIP4-midPliocene-eoi400.001.'+
                           'cam.h0.LANDFRAC.1101.1200.nc')
            else:
                filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' +
                           'b40.B1850.f09_g16.PMIP4-pliomip2.' +
                           'LANDFRAC.1001.1100.nc')
        if exptnamein == 'E280':
            #if modelname == 'CESM2':
            #    filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' +
           #                'b.e21.B1850.f09_g17.' +
           #                'CMIP6-piControl.001.cam.h0.'+
           #                'LANDFRAC.1300.1399.nc')
           # else:
            filelsm = ('/nfs/hera1/pliomip2/data/NCAR/' +
                       'b.e12.B1850.f09_g16.preind.' +
                       'cam.h0.LANDFRAC.0701.0800.nc')
        lsmcube = get_cesm12_singlecube(filelsm, 'LANDFRAC')

        nt, nx, ny = np.shape(cube.data)
        mymask_2d = lsmcube[0].data
        mymask_2d = np.where(mymask_2d > 0.01, 1.0, 0.0)
        mymask = np.vstack([mymask_2d] * nt)

        cube.data = np.ma.array(cube.data, mask=mymask)


    return cube



######################################################
def correct_start_month(cube):
    """
    parameters: cube
    returns: the same cube

    if month doesn't start on january we will have to change some of the
    years from the end to match those at the start to give 100 full years


    """

    if modelname in ('CCSM4', 'CESM1.2', 'CESM2','NorESM1-F','NorESM-L'):
        print('CCSM',cube.coord('month').points)
        months = cube.coord('month').points
        months = months -1
        for i, month in enumerate(months):
            if month == 0: months[i] = 12
        print(months)
        cube.coord('month').points = months
 
   


    return cube

######################################################
def cube_avg(cube):
    """
    Extract global annual averaged data from an array

    Parameters:
    cube (iris cube): A cube with montly data that we average

    Returns:
    meanyear (numpy array): the global mean of the field

    """

    meanyearcube = cube.aggregated_by('year', iris.analysis.MEAN)
    
    meanyearcube.coord('latitude').guess_bounds()
    meanyearcube.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(meanyearcube)
    
    yearglobavg_cube = (meanyearcube.collapsed(['longitude', 'latitude'],
                                       iris.analysis.MEAN, weights = grid_areas))
    global_avg = yearglobavg_cube.data
  

    return yearglobavg_cube.data




##############################################
def get_timeseries_data(fieldnamein, exptnamein):
    """
    regrid the data
    """

    # outfile
    if linux_win == 'l':
        outstart = ('/nfs/hera1/earjcti/regridded/' + modelname +
                    '/timeseries/' + exptnameout + '.' + fieldnameout + '.')
        lsmstart = '/nfs/hera1/earjcti/regridded/'
    else:
        outstart = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
                    + modelname + '\\timeseries\\' + exptnameout 
                    + '.' + fieldnameout + '.')
        lsmstart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'


    #####################################
    # get all data in a single cube
    if modelname in ('EC-Earth3.1', 'EC-Earth3.3'): # all fields in one file
        cube = get_ecearth_cube(exptnamein, lsmstart)
    elif modelname in ('HadCM3', 'MRI2.3'):
        cube = get_hadcm3_cube(modelname)
    elif modelname in ('IPSLCM5A', 'IPSLCM5A2') and fieldnamein != 'tos':
        cube = get_ipslcm5_atm(exptnamein, fieldnamein)
    elif modelname in ('NorESM1-F', 'NorESM-L') and fieldnamein == 'tos':
        cube = get_noresm_ocn(exptnamein, fieldnamein)
    elif modelname == 'IPSLCM6A':
        cube = get_ipslcm6()
    elif modelname in ('MIROC4m', 'tos'):
        cube = get_miroc_tos()
    elif modelname == 'GISS2.1G':
        cube = get_giss()
    elif modelname == 'CCSM4-Utr':
        cube = get_ccsm4_2deg()
    elif modelname in ('CESM1.2', 'CCSM4', 'CESM2'):
        cube = get_cesm12(exptnamein)
    elif (modelname == 'CCSM4-UoT'):
        cube = get_ccsm4_uot(fieldnamein)
    else:
        cube = iris.load_cube(filename)

    

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
    mean_year_array = cube_avg(regridded_cube)
    
    return mean_year_array




#############################################################################
def getnames(modelname, filestart, fieldnamein, exptnamein):

# this program will get the names of the files and the field for each
# of the model

    # set up model specific dictionaries
    MIROC_FIELDS = {"pr" : "pr",
                    "tas" : "tas",
                    "sic" : "SeaIceAreaFraction",
                    "tos" : "tos"
                    }

    COSMOS_FIELDS = {"pr" : "TotalPrecip",
                     "tas" : "NearSurfaceAirTemp",
                     "sic" : "SeaIceAreaFraction",
                     "tos" : "SeaSurfaceTemp"
                     }

    ECearth_FIELDS = {"pr" : "totp",
                      "tas" : "tas",
                      "tos" : "sst",
                      "sic" : "SeaIceAreaFraction"
                      }

    IPSLCM5A_FIELDS = {"pr" : "TotalPrecip_pr",
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

    CESM12_EXTRA = {"Eoi400": "b.e12.B1850.f09_g16.PMIP4-pliomip2.cam.h0.",
                    "E280": "b.e12.B1850.f09_g16.preind.cam.h0.",
                    }

    CESM2_EXTRA = {"Eoi400": "b.e21.B1850.f09_g17.PMIP4-midPliocene-eoi400.001.cam.h0.",
                   "E280": "b.e21.B1850.f09_g17.CMIP6-piControl.001.cam.h0.",
                   }

    CCSM4_EXTRA = {"Eoi400": "b40.B1850.f09_g16.PMIP4-pliomip2.",
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
                   "Eoi400" : ".1101.1200"
                   }

    CESM2_TIME = {"E280" : ".110001-120012",
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
                   "tos":"Omon"}
    cosmos_version = {"tas": "",
                      "pr": "",
                      "tos":"_remapbil"}

    # get names for each model
    if modelname == 'MIROC4m':
        filename = filestart + modelname + '/'
        fielduse = MIROC_FIELDS.get(fieldnamein)
        filename = (filename + fielduse + '/MIROC4m_'+exptnamein
                    + '_' + atm_ocn_ind.get(fieldnamein) + '_' + fielduse + '.nc')
    if modelname == 'COSMOS':
        if linux_win == 'l':
            filename = filestart+'/AWI/COSMOS/'
            filename = filename+exptnamein+'/'
        else:
            filename = filestart+'/COSMOS/'
        fielduse = COSMOS_FIELDS.get(fieldnamein)
        filename = (filename + exptnamein + '.' + fielduse +
                    '_CMIP6_name_' + fieldnamein +
                    '_2650-2749_monthly_mean_time_series' +
                    cosmos_version.get(fieldnamein) + '.nc')
    if modelname == 'CCSM4-UoT':
        if linux_win == 'l':
            filename = filestart + 'UofT/'
            filename = (filename + 'UofT-CCSM4/for_julia/' +
                        exptnamein + '/' + atm_ocn_ind.get(fieldnamein) + '/')
        else:
            filename = filestart + 'UofT-CCSM4\\' + exptnamein + '\\'
        fielduse = MIROC_FIELDS.get(fieldnamein)

        filename = (filename + fielduse +
                    '_' + atm_ocn_ind.get(fieldnamein) +
                    '_' + exptnamein + '_UofT-CCSM4_gr.nc')
    if modelname == 'HadCM3':
        exptuse = exptname_l.get(exptnamein)
        fielduse = fieldname.get(fieldnamein)
        filename = (filestart + 'LEEDS/HadCM3/' + exptuse + '/' + fielduse + '/'
                    + exptuse + '.' + fielduse + '.')
    if modelname == 'MRI2.3':
        exptuse = exptname_l.get(exptnamein)
        fielduse = MIROC_FIELDS.get(fieldnamein)
        filename = (filestart+'MRI-CGCM2.3/'+fielduse+'/'
                    +exptuse+'.'+fielduse+'.')
    if modelname == 'EC-Earth3.1' or modelname == 'EC-Earth3.3':
        fileend = '_surface.nc'
        if fieldnamein == 'tos':
            fileend = '_ci-sst.nc'
        exptuse = exptname_l.get(exptnamein)
        fielduse = ECearth_FIELDS.get(fieldnamein)
        filename = (filestart + modelname + '/'
                    + modelname
                    + '_'
                    + ECearth_EXPT.get(exptnamein)
                    + fileend)
    if modelname == 'IPSLCM5A' or modelname == 'IPSLCM5A2':
        exptuse = exptname_l.get(exptnamein)
        if modelname == 'IPSLCM5A':
            timeuse = IPSLCM5A_TIME.get(exptnamein)
        if modelname == 'IPSLCM5A2':
            timeuse = IPSLCM5A21_TIME.get(exptnamein)
        fielduse = IPSLCM5A_FIELDS.get(fieldnamein)
        if fieldnamein == 'tos':
            filename = (filestart + modelname + '/'
                        + IPSLCM5A_EXPT.get(exptnamein) + '.'
                        + fielduse + '_' + timeuse 
                        + '_monthly_TS_rectilinear.nc')
        else:
            filename = (filestart+modelname + '/'
                        + IPSLCM5A_EXPT.get(exptnamein) + '.'
                        + fielduse + '_' + timeuse + '_monthly_TS.nc')

    if modelname == 'NorESM1-F' or modelname == 'NorESM-L':
        fielduse = NorESM_FIELDS.get(fieldnamein)
        if fieldnamein == 'tos':
            filename = (filestart + modelname + '/' +fielduse + '/'
                    + modelname + '_' + exptnamein + '.' + fielduse + '.')
        else:
            filename = (filestart + modelname + '/' + modelname + '_' +
                    exptnamein + '_' + fielduse + '.nc')
    if modelname == 'IPSLCM6A':
        fielduse = MIROC_FIELDS.get(fieldnamein)
        if fieldnamein == 'tos':
            filename = (filestart + modelname + '/' + fielduse +
                        '_Omon_IPSL-CM6A-LR_' + IPSLCM6A_TIME_ALT.get(exptnamein) 
                        + '_rectilinear.nc')
        else:
            filename = (filestart + modelname + '/' + fielduse +
                        '_' + atm_ocn_ind.get(fieldnamein) + '_IPSL-CM6A-LR_' 
                        + IPSLCM6A_TIME.get(exptnamein) + '.nc')
    if modelname == 'GISS2.1G':
        fielduse = fieldnamein
        exptuse = exptname_l.get(exptnamein)
        filename = []
        filename.append(filestart + modelname + '/' + exptuse + '/' + fielduse +
                        '_' + atm_ocn_ind.get(fieldnamein) +
                        '_GISS-E2-1-G_' + GISS_TIME1.get(exptnamein)
                        + '.nc')
        filename.append(filestart + modelname + '/' + exptuse + '/' 
                        + fielduse +
                        '_' + atm_ocn_ind.get(fieldnamein) +
                        '_GISS-E2-1-G_' + GISS_TIME2.get(exptnamein) + '.nc')

    if modelname == 'CCSM4-Utr':
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

    if modelname == 'CESM2':
        if fieldnamein == 'pr':
            # this has been passed in two files.
            # for convective and large scale precipitaiton put filename in a list
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
        if fieldnamein =='tas':
            filename=(filestart + 'NCAR/' +
                      CESM2_EXTRA.get(exptnamein) +
                      CESM12_FIELDS.get(fieldnamein) +
                      CESM2_TIME.get(exptnamein) + '.nc')
            fielduse = CESM12_FIELDS.get(fieldnamein)

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


    retdata = [fielduse, filename]
    return retdata

def checkdrift(plio_ts, pi_ts):
    """
    this program needs a pliocene timeseries and a preindustrial timeseries
    it will calculate the linear regression and use the slope to find
    the expected drift over 100 years (in the format ??degC / century)
    
    it will do this for the pliocene the preindustrial and the anomaly
    """
    
    nyears = len(plio_ts)
    allyears = np.arange(0.0, nyears, 1.0)
    nyears_pi = len(pi_ts)
    allyears_pi = np.arange(0.0, nyears_pi, 1.0)
    nyears_min = np.min([nyears, nyears_pi])
    allyears_min = np.arange(0.0, nyears_min, 1.0)
     
    print(nyears, nyears_pi, nyears_min)
  
    print(modelname)
    print('=======')
    (slope, intercept, r_value, 
            p_value, std_err) = stats.linregress(allyears_min, plio_ts[0:nyears_min])
    print('pliocene drift = ' + np.str(np.around((slope * 100.),2)) 
          + 'dec C / centuary')
    
    
    (slope, intercept, r_value, 
            p_value, std_err) = stats.linregress(allyears_min, pi_ts[0:nyears_min])
    print('pi drift = ' + np.str(np.around((slope * 100.), 2))
          + 'dec C / centuary')
    
    anomaly = plio_ts[0:nyears_min] - pi_ts[0:nyears_min]
    allyears = np.arange(0.0, nyears_min, 1.0)
    
    (slope, intercept, r_value, 
            p_value, std_err) = stats.linregress(allyears, anomaly)
    print('anomaly drift = ' +  np.str(np.around((slope * 100.),2) )
              + 'dec C / centuary')


##########################################################
# main program

filename =  ' '
linux_win =  'l'
modelname = "CESM2" # MIROC4m  COSMOS CCSM4-UoT EC-Earth3.1
                   # HadCM3 MRI2.3
                   # IPSLCM5A,  IPSLCM5A2
                   # NorESM1-F NorESM-L
                   # IPSLCM6A GISS2.1G
                   # CCSM4-Utr, CESM1.2
                   # CCSM4
                   # new to this version
                   # EC-Earth3.3 CESM2 (b.e21)

exptname = {
        "E280" : "E280",
        "Eoi400" : "EOI400",
        "E400":"E400",
        "E560": "E560"
            }

exptname_l = {
        "E280" : "e280",
        "Eoi400" : "eoi400",
        "E400":"e400",
        "E560": "e560"
            }

fieldname = {
        "pr" : "TotalPrecipitation",
        "tas" : "NearSurfaceTemperature",
        "sic" : "SeaIceConcentration",
        "tos": "SST"
        }


# this is regridding where all results are in a single file
#fieldnamein = ['pr','tas','tos']

fieldnamein = 'tas'
#fieldnamein = ['tos'] # ocean tempeature or sst
#exptnamein = ['Eoi400']

#fieldnamein = ['tos','pr','tas']
#exptnamein = ['Eoi400', 'E280']
#exptnamein = ['E560']
if linux_win == 'l':
    filestart = '/nfs/hera1/pliomip2/data/'
else:
    filestart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'


if (modelname in ('IPSLCM5A', 'IPSLCM5A2', 'CCSM4-Utr') 
    and fieldnamein == 'tos'):
        filestart = '/nfs/hera1/earjcti/PLIOMIP2/'
if modelname in ('IPSLCM6A', 'GISS2.1G'):
    filestart = '/nfs/hera1/earjcti/PLIOMIP2/'

# call program to get model dependent names
# fielduse,  and  filename and process for eoi400
fielduse, filename = getnames(modelname, filestart, fieldnamein, 'Eoi400')
fieldnameout = fieldname.get(fieldnamein)
exptnameout = exptname.get('Eoi400')
plio_timeseries = get_timeseries_data(fieldnamein, 'Eoi400')

# call program to get model dependent names
# fielduse,  and  filename and process for e280
fielduse, filename = getnames(modelname, filestart, fieldnamein, 'E280')
fieldnameout = fieldname.get(fieldnamein)
exptnameout = exptname.get('E280')
pi_timeseries = get_timeseries_data(fieldnamein, 'E280')

checkdrift(plio_timeseries, pi_timeseries)




