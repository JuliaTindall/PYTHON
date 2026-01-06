#NAME
#    average_d18o.py
#PURPOSE 
#
#  This program will create netcdf files containing the average precipitation
# and the average d18op.  For annual mean and monthly mean
#
#  before running you will need to create files in
# /nfs/hera1/earjcti/um/'expt'/d18op/means 
# using program PlioMIP_new/large_scale_features/regrid_HCM3_50_year_avg.py
# which itself uses data from extract_d18op.py

# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import sys
#from netCDF4 import Dataset, MFDataset
#from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



expt='xozzb'
startyear=2550
endyear=2600
filestart = '/nfs/hera1/earjcti/um/' + expt + '/d18op/means/'


# annual average
# read in 18o and 16o

fileann16o = (filestart + 'd18op_16o_'+expt + '_Annual_Average_' + str(startyear) + '_' + str(endyear) + '.nc')
fileann18o = (filestart + 'd18op_18o_'+expt + '_Annual_Average_' + str(startyear) + '_' + str(endyear) + '.nc')

cubeann16o = iris.load_cube(fileann16o)
cubeann18o = iris.load_cube(fileann18o)

# calculate d18o and convert 16o to mm/day
cubeannd18o = ((cubeann18o / cubeann16o) - 2005.2E-6) / 2005.2E-9
cubeannd18o.long_name = 'd18o_p'
cubeann16omn = cubeann16o * 60. * 60. *24.
cubeann16omn.long_name = 'Precipitation (mm/day)'

outfile_ann= (filestart + '/d18op_precip_'+expt + '_Annual_Average_' + str(startyear) + '_' + str(endyear) + '.nc')

iris.save([cubeannd18o,cubeann16omn], outfile_ann, netcdf_format='NETCDF3_CLASSIC',fill_value=1.0E20)
       

# monthly average
# read in 18o and 16o

fileann16o = (filestart + 'd18op_16o_'+expt + '_Monthly_Average_' + str(startyear) + '_' + str(endyear) + '.nc')
fileann18o = (filestart + 'd18op_18o_'+expt + '_Monthly_Average_' + str(startyear) + '_' + str(endyear) + '.nc')

cubeann16o = iris.load_cube(fileann16o)
cubeann18o = iris.load_cube(fileann18o)

# calculate d18o and convert 16o to mm/day
cubeannd18o = ((cubeann18o / cubeann16o) - 2005.2E-6) / 2005.2E-9
cubeannd18o.long_name = 'd18o_p'
cubeann16omn = cubeann16o * 60. * 60. *24.
cubeann16omn.long_name = 'Precipitation (mm/day)'

outfile_ann= (filestart + '/d18op_precip_'+expt + '_Monthly_Average_' + str(startyear) + '_' + str(endyear) + '.nc')

iris.save([cubeannd18o,cubeann16omn], outfile_ann, netcdf_format='NETCDF3_CLASSIC',fill_value=1.0E20)
       

