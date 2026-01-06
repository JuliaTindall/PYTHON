# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
@author: earjcti

This program will create a netcdf file which contains the vegetation data for PlioMIP2
"""

import iris
import numpy as np


def overwrite_data(arr1, arr2):
    """
    overwrites the data in arr1 with arr2 where the following conditions are met
    1. arr1 is not a sea point (level 1 lt 2.0)
    2. arr1 is not ice (level 9/8 not equal 1)
    3. arr1 is not a partial lake (level 7/6 = 0)
    """
    
    nt, nz, ny, nx = np.shape(arr1)
    for j in range(0, ny):
        for i in range(0, nx):
            if arr1[0, 0, j, i] < 2.0: # not sea
                if arr1[0, 8, j, i] != 1.0: # not ice
                    if arr1[0, 6, j, i] == 0.0: # is not lake
                        vegdata = overwrite_nonlake(arr2[0, :, j, i])
                    else: # partial lake
                        vegdata = overwrite_lake(arr1[:, :, j, i], arr2[:, :, j, i])



    
def main():
    """
    1. loads data in from SJH file
    2. loads data in from my file
    3. overwrite SJH data with my data where appropriate
    4. write out to a new file
    """

    sjhcube = iris.load_cube(STEVEH_FILE, STEVEH_FIELD)
    sjhdata = sjhcube.data

    jctcube = iris.load_cube(PRISM3_DUMP, PRISM3_FIELD)
    jctdata = jctcube.data
    
    if np.shape(sjhdata) != np.shape(jctdata):
        print('shapes of data do not match')
        sys.exit(0)
        
    overwrite_data(sjhdata, jctdata)
    
    


"""
data required
"""

STEVEH_FILE = '/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrfrac.type_sjh.nc'
STEVEH_FIELD = 'FRACTIONS OF SURFACE TYPES'
PRISM3_DUMP = '/nfs/hera1/earjcti/um/xibol/netcdf/xibola@daz00c1.nc'
PRISM3_FIELD = 'TILE FRACTIONS (B.LAYER)'

main()