#!/usr/bin/env python3
# created 08/12/2022 by Julia
#
# This program will put all the regridded data for pliomip2 in one file 
# so they are easier for sharing.

import iris
from iris.cube import CubeList
import numpy as np

MODELNAMES = ['CESM2', 'IPSLCM6A', 'COSMOS', 
            'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
            'MIROC4m', 'IPSLCM5A2', 'HadCM3',
            'GISS2.1G', 'CCSM4', 
            'CCSM4-Utr', 'CCSM4-UoT', 
            'NorESM-L', 'MRI2.3', 'NorESM1-F','HadGEM3'
            ]

FIELD = 'SST'

EXPTNAME = 'E280'

FILEINSTART = '/nfs/hera1/earjcti/regridded/'

cubelist = CubeList([])
for model in MODELNAMES:
    filename = FILEINSTART + model + '/' + EXPTNAME + '.SST.allmean.nc'
    cube = iris.load_cube(filename)
    print(cube)
    for coord in cube.coords():
        coord.bounds = None
    cube.long_name = model + '_' + FIELD
    cube.data = np.where(cube.data.mask == True, -99999., cube.data)
    cubelist.append(cube)
    


fileout = FILEINSTART + FIELD + '_' + EXPTNAME + '_allmodels.nc'
iris.save(cubelist,fileout,fill_value = -99999.)
