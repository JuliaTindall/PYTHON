#!/usr/bin/env python2.7
# julia 23/06/2021
# weighted average of r_dust_ancil.nc


import os
import iris
import numpy as np
import numpy.ma as ma
from iris.analysis.cartography import cosine_latitude_weights
from iris.cube import CubeList
import sys


origcube = iris.load_cube('/nfs/see-fs-01_users/ee14s2r/MIROC/woa18_decav_s00_01.nc')

gridcube = iris.load_cube('/nfs/see-fs-01_users/ee14s2r/MIROC/to.nc')
gridcube.attributes=None
gridcube.cell_methods=None

print(origcube)
print(gridcube)

# manually equalise attributes on latitude
latorig = origcube.coord('latitude')
latgrid = gridcube.coord('latitude')
latorig.points = latorig.points.astype('float64')
latgrid.attributes = latorig.attributes
latgrid.units='degrees'


# manually equalise attributes on latitude
lonorig = origcube.coord('longitude')
longrid = gridcube.coord('longitude')
lonorig.points = lonorig.points.astype('float64')
longrid.attributes = lonorig.attributes
longrid.units='degrees'

# check the depth coord
depthorig = origcube.coord('depth')
depthgrid = gridcube.coord('depth')
depthorig.points = depthorig.points.astype('float64')
depthgrid.attributes = depthorig.attributes
depthgrid.units='meters'
print(depthorig)
print(depthgrid)

newcube = origcube.regrid(gridcube,iris.analysis.Linear())
iris.save(newcube,'woa18_decav_s00_01_regrid.nc')
