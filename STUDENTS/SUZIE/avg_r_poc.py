#!/usr/bin/env python2.7
# julia 23/06/2021
# weighted average of r_dust_ancil.nc


import os
import iris
import numpy as np
import numpy.ma as ma
from iris.analysis.cartography import cosine_latitude_weights
import sys

cube = iris.load_cube('r_opal.ancil.nc')
new_cube_data = ma.masked_where(cube.data == 0, cube.data)

newcube = cube.copy(data=new_cube_data)

iris.save(newcube, 'r_opal_ancil_masked.nc', fill_value=-9999)
coslat = np.cos(np.deg2rad(newcube.coord('latitude').points))
depths = newcube.coord('unspecified_1').points


newcube.coord('latitude').guess_bounds()
newcube.coord('longitude').guess_bounds()
newcube.coord('unspecified_1').guess_bounds()


grid_areas = iris.analysis.cartography.area_weights(newcube)
grid_areas2 = np.zeros(np.shape(grid_areas))
for k, bound in enumerate(newcube.coord('unspecified_1').bounds):
    print(bound, bound[1] - bound[0],k, grid_areas.shape)
    grid_areas2[:, k, :, :] = grid_areas[:, k, :, :] * (bound[1] - bound[0])
  #  print('altered', grid_areas2[0, k, :, 0])

weighted_mean = newcube.collapsed(['longitude','latitude','unspecified_1'],
                                  iris.analysis.MEAN, weights=grid_areas2)
print('means',weighted_mean.data, np.mean(newcube.data), np.mean(cube.data))
weighted_mean_depth = newcube.collapsed(['longitude','latitude'],
                                  iris.analysis.MEAN, weights=grid_areas2)
print('weighted mean by depth',weighted_mean_depth.data, np.mean(weighted_mean_depth.data))

gridcube = newcube.copy(data = grid_areas2)
iris.save(gridcube,'grid_areas.nc')
sys.exit(0)

#
totval = 0
weighting = 0

for coord in newcube.coords():
    print(coord.name)
sys.exit(0)


print(newcube.coord('latitude').points)
print(coslat)


# get average of newcube
sys.exit(0)

cosweights = cosine_latitude_weights(newcube)

weights = cosweights


print(cosweights)
