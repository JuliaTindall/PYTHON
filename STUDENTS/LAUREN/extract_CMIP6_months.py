#/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08/09/2020

@author: earjcti

This will plot global mean 

"""

import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import iris
import iris.quickplot as qplt
import iris.coord_categorisation

filename = '/nfs/hera1/eeleb/CMIP6_data/cmip6_data/tas/tas_mon_UKESM1-0-LL_ssp585_r1i1p1f2_g025.nc'

cube = iris.load_cube(filename,'air_temperature')
iris.coord_categorisation.add_season(cube, 'time', name='clim_season')

first_point = cube.coord('clim_season').points[0]

if first_point != 'jja':
    print('check first month in file',filename,'firstpoint=',first_point)
    sys.exit(0)


seasmean_cube = cube.aggregated_by(['clim_season'],iris.analysis.MEAN)
print(seasmean_cube)
DJF_cube = seasmean_cube[0,:,:]
JJA_cube = seasmean_cube[2,:,:]

#sys.exit(0)
#print(seasmean_cube.coord('clim_season'))
#iris.util.promote_aux_coord_to_dim_coord(seasmean_cube, "clim_season")


iris.save(seasmean_cube,'seas_mean.nc')
iris.save(DJF_cube,'DJF_mean.nc')
iris.save(JJA_cube,'JJA_mean.nc')
