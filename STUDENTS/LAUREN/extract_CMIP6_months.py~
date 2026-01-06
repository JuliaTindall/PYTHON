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

def individual_model(model):
    period = 'EOI400'
    FILENAME = ('/nfs/hera1/earjcti/regridded/' 
                + model + 
            '/' + period + '.' + FIELD + '.allmean.nc')


    cube_pliocene = iris.load_cube(FILENAME, FIELD)

    period = 'E280'
    FILENAME = ('/nfs/hera1/earjcti/regridded/' + model + 
            '/' + period + '.' + FIELD + '.allmean.nc')

    cube_pi = iris.load_cube(FILENAME, FIELD)

    cube_anomaly = ((cube_pliocene - cube_pi)/cube_pi)*100

    V=np.arange(-30.0, 30.0, 5.0)
    
    cubedata = cube_anomaly.data
    new_array = np.ma.where(np.abs(cubedata) > 20, cubedata, -999)
    new_array2 = ma.masked_values(new_array, -999)
    print(new_array)
    newcube = cube_anomaly.copy(data=new_array2)

    #qplt.contourf(cube_anomaly, V, extend='both', cmap='coolwarm')
    qplt.contourf(newcube, V, extend='both', cmap='coolwarm')
             #should new_array replace cube_anomaly above??
    plt.title(model + ' % precip change EOI400-E280')
    plt.gca().coastlines()
    #plt.savefig(model + '_%_precip_change.jpg')
    plt.show()
    sys.exit(0)
    #plt.close()

    return cube_anomaly

def get_means(model, cube):
    """
    get global mean
    """
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    
    globalmean = cube.collapsed(['longitude','latitude'],
                                iris.analysis.MEAN,
                                weights = grid_areas)

    print(model, globalmean.data)

##########################################

#MODELNAMES = ['COSMOS', 'HadCM3', 'MIROC4m', 'CCSM4-UoT']
MODELNAMES = ['COSMOS']
FIELD = 'TotalPrecipitation'

for model in MODELNAMES:
    cube_anomaly = individual_model(model)
    global_means = get_means(model, cube_anomaly)

#################################################
# 1. COSMOS, HadCM3, MIROC  Eoi450 - Eoi400 anomaly
#    (optional changing limits, changing title, 
#     changing projection, changing colours (cmap - look in matplotlib documentation or iris))
#     (E400 - E280)  vs (Eoi400 - E280)
#
# 2. maybe tricky or maybe not.  Try plotting the global mean from all models on one figure.
#
# 3. latitudinal means and plot
#
# 4. [don't do this unless you can] as 1. but focus in on the north atlantic put.  Plotting SST. 
#
# 5. [very hard]  plotting some winds using quiverplot.  

# raw data /nfs/hera1/pliomip2/data/
# E400: COSMOS, LEEDS, MIROC4m, UoT
 
