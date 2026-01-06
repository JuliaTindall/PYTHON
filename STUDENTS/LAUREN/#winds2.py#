#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08/09/2020

@author: earjcti

This will plot winds as a quiverplot

"""

import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.quickplot as qplt
import cartopy
import cartopy.crs as ccrs


def individual_model(model):
    # eastward winds for pliocene and preindustrial
    period = 'EOI400'
    FILENAME_U = ('/nfs/hera1/earjcti/regridded/' 
                + model + '/EOI400.ua_850.0.allmean.nc')

    print(FILENAME_U)
    cubeU_plio = iris.load_cube(FILENAME_U, 'eastward_wind')


#  JULIA NOTE:  IF YOU XCONV THIS PI FILE YOU WILL SEE THAT THE LONG FIELD
# NAME IS UA_850.0
    period = 'E280'
    FILENAME_U = ('/nfs/hera1/earjcti/regridded/'
                  + model + '/E280.ua_850.0.allmean.nc')

    print(FILENAME_U)
    cubeU_pi = iris.load_cube(FILENAME_U, 'ua_850.0')


# JULIA NOTE YOU CANNOT DO A SUBTRACT IF THE UNITS ARE DIFFERENT SO WE
# WILL SET THE UNITS TO THE SAME
    cubeU_pi.units = 'm s-1'
    cubeU_anom = cubeU_plio - cubeU_pi

    # northward winds
    period = 'EOI400'
    FILENAME_V = ('/nfs/hera1/earjcti/regridded/' 
                + model + '/EOI400.va_850.0.allmean.nc')

    print(FILENAME_V)
    cubeV_plio = iris.load_cube(FILENAME_V, 'northward_wind')

    period = 'E280'
    FILENAME_V = ('/nfs/hera1/earjcti/regridded/'
                  + model + '/E280.va_850.0.allmean.nc')

    print(FILENAME_V)
    cubeV_pi = iris.load_cube(FILENAME_V, 'va_850.0')
    cubeV_pi.units = 'm s-1'
    

    cubeV_anom = cubeV_plio - cubeV_pi

    # plot
    ax1 = plt.axes(projection = ccrs.PlateCarree())
    ax1.coastlines()
    #ax1.set_extent([-35, 75, 20, 90], ccrs.PlateCarree())
    ulon = cubeU_anom.coord('longitude')
    x = ulon.points
    y = cubeU_anom.coord('latitude').points
    U_anom = cubeU_anom.data
    V_anom = cubeV_anom.data
    windspeed = (cubeU_anom ** 2 + cubeV_anom **2) ** 0.5
     
    # plot strength of winds
    V = np.arange(0, 5, 0.5)
    cs = qplt.contourf(windspeed,V, cmap='Blues')
    plt.title("Wind speed at 850mb in " + model)
        
  
    # Add arrows to show the wind vectors
    n=12 #plot every nth row
    scale=50
    Q = plt.quiver(x[::n], y[::n], U_anom[::n, ::n], V_anom[::n, ::n], pivot='middle', 
               scale=scale)
    qk = ax1.quiverkey(Q, 0.8, 1.05, 1, ' 1 m/s', labelpos='E',
                          fontproperties={'size':10})

    plt.show()
    #plt.savefig(model + '_850mb_winds'.jpg)
    #plt.close()

##########################################

MODELNAMES = ['CESM2']

for model in MODELNAMES:
    cube_anomaly = individual_model(model)

