#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08/09/2020

@author: earjcti

This will plot global mean 

"""

import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.quickplot as qplt
import cartopy
import cartopy.crs as ccrs




def individual_model(model):
    # eastward winds
    FILENAME_U = ('/nfs/hera1/earjcti/regridded/' 
                + model + '/EOI400.ua_850.0.allmean.nc')

    print(FILENAME_U)
    cubeu = iris.load_cube(FILENAME_U, 'eastward_wind')


    # northward winds
    FILENAME_V = ('/nfs/hera1/earjcti/regridded/' 
                + model + '/EOI400.va_850.0.allmean.nc')

    cubev = iris.load_cube(FILENAME_V, 'northward_wind')

    # plot
    ax1 = plt.axes(projection = ccrs.PlateCarree())
    ax1.coastlines()
    #ax1.set_extent([-35, 75, 20, 90], ccrs.PlateCarree())
    ulon = cubeu.coord('longitude')
    x = ulon.points
    y = cubeu.coord('latitude').points
    u = cubeu.data
    v = cubev.data
    windspeed = (cubeu ** 2 + cubev **2) ** 0.5
     
    # plot strength of winds
    cs = qplt.contourf(windspeed,20, cmap='Blues')
    plt.title("Wind speed at 850mb")
        
  
    # Add arrows to show the wind vectors
    n=12 #plot every nth row
    scale=400
    Q = plt.quiver(x[::n], y[::n], u[::n, ::n], v[::n, ::n], pivot='middle', 
               scale=scale)
    qk = ax1.quiverkey(Q, 0.8, 1.05, 10, ' 1 m/s', labelpos='E',
                          fontproperties={'size':10})

    plt.show()    

    


##########################################

MODELNAMES = ['CESM2']

for model in MODELNAMES:
    cube_anomaly = individual_model(model)

 
