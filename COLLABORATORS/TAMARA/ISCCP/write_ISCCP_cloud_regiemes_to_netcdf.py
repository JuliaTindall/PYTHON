#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29.03.2022 by Julia

this is an experimental program to see what is in the MOSIS-CR***nc4 file
"""
import numpy as np
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import sys

  

TESTFILE = '2016.nc'
allfile = iris.load(TESTFILE)
ISCCP_cube = iris.load_cube(TESTFILE,'ws')


# read modis data because the bounds should be the same
MODIS_CR = iris.load_cube('../MODIS/MODIS_cloud_regiemes.nc')
cube_template = MODIS_CR[:,:,0:10]
cube_template.attributes =None
cube_template.coord('index of Cloud Regime').rename('index of Weather State')
print(cube_template)

# put weather states in an array
ISCCP_data = ISCCP_cube.data
print(np.shape(ISCCP_data))
print(np.shape(cube_template.data))
ws_data = np.zeros((7,6,10))
for i in range(0,10):
    ws = ISCCP_data[i,:]
    wsr =ws.reshape(6,7).T
    ws_data[:,:,i] = wsr[::-1,:]

ws_cube = cube_template.copy(data = ws_data)
iris.save(ws_cube, "ISCCP_weather_states.nc", netcdf_format="NETCDF4")


# plot all cloud regiemes to screen to see what they look like
fig = plt.figure(figsize=[12,12])
V = np.arange(0,100,10)
for i in range(0,10):
    cube = ws_cube[:,:,i]
    fig.add_subplot(4,3,i+1)
    cs = plt.pcolormesh(cube.data)
    cbar = plt.colorbar(cs)
    print(cube)
    plt.yticks(ticks=np.arange(7),labels=cube.coord('cloud top pressure dimension').points)
    plt.xticks(ticks = np.arange(6),labels=cube.coord('cloud optical thickness dimension').points,rotation=45)
    plt.title('WS:' + str(i))
plt.tight_layout()
plt.show()
              
    
