#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08/09/2020

@author: earjcti

This will plot the UCEI index 

"""


import sys
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import iris.analysis.cartography
import xarray as xr
import math as math


def low_pass_weights(window, cutoff):
    """
    # calculate weights for a running mean filter

    """
    # window: int -  The length of the filter window.

    # cutoff: float - The cutoff frequency in inverse time steps.

    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    return w[1:-1]


def get_nino_timeseries(cube, latmin, latmax, lonmin, lonmax):
   """
   gets the average timeseries from 'cube' in the region and apply a low pass filter
   """

   nino_slice = cube.extract(iris.Constraint(latitude = lambda cell: latmin <= cell <= latmax,
                                         longitude = lambda cell: lonmin <= cell <= lonmax,))
   nino_slice.coord('latitude').guess_bounds()
   nino_slice.coord('longitude').guess_bounds()
   grid_areas = iris.analysis.cartography.area_weights(nino_slice)

   nino_mean = nino_slice.collapsed(['longitude','latitude'], 
                                    iris.analysis.MEAN, weights = grid_areas)

   # apply running mean (Alice wasn't sure why you did this so have temporarily removed it))
   #window = 5
   # Construct 5 month low pass filters for the monthly mean sst data
   #wgts5 = low_pass_weights(window, 1. / 5.)
   # Apply  filter
   #Nino = nino_mean.rolling_window('time', iris.analysis.SUM, len(wgts5), weights = wgts5)
   Nino = nino_mean


   # sort out x axis
   nino_sst = Nino.data
   nmonths = len(nino_sst)
   nino_plio = np.zeros(nmonths)

   for i in range(1, nmonths-1):
        nino_plio[i] = (nino_sst[i-1] + nino_sst[i] + nino_sst[i+1]) / 3.0


 
   return nino_plio

def get_elnino_events(nino34_oni):
    """
    find out where we have el nino or la nina events
    returns event-array which is 1 for el nino and -1 for la nina
    and 0 for neutral
    """
  
    n = len(nino34_oni)
    event_array = np.zeros(n)

    for i in range(0, n):
        if nino34_oni[i] > 0.5:  # possible el nino event
            # if previous was an el nino than this one is
            if i > 0:
                if event_array[i-1] == 1.0: event_array[i] = 1.0
                pass # next loop
            # check next 4
            if i < n-5:
                if (nino34_oni[i+1] > 0.5 and nino34_oni[i+2] > 0.5
                    and nino34_oni[i+3] > 0.5
                    and nino34_oni[i+4] > 0.5):
                    event_array[i] = 1.0
        if nino34_oni[i] < -0.5:  # possible la nino event
            # if previous was an el nino than this one is
            if i > 0:
                if event_array[i-1] == -1.0: event_array[i] = -1.0
                pass # next loop
            # check next 4
            if i < n-5:
                if (nino34_oni[i+1] < -0.5 and nino34_oni[i+2] < -0.5
                    and nino34_oni[i+3] < -0.5
                    and nino34_oni[i+4] < -0.5):
                    event_array[i] = -1.0


    return event_array
                  
                
            



# define where to find data
#FILESTART = 'C:\\Users\\Team Knowhow\\OneDrive\\PlioMIP data\\CCSM4\\timeseries'
FILESTART = '/nfs/hera1/earjcti/regridded/HadCM3/timeseries/'
PERIOD = 'E280'
FIELDNAME = 'SST'
MODELNAME = 'HadCM3'
# define filename
FILENAME = (FILESTART + '/' + PERIOD + '.SST.timeseries_no_ann_cycle.nc')

# load cube
cube = iris.load_cube(FILENAME)

## need nino3 nino34 and nino4 timeseries

# limit domain to nino4 region
nino4_oni = get_nino_timeseries(cube, -5.5, 5.5, 160., 210.)
nino3_oni = get_nino_timeseries(cube, -5.5, 5.5, 210., 270.)
nino34_oni = get_nino_timeseries(cube, -5.5, 5.5, 190., 240.)

# get el nino array
elnino_array = get_elnino_events(nino34_oni) # e; nino array
                                             # will be 1 for 
                                             # el nino or -1 for la nina

# sort out x axis
nmonths = len(nino4_oni)
years = np.arange(0, nmonths, 1) / 12

# check that index looks correct
#plt.plot(years,nino34_oni)
#plt.plot(years,elnino_array, linestyle='dotted')
#plt.hlines(0.5,0,100)
#plt.hlines(-0.5,0,100)
#plt.show()
#sys.exit(0)

# normalise by standard deviation
mean_plio4 = np.mean(nino4_oni)
std_plio4 = np.std(nino4_oni)
standard_nino4_plio = ((nino4_oni - mean_plio4)/ std_plio4)

mean_plio3 = np.mean(nino3_oni)
std_plio3 = np.std(nino3_oni)
standard_nino3_plio = ((nino3_oni - mean_plio3)/ std_plio3)

# UCEI calculation


plio_N3 = standard_nino3_plio.data
plio_N4 = standard_nino4_plio.data

plio_N3_sq = np.square(plio_N3)
plio_N4_sq = np.square(plio_N4)

plio_r = np.sqrt( 2 * (plio_N3_sq + plio_N4_sq))
 
plio_real = np.add(plio_N3, plio_N4)
plio_imagine = np.subtract(plio_N3, plio_N4)
           
plio_var = np.divide(plio_imagine, plio_real)

plio_theta = np.zeros(np.shape(plio_var))

for i, x in enumerate(plio_real):
    if x >= 0:
         plio_theta[i] = np.arctan(plio_var[i])
    else:
         plio_theta[i] = np.arctan(plio_var[i]) - math.pi
       
plio_theta_deg = np.rad2deg(plio_theta)


# only plot those which are el ninos
plio_r = nino34_oni / elnino_array # may be a better indication of strength
plio_theta_deg = plio_theta_deg / np.abs(elnino_array)

## plot pliocene UCEI

plt.figure(figsize=(10,5))
plt.scatter(plio_theta_deg, plio_r, c ='k', s =1.5, marker = 'o')

plt.plot([-280, 100], [0.5, 0.5], c ='k', linestyle='--', linewidth=1.0)
plt.plot([-280, 100], [1.0, 1.0], c ='k', linestyle='--', linewidth=1.0)
plt.plot([-280, 100], [2.0, 2.0], c ='k', linestyle='--', linewidth=1.0)
plt.plot([-195, -195], [-0.0, 8.5], c ='k', linestyle='--', linewidth=1.0)
plt.plot([-165, -165], [-0.0, 8.5], c ='k', linestyle='--', linewidth=1.0)
plt.plot([-90, -90], [-0.0, 8.5], c ='k', linestyle='-', linewidth=1.0)
plt.plot([-15, -15], [-0.0, 8.5], c ='k', linestyle='--', linewidth=1.0)
plt.plot([15, 15], [-0.0, 8.5], c ='k', linestyle='--', linewidth=1.0)

plt.xlim(-270,90)
plt.ylim(0,8)

plt.title(PERIOD, fontsize = 14)
plt.xlabel('θ', fontsize = 12)
plt.ylabel('ENSO Strength', fontsize = 12)
plt.show()
