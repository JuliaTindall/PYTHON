#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29.03.2022 by Julia

We are trying to plot cloud regiemes based on Williams and Webb 2008.

We should have extracted:
1. cloud albedo
2. cloud top pressure
3. total cloud cover

use this to get cloud regiemes at a site
"""
import numpy as np
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
import matplotlib.pyplot as plt
from numpy.linalg import norm  # used in calculating euclidian distance
import sys

def get_data():
    """
    extracts the data at the nearest gridbox
    """

    ctp_cube = iris.load_cube(FILEIN,'TOTAL CLOUD TOP HEIGHT (KFT)')
    tcc_cube = iris.load_cube(FILEIN,'TOTAL CLOUD AMOUNT - RANDOM OVERLAP')

    # find nearest longitude and latitude
    lonix = np.abs(ctp_cube.coord('longitude').points - SITELON).argmin()
    latix = np.abs(ctp_cube.coord('latitude').points - SITELAT).argmin()
   
    tcc = iris.util.squeeze(tcc_cube[:, :, latix,lonix])
    ctp = iris.util.squeeze(ctp_cube.data[:, :, latix,lonix])

    # calculate albedo
    SW_IN_TOP_cube = iris.util.squeeze(iris.load_cube(
            FILEIN, 'INCOMING SW RAD FLUX (TOA): ALL TSS'))
   
    SW_DOWN_SURF_cube = iris.util.squeeze(iris.load_cube(
            FILEIN, 'TOTAL DOWNWARD SURFACE SW FLUX'))
    
    cloud_albedo_cube_data = 1.0 - (SW_DOWN_SURF_cube.data / SW_IN_TOP_cube.data)
    cloud_albedo_cube = SW_DOWN_SURF_cube.copy(data=cloud_albedo_cube_data)
    albedo = cloud_albedo_cube.data[:, latix, lonix]

    swintop = SW_IN_TOP_cube[:, latix, lonix]
    swdownsurf = SW_DOWN_SURF_cube[:, latix, lonix]

  
    return albedo, tcc, ctp, swintop, swdownsurf
##########################################################################
def plot_timeseries(albedo, tcc, ctp, swintop, swdownsurf):
    """
    plots the timeseries of the things we are using to diagnose cloud
    """
    tccdata = tcc.data
    time = np.arange(0, len(tccdata)) / 8.0
    plt.figure(figsize=(8.0, 12.0))
    plt.subplot(4,1,1)
    plt.plot(time,tccdata)
    plt.title(MONTHNAME + ' Total cloud cover')
    plt.ylabel('frac')
    
    plt.subplot(4,1,2)
    ctpdata = np.where(ctp.data > 1E10, np.nan, ctp.data)
    plt.plot(time,ctpdata/1000)
    plt.ylim(1000, 10)
    plt.title(MONTHNAME + ' Cloud Top Pressure')
    plt.xlabel('day')
    plt.ylabel('mb')
    plt.yscale('log')
   
    albedodata = albedo.data
    plt.subplot(4,1,3)
    plt.plot(time,swintop.data, color='red')
    plt.plot(time + 1.0/8.0,swdownsurf.data, color='blue')
    plt.plot(time,swdownsurf.data, color='cyan')
    #plt.plot(time,albedo.data)
    plt.title(MONTHNAME + ' radiation')
    #plt.ylim(0,1)
    plt.xlim(0,5)
    plt.ylabel('w/m2')
   
    albedodata = np.where(swintop.data < 50, np.nan, albedodata)
    plt.subplot(4,1,4)
    plt.plot(time,albedodata)
    plt.title(MONTHNAME + ' Albedo')
    plt.ylim(0,1)
    plt.ylabel('0-1')
   
    plt.subplots_adjust()
    plt.show()
    sys.exit(0)
##########################################################################
def get_ice_free_et_clouds(observations):
    """
    input observations are [albedo, total cloud cover, cloud top pressure]
    normalised on scale 0.1 
    gets the cloud types for an ice free extra tropics site
    """
    REGIEME_NAMES_ET = {0:"Shallow cumulus", 1:"Congestus",
                  2:"Stratocu./Cu. Transition", 3:"cirrus",
                  4:"Stratocumulus", 5: "Frontal", 6:"Thin Cirrus"}
                 
    REGIEME_CHARACTERISTICS = {0:np.array([0.286, 0.643, 0.473]), 
                           1:[0.457, 0.607, 0.932],
                           2:[0.375, 0.799, 0.802],
                           3:[0.325, 0.430, 0.914],
                           4:[0.438, 0.723,0.900],
                           5:[0.581, 0.393, 0.978],
                           6:[0.220, 0.389, 0.713] }

   # norms = np.zeros(7)
   # for cloud_type in range(0, 7):
   #     chars = REGIEME_CHARACTERISTICS.get(cloud_type)
   #     norms(cloud_type) = norm(observations - chars))
                                   
         



   
#########################################################################   
#SITENAME = 'Beaver Pond'
#SITELAT = 79.0   # 79N
#SITELON = 278.0  # 82W 278E

SITENAME = 'York'
SITELAT = 54.0
SITELON = 0.0
MONTHNAME = 'July' 
shortmonth = {'July' : 'jl'}
FILEIN = 'xozzaa@pap01' + shortmonth.get(MONTHNAME) +  '.nc'

albedo, tcc, ctp, swintop, swdownsurf = get_data()
plot_timeseries(albedo, tcc,ctp, swintop, swdownsurf)

cloud_types = get_ice_free_et_clouds(np.array([albedo,tcc,ctp/100000.]))

