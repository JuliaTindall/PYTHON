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

  

TESTFILE = 'MODIS_CR.EqualAngle_daily.C61.V1.0.L3.2017.nc4'
modis_cube = iris.load_cube(TESTFILE,'centroids_of_cloud_regimes')


print(modis_cube.coord('cloud top pressure dimension'))
cube2 = iris.util.reverse(modis_cube,'cloud top pressure dimension')

print(cube2)
print(cube2.coord('cloud top pressure dimension'))
print(cube2.coord('cloud optical thickness dimension'))
print(cube2.coord('index of Cloud Regime'))
cube2.transpose([1,2,0])
print(cube2)

iris.save(cube2, "MODIS_cloud_regiemes.nc", netcdf_format="NETCDF4")
