#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 18 14:13:50 2019

#@author: earjcti1
#
#  This program will extract fileds and put in a timeseries file
#
#

import os
import numpy as np
import scipy as sp
#import cf
import iris
from iris.cube import CubeList
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import sys
import warnings

     
    
   

##########################################################
# main program

# this is regridding where all results are in a single file
# create a dictionary with the long field names in and the field names we want
# we are also using dictionaries so that we only have to change timeperiod name
# when rerunning
            

expts = ['xqfmf','xqfmb','xqfmc','xqfmd']
#expt = 'xqfmg'
alt_expt = {'xqfmg': 'F_EP280',
            'xqfmh': 'F_EP',
            'xqfme' : 'F_LP280',
            'xqfmf' : 'F_LP',
            'xqfmb' : 'F_PI280',
            'xqfmc' : 'F_PI400',
            'xqfmd' : 'F_PI490'}
            
field = "OUTGOING SW RAD FLUX (TOA)"

for expt in expts:
    print(field)
    allcubes = iris.load('/uolstore/home/users/earjcti/hera1/um/' + expt + '/netcdf/' + expt + 'a@pd*',field)
    print('allcubes',allcubes)
    iris.util.equalise_attributes(allcubes)
    cubes = allcubes.concatenate_cube()
    print(cubes)
    

    iris.save(cubes,alt_expt.get(expt) + '_' + expt + '_' + field + '.nc')
