1#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.09.2019 by Julia
This program will plot the annual mean difference between two models
for a given field
Note the data must have been preprocessed by CEMAC/PLIOMIP2/regrid_HCM3_50_year_avg.py
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
from iris.cube import CubeList
import iris.quickplot as qplt
import iris.analysis.cartography
import iris.coord_categorisation
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap,Normalize


class PiecewiseNorm(Normalize):
    def __init__(self, levels, clip=False):
        # the input levels
        self._levels = np.sort(levels)
        # corresponding normalized values between 0 and 1
        self._normed = np.linspace(0, 1, len(levels))
        Normalize.__init__(self, None, None, clip)

    def __call__(self, value, clip=None):
        # linearly interpolate to get the normalized value
        return np.ma.masked_array(np.interp(value, self._levels, self._normed))

def diff_two_experiments():
    """
    plots the difference in annual mean between experiment1 and experiment2
    """
    annmean_exp1_cube = iris.load_cube(FILESTART + EXPT1 + FILEEND)
    annmean_exp2_cube = iris.load_cube(FILESTART + EXPT2 + FILEEND)
    annmean_exp1_cube = iris.util.squeeze(annmean_exp1_cube)
    annmean_exp2_cube = iris.util.squeeze(annmean_exp2_cube)

    anom_cube = annmean_exp2_cube - annmean_exp1_cube
   
    if type == 'PliominPi':
        V = np.arange(-10.0,11.0, 1.0)
    else:
        V = [-30,-15.,-10., -5., -2., -1., -0.5, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 30.]
   
    #plt.subplot(2,1,1)
    
    if FIELDNAME == 'TotalPrecipitation':
        for i,vind in enumerate(V):
            V[i] = vind / 10.
        mycmap = cm.get_cmap('rainbow_r', len(V)+2)
    else:
        mycmap = cm.get_cmap('rainbow', len(V)+2)
    newcolors = mycmap(np.linspace(0,1,len(V)+2))
    white = ([1,1,1,1])
    print('julia',int(len(V)/2),int(len(V)/2+1))
    newcolors[int((len(V)/2)):int(len(V)/2+3),:] = white
    mycmap = ListedColormap(newcolors)

    cs=iplt.contourf(anom_cube, levels=V, extend='both', cmap=mycmap,
                  norm=PiecewiseNorm(V))
    cbar = plt.colorbar(cs,orientation='horizontal',ticks=V)
    if FIELDNAME == 'TotalPrecipitation':
        cbar.set_label('mm/day')
    plt.title(FIELDNAME + ':' + EXPT2 + '-' + EXPT1 + ' annmean')
    plt.gca().coastlines()

   # plt.subplot(2,1,2)
   # qplt.contourf(anom_cube, levels=np.arange(-3,3.5,0.5), extend='both', cmap='Rd#Bu_r')
   # plt.title(FIELDNAME + ':' + EXPT2 + '-' + EXPT1 + ' annmean')
   # plt.gca().coastlines()

    fileout = ('/nfs/hera1/earjcti/HadCM3_plots/' + FIELDNAME + '/' + EXPT2 + '-' + EXPT1 + '_annmean')
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()


def diff_two_anomalies():
    """
    plots the difference in annual mean between experiment1 and experiment2
    """
    expt1e = EXPT1[0:5]
    expt1c = EXPT1[6:11]
    expt2e = EXPT2[0:5]
    expt2c = EXPT2[6:11]
  
    annmean_exp1e_cube = iris.util.squeeze(iris.load_cube(FILESTART + expt1e + FILEEND))
    annmean_exp1c_cube = iris.util.squeeze(iris.load_cube(FILESTART + expt1c + FILEEND))
    annmean_exp2e_cube = iris.util.squeeze(iris.load_cube(FILESTART + expt2e + FILEEND))
    annmean_exp2c_cube = iris.util.squeeze(iris.load_cube(FILESTART + expt2c + FILEEND))

    anom_cube = ((annmean_exp2e_cube - annmean_exp2c_cube) -
                 (annmean_exp1e_cube - annmean_exp1c_cube))

    V = [-15.,-10., -5., -2. -1., -0.5, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]

    if FIELDNAME == 'TotalPrecipitation':
        for i,vind in enumerate(V):
            V[i] = vind / 10.
        mycmap = cm.get_cmap('rainbow_r', len(V)+2)
    else:
        mycmap = cm.get_cmap('rainbow', len(V)+2)
    newcolors = mycmap(np.linspace(0,1,len(V)+2))
    white = ([1,1,1,1])
    newcolors[int((len(V)/2)):int(len(V)/2+2),:] = white
    mycmap = ListedColormap(newcolors)


    cs=iplt.contourf(anom_cube, levels=V, extend='both', cmap=mycmap,
                  norm=PiecewiseNorm(V))
    cbar = plt.colorbar(cs,orientation='horizontal',ticks=V)
   
    plt.title(FIELDNAME + ':' + EXPT2 + '-' + EXPT1 + ' annmean')
    plt.gca().coastlines()

    fileout = ('/nfs/hera1/earjcti/HadCM3_plots/' + FIELDNAME + '/' + EXPT2 + '-' + EXPT1 + '_annmean')
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

##########################################################
# main program

EXPT1 = 'tenvj-xozzz'
EXPT2 = 'xpkmc-xpkma'
#type = 'PliominPi' # type is PliominPi, PiminPi, PliominPlio
type = 'PliominPi'

FIELDNAME = 'NearSurfaceTemperature'

FILESTART = '/nfs/hera1/earjcti/um/'
FILEEND = '/' + FIELDNAME + '/means/allmean.nc'

if len(EXPT1) > 10 and len(EXPT2) > 10:
    diff_two_anomalies()
else:
     diff_two_experiments()

