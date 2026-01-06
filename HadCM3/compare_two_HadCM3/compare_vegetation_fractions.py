1#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.09.2019 by Julia
This program will plot the triffid plots for each vegetation type.
It will also plot the difference between two fields
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

def get_veg(fileend):
    """
    gets the iris vegetation cube for the file
    """
    filestart = '/nfs/hera1/earjcti/um/'
    
    fileuse = filestart + fileend[0:5] + '/pi/' + fileend
    cube = iris.load_cube(fileuse,'TILE FRACTIONS (B.LAYER)')
   
    return iris.util.squeeze(cube)

def plot_veg(cube, anom_ind, fileoutend):
    """
    plots the vegetation fraction for cube.
    If anom_ind = 'n' then the cube is vegetation
    If anom_ind = 'y' then the cube is an anomaly
    """
   
    vegtype = ['Broadleaf tree','needleleaf tree','c3 grass','c4 grass','shrub','urban','water','soil','ice']
    if anom_ind == 'Y' or anom_ind == 'y':
        V = np.arange(-0.3, 0.35, 0.05)
        mycmap = cm.get_cmap('RdBu_r', len(V)+2)
        newcolors = mycmap(np.linspace(0,1,len(V)+2))
        white = ([1,1,1,1])
        print('julia',int(len(V)/2),int(len(V)/2+1))
        newcolors[int((len(V)/2)):int(len(V)/2+3),:] = white
        mycmap = ListedColormap(newcolors)
        extendreq='both'
    else:
        V = np.arange(0.0, 1.1, 0.1)
        mycmap = 'Greens'
        extendreq='neither'
   
    fig = plt.figure(figsize=[8.0, 8.0])
  
    count=0
    for i in range(0,9):
        if i in (0,1,2,3,4,7):
            count=count+1
            cubeveg = cube[i,:,:]
            plt.subplot(3,2,count)
            cs=iplt.contourf(cubeveg, levels=V, extend=extendreq, cmap=mycmap)
            plt.gca().coastlines()
            plt.title(vegtype[i])

    fig.subplots_adjust(bottom=0.2, top=0.9, left=0.1, right=0.9,wspace=0.05,hspace=0.05)
    cb_ax=fig.add_axes([0.1,0.1,0.8,0.05])
    cbar=fig.colorbar(cs,cax=cb_ax,orientation='horizontal')

    fileout = ('/nfs/hera1/earjcti/HadCM3_plots/vegetation/' + fileoutend)
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
  
    annmean_exp1e_cube = iris.load_cube(FILESTART + expt1e + FILEEND)
    annmean_exp1c_cube = iris.load_cube(FILESTART + expt1c + FILEEND)
    annmean_exp2e_cube = iris.load_cube(FILESTART + expt2e + FILEEND)
    annmean_exp2c_cube = iris.load_cube(FILESTART + expt2c + FILEEND)

    anom_cube = ((annmean_exp2e_cube - annmean_exp2c_cube) -
                 (annmean_exp1e_cube - annmean_exp1c_cube))

    V = [-15.,-10., -5., -2. -1., -0.5, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0]
    qplt.contourf(anom_cube, levels=V, extend='both', cmap="rainbow",
                  norm=PiecewiseNorm(levels))
    plt.title(FIELDNAME + ':' + EXPT2 + '-' + EXPT1 + ' annmean')
    plt.gca().coastlines()

    fileout = ('/nfs/hera1/earjcti/HadCM3_plots/' + FIELDNAME + '/' + EXPT2 + '-' + EXPT1 + '_annmean')
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

##########################################################
# main program

FILE1 = 'xozzza@pip00ja.nc'
FILE2 = 'tenvoa@piu00ja.nc'
#FILE2 = 'xozzza@piq00ja.nc'
veg_cube1 = get_veg(FILE1)
veg_cube2 = get_veg(FILE2)
veg_diff = veg_cube2 - veg_cube1

plot_veg(veg_cube1,'n',FILE1)
plot_veg(veg_cube2,'n',FILE2)
plot_veg(veg_diff,'y',FILE2 + '-' + FILE1)
