#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.09.2019 by Julia

This program will plot a given field from the individual models
for either the Pliocene or the preindustrail or the difference between them

It will subtract the multimodel mean so that the differences are very clear
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.quickplot as qplt
import iris.plot as iplt
from iris.cube import CubeList


def plotfig(fig, figno, cube, figtitle, anom_ind):
     """
     plots the subfigure
     """
     if anom_ind == "n":
         cmapname = "terrain_r"
         Vuse=V
         if FIELDNAME == "TotalPrecipitation":
             Vuse = np.arange(0,21,1)
         
     else:
         Vuse=V
         if FIELDNAME == "TotalPrecipitation":
             cmapname = "RdBu"
         else:
             cmapname = "RdBu_r"
   
     fig.add_subplot(figno[0],2,figno[1])
     cs = iplt.contourf(cube, levels=Vuse,extend='both', cmap=cmapname)
     plt.title(FIELDNAME + ':' + figtitle)
     plt.gca().coastlines()
     cbar = plt.colorbar(cs,orientation='horizontal')

     return


#=====================================================================
def diff_two_experiments():
    """
    plots the difference in seasonal mean between experiment1 and experiment2
    """
    monmean_exp1_cube = iris.load_cube(FILESTART + EXPT1 + FILEEND)
    monmean_exp2_cube = iris.load_cube(FILESTART + EXPT2 + FILEEND)
    if FIELDNAME == 'TotalPrecipitation':
        monmean_exp1_cube.units = 'mm day'
        monmean_exp2_cube.units = 'mm day'
     

  
    anom_data = monmean_exp2_cube.data - monmean_exp1_cube.data
    anom_cube = monmean_exp1_cube.copy(data=anom_data)
    print(anom_cube)
  
    djf_anom_cube = (anom_cube[0,:,:] + anom_cube[1,:,:] + anom_cube[11,:,:]) / 3.0
    mam_anom_cube = (anom_cube[2,:,:] + anom_cube[3,:,:] + anom_cube[4,:,:]) / 3.0
    jja_anom_cube = (anom_cube[5,:,:] + anom_cube[6,:,:] + anom_cube[7,:,:]) / 3.0
    son_anom_cube = (anom_cube[8,:,:] + anom_cube[9,:,:] + anom_cube[10,:,:]) / 3.0


    djf_e1_cube = (monmean_exp1_cube[0,:,:] + monmean_exp1_cube[1,:,:] + monmean_exp1_cube[11,:,:]) / 3.0
    mam_e1_cube = (monmean_exp1_cube[2,:,:] + monmean_exp1_cube[3,:,:] + monmean_exp1_cube[4,:,:]) / 3.0
    jja_e1_cube = (monmean_exp1_cube[5,:,:] + monmean_exp1_cube[6,:,:] + monmean_exp1_cube[7,:,:]) / 3.0
    son_e1_cube = (monmean_exp1_cube[8,:,:] + monmean_exp1_cube[9,:,:] + monmean_exp1_cube[10,:,:]) / 3.0

    if FIELDNAME == 'NearSurfaceTemperature':
        djf_e1_cube.convert_units('degC')
        mam_e1_cube.convert_units('degC')
        jja_e1_cube.convert_units('degC')
        son_e1_cube.convert_units('degC')
  
    fig = plt.figure(figsize=[7,10])

    # plot them all
    plotfig(fig,[4,1],djf_e1_cube, EXPT1 + ' djf mean','n')
    plotfig(fig,[4,2],djf_anom_cube, EXPT2 + '-' + EXPT1 + ' djf mean','y')
    plotfig(fig,[4,3],mam_e1_cube, EXPT1 + ' mam mean','n')
    plotfig(fig,[4,4],mam_anom_cube, EXPT2 + '-' + EXPT1 + ' mam mean','y')
    plotfig(fig,[4,5],jja_e1_cube, EXPT1 + ' jja mean','n')
    plotfig(fig,[4,6],jja_anom_cube, EXPT2 + '-' + EXPT1 + ' jja mean','y')
    plotfig(fig,[4,7],son_e1_cube, EXPT1 + ' son mean','n')
    plotfig(fig,[4,8],son_anom_cube, EXPT2 + '-' + EXPT1 + ' son mean','y')
   
    fileout = ('/nfs/hera1/earjcti/HadCM3_plots/' + FIELDNAME + '/' + EXPT2 + '-' + EXPT1 + '_seasmean')
    plt.tight_layout()
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()




#=====================================================================
def diff_two_anomalies():
    """
    plots the difference in seasonal mean between experiment1 and experiment2
    """

    expt1e = EXPT1[0:5]
    expt1c = EXPT1[6:11]
    expt2e = EXPT2[0:5]
    expt2c = EXPT2[6:11]
  
    monmean_exp1e_cube = iris.load_cube(FILESTART + expt1e + FILEEND)
    monmean_exp1c_cube = iris.load_cube(FILESTART + expt1c + FILEEND)
    monmean_exp2e_cube = iris.load_cube(FILESTART + expt2e + FILEEND)
    monmean_exp2c_cube = iris.load_cube(FILESTART + expt2c + FILEEND)

    anom_data = ((monmean_exp2e_cube.data - monmean_exp2c_cube.data) - 
                 (monmean_exp1e_cube.data - monmean_exp1c_cube.data))

    anom_cube = monmean_exp1e_cube.copy(data=anom_data)

    djf_anom_cube = (anom_cube[0,:,:] + anom_cube[1,:,:] + anom_cube[11,:,:]) / 3.0
    mam_anom_cube = (anom_cube[2,:,:] + anom_cube[3,:,:] + anom_cube[4,:,:]) / 3.0
    jja_anom_cube = (anom_cube[5,:,:] + anom_cube[6,:,:] + anom_cube[7,:,:]) / 3.0
    son_anom_cube = (anom_cube[8,:,:] + anom_cube[9,:,:] + anom_cube[10,:,:]) / 3.0

    
    fig = plt.figure(figsize=[7,5])

    # plot them all
    plotfig(fig,[2,1],djf_anom_cube, EXPT2 + '-' + EXPT1 + ' djf','y')
    plotfig(fig,[2,2],mam_anom_cube, ' mam anom','y')
    plotfig(fig,[2,3],jja_anom_cube, ' jja anom','y')
    plotfig(fig,[2,4],son_anom_cube, ' son anom','y')
    
    fileout = ('/nfs/hera1/earjcti/HadCM3_plots/' + FIELDNAME + '/' + EXPT2 + '-' + EXPT1 + '_seasmean')
    plt.tight_layout()
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()


#
##########################################################
# main program

EXPT1 = 'xozzb'
EXPT2 = 'xpkmb'
#type = 'PliominPi' # type is PliominPi, PiminPi, PliominPlio
type = 'PliominPlio'

FIELDNAME = 'NearSurfaceTemperature'

FILESTART = '/nfs/hera1/earjcti/um/'
FILEEND = '/' + FIELDNAME + '/means/mean_month.nc'


if FIELDNAME == 'TotalPrecipitation':
    if type == 'PliominPi':
        V = np.arange(-3.0,3.5, 0.5)
    else:
        V = np.arange(-3.0, 3.5, 0.5)

if FIELDNAME == 'NearSurfaceTemperature':
    if type == 'PliominPi':
        V = np.arange(-10.0,11.0, 1.0)
    else:
        V = np.arange(-3.0, 3.5, 0.5)


if len(EXPT1) > 10 and len(EXPT2) > 10:
    diff_two_anomalies()
else:
    diff_two_experiments()

