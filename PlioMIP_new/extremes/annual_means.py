#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia
Used for comparing with diags 6-9
Plots the difference between the pliocene and the preindustrial mean temperature for each month
"""
import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import sys




   
#########################################################
def plot_anomalies(cube, field_ind, plotname, title):
    """
    this is for plotting all the anomalies
    """ 

    cube.units = 'degC'
    if field_ind == 't':
        vals = np.arange(0,14,2)
        cs = qplt.contourf(cube, levels=vals, cmap='Reds', extend='both')
    if field_ind == 'p':
        vals = np.arange(-1.5,1.5,0.2)
        cs = qplt.contourf(cube, levels=vals, cmap='BrBG', extend='both')
  
    plt.gca().coastlines()
    plt.title(title) 
    
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/means/' + plotname)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()


def plot_textremes_toplevel():
    """
    for each month plots the difference between the Tmean for the pliocene and the preindustrial
    """

    # temperature
    filestart = '/nfs/hera1/earjcti/regridded/HadCM3/'
    plio_cube_nsat = iris.load_cube(filestart + 'EOI400.NearSurfaceTemperature.allmean.nc')
    pi_cube_nsat = iris.load_cube(filestart + 'E280.NearSurfaceTemperature.allmean.nc')

    print(plio_cube_nsat.data)
    print(pi_cube_nsat.data)
    
    anom_cube_nsat = plio_cube_nsat - pi_cube_nsat
    print(anom_cube_nsat.data)
    plotname = 'HadCM3_pliomip2_NSAT_annmean_' 

    plot_anomalies(anom_cube_nsat,'t', plotname,'Plio_PI: SAT')
  
    # precipitation
    filestart = '/nfs/hera1/earjcti/regridded/HadCM3/'
    plio_cube_prec = iris.load_cube(filestart + 'EOI400.TotalPrecipitation.allmean.nc')
    pi_cube_prec = iris.load_cube(filestart + 'E280.TotalPrecipitation.allmean.nc')

    
    anom_cube_prec = plio_cube_prec - pi_cube_prec
    plotname = 'HadCM3_pliomip2_prec_annmean_' 

    plot_anomalies(anom_cube_prec,'p', plotname,'Plio_PI: prec')
  
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'

EXPTNAME = 'tenvj'
CNTLNAME = 'xozza'
plot_textremes_toplevel()
 

