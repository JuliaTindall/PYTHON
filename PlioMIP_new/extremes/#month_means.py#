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


def get_data(expt,extra):
    """
    reads in the maximum  and minimum temperature for the year and puts it in 
    a single cube
    """
    filestart = '/nfs/hera1/earjcti/um/' + expt + '/pd/' + expt + 'a@pd'+ extra
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']

    allmonths_cubes = iris.cube.CubeList([])

    for i, month in enumerate(months):
        meanT_cubelist = iris.cube.CubeList([])
  
        for year in range(0,80):
            filename = filestart + np.str(year).zfill(2) + month + '.nc'
            cube=iris.load_cube(filename,'TEMPERATURE AT 1.5M')
            meanT_cubelist.append(cube)
      
        iris.util.equalise_attributes(meanT_cubelist)
        allcubes = meanT_cubelist.concatenate_cube()
        meanmoncube = allcubes.collapsed('t',iris.analysis.MEAN)
        meanmoncube.coord('ht').points=i+1
       
        allmonths_cubes.append(meanmoncube)

    return allmonths_cubes 
   



   
#########################################################
def plot_anomalies(cube, monthname, ocn_mask, plotname, plottype):
    """
    this is for plotting all the anomalies
    """ 

    cube.units = 'degC'
    maskcube=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
    vals = np.arange(-12,14,2)
    # process and plot
    if ocn_mask == 'y':
        cube.data.mask = (maskcube.data - 1.0) * (-1.0)
    cs = qplt.contourf(cube, levels=vals, cmap='RdBu_r', extend='both')
    plt.gca().coastlines()
    plt.title(plottype + ' Plio-Pi ' + monthname) 
    
   
    #cb_ax = fig.add_axes([0.35, 0.05, 0.30, 0.02])
           
    #cbar = fig.colorbar(cs, cax=cb_ax, orientation='horizontal')
    #cbar.set_label('degC', fontsize=10)
         
    #cbar.ax.tick_params(labelsize=10)
  
    # save to file
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/means/' + plotname +  monthname)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()


def plot_textremes_toplevel(expt,exptextra,cntl,cntlextra):
    """
    for each month plots the difference between the Tmean for the pliocene and the preindustrial
    """
    # mean_pliocubes is a cubelist one per month
    mean_pliocubes = get_data(expt,exptextra)
    mean_picubes = get_data(cntl,cntlextra)

    jan_anom = iris.util.squeeze(mean_pliocubes[0] - mean_picubes[0])
    feb_anom = iris.util.squeeze(mean_pliocubes[1] - mean_picubes[1])
    mar_anom = iris.util.squeeze(mean_pliocubes[2] - mean_picubes[2])
    apr_anom = iris.util.squeeze(mean_pliocubes[3] - mean_picubes[3])
    may_anom = iris.util.squeeze(mean_pliocubes[4] - mean_picubes[4])
    jun_anom = iris.util.squeeze(mean_pliocubes[5] - mean_picubes[5])
    jul_anom = iris.util.squeeze(mean_pliocubes[6] - mean_picubes[6])
    aug_anom = iris.util.squeeze(mean_pliocubes[7] - mean_picubes[7])
    sep_anom = iris.util.squeeze(mean_pliocubes[8] - mean_picubes[8])
    oct_anom = iris.util.squeeze(mean_pliocubes[9] - mean_picubes[9])
    nov_anom = iris.util.squeeze(mean_pliocubes[10] - mean_picubes[10])
    dec_anom = iris.util.squeeze(mean_pliocubes[11] - mean_picubes[11])
    plotname = 'Tanom_' + EXPTNAME +  '-' + CNTLNAME + '_' 

    # plot monthly mean temperature anomaly for this month
    plot_anomalies(jan_anom,'jan','y', plotname, 'Tmean')
    plot_anomalies(feb_anom,'feb','y', plotname, 'Tmean')
    plot_anomalies(mar_anom,'mar','y', plotname, 'Tmean')
    plot_anomalies(apr_anom,'apr','y', plotname, 'Tmean')
    plot_anomalies(may_anom,'may','y', plotname, 'Tmean')
    plot_anomalies(jun_anom,'jun','y', plotname, 'Tmean')
    plot_anomalies(jul_anom,'jul','y', plotname, 'Tmean')
    plot_anomalies(aug_anom,'aug','y', plotname, 'Tmean')
    plot_anomalies(sep_anom,'sep','y', plotname, 'Tmean')
    plot_anomalies(oct_anom,'oct','y', plotname, 'Tmean')
    plot_anomalies(nov_anom,'nov','y', plotname, 'Tmean')
    plot_anomalies(dec_anom,'dec','y', plotname, 'Tmean')

    # get max Tmax for month
    filename = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag6-9/' + 
                EXPTNAME + '_' + CNTLNAME + '_maxTmax_jan.nc')
    jancube = iris.load_cube(filename, 'Plio - PI: Mean max of Tmax')
    plotname = 'maxTmax_anom_' + EXPTNAME +  '-' + CNTLNAME + '_' 

    plot_anomalies(jancube,'jan','y', plotname, 'maxTmax')

     # get max Tmax for month
    filename = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag6-9/' + 
                EXPTNAME + '_' + CNTLNAME + '_maxTmax_jul.nc')
    jancube = iris.load_cube(filename, 'Plio - PI: Mean max of Tmax')
    plotname = 'maxTmax_anom_' + EXPTNAME +  '-' + CNTLNAME + '_' 

    plot_anomalies(jancube,'jul','y', plotname, 'maxTmax')

  
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'

EXPTNAME = 'tenvj'
CNTLNAME = 'xozza'
plot_textremes_toplevel(EXPTNAME,'o',CNTLNAME,'o')
 

