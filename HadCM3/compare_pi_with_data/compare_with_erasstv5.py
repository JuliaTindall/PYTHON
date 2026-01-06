1#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.09.2019 by Julia
This program will plot the difference between the model and the NOAA ERSSTv5
The run we want to use is a preindustrial run
The run should have been preprocessed by CEMAC/PLIOMIP2/regrid_HCM3_50_year_avg.py
"""

import os 
import sys
import numpy as np
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import iris.analysis.cartography
import iris.coord_categorisation



def diffmonmean():
    """
    difference between ERSSTv5 and the experiment for each month mean
    """
    monmean_exp_cube = iris.load_cube(FILESTART + EXPT + 
                                       '/SST/means/mean_month.nc')
    monmean_ERSST_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/NOAAERSST5/E280.SST.mean_month.nc')

    # regrid the experiment cube onto the same grid as ersst
    monmean_expr_cube = monmean_exp_cube.regrid(monmean_ERSST_cube,iris.analysis.Linear())

    # find the difference
    monmean_expr_cube.attributes = None
    monmean_ERSST_cube.attributes = None
    monmean_expr_cube.remove_coord('year')
    monmean_ERSST_cube.remove_coord('year')
  
    print(monmean_expr_cube)
    print(monmean_ERSST_cube)
    diffcube = monmean_expr_cube - monmean_ERSST_cube
   
    monthnames=['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nor','dec']
    fig = plt.figure(figsize=[15.0,15.0])
    for mon in range(0,12):
        fig.add_subplot(4,3,mon+1)
        V = np.arange(-5,5,1)
        qplt.contourf(diffcube[mon,:,:],levels=V,extend='both',cmap='RdBu_r')
        plt.gca().coastlines()
   
    
        squares = np.square(diffcube.data[mon,:,:])
        avgsquares = np.mean(squares)
        rmse = np.sqrt(avgsquares)
        mae = np.mean(diffcube.data[mon,:,:])

        plt.title(monthnames[mon]+ ': RMSE:' + str(np.around(rmse,1))+ 'mae:' + str(np.around(mae,1))  + 'degC- unweighted')
        plt.savefig('/nfs/hera1/earjcti/HadCM3_plots/SST/' + EXPT + '-ERSSTv5_monmean.eps')
        plt.savefig('/nfs/hera1/earjcti/HadCM3_plots/SST/' + EXPT + '-ERSSTv5_monmean.png')

def diffannmean():
    """
    difference between ERSSTv5 and the experiment in the annual mean
    """
    annmean_exp_cube = iris.load_cube(FILESTART + EXPT + 
                                       '/SST/means/allmean.nc')
    annmean_ERSST_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/NOAAERSST5/E280.SST.allmean.nc')

    # regrid the experiment cube onto the same grid as ersst
    annmean_expr_cube = annmean_exp_cube.regrid(annmean_ERSST_cube,iris.analysis.Linear())

    # find the difference
    diffcube = annmean_expr_cube - annmean_ERSST_cube
    
    V = np.arange(-5,5,1)
    qplt.contourf(diffcube,levels=V,extend='both',cmap='RdBu_r')
    plt.gca().coastlines()
    
    squares = np.square(diffcube.data)
    avgsquares = np.mean(squares)
    rmse = np.sqrt(avgsquares)
    mae = np.mean(diffcube.data)


    plt.title(EXPT + ' - ERSST, RMSE:' + str(np.around(rmse,2))+ ', mae=' + str(np.around(mae,1)) + 'degC - unweighted')
    plt.savefig('/nfs/hera1/earjcti/HadCM3_plots/SST/' + EXPT + '-ERSSTv5_annmean.eps')
    plt.savefig('/nfs/hera1/earjcti/HadCM3_plots/SST/' + EXPT + '-ERSSTv5_annmean.png')


##########################################################
# main program

EXPT = 'xpkma'
FILESTART = '/nfs/hera1/earjcti/um/'

diffannmean()
diffmonmean()



