#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This program will difference some of the ancil files
between the preindustrial and the pliocene experiments

This is to show the changes we have made between the paper.
Created on Fri Dec 20 10:33:08 2019

@author: earjcti
"""

import iris
import numpy as np
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt

import sys


def get_data(filename, fieldreq):
    """

    Parameters
    ----------
    filename : name of the file (start dump containing the data
    fieldreq : which field we want to get

    Returns
    -------
    an iris cube containing the required data

    """
    
    if fieldreq == 'orog':
        cubefull = iris.load_cube(filename, 'OROGRAPHY (/STRAT LOWER BC)')
        cube = cubefull.collapsed(['surface','t'],iris.analysis.MEAN)
   
    if fieldreq == 'ice' and MODEL == 'HadCM3':
        cubefull = iris.load(filename, 'FRACTIONS OF SURFACE TYPES')
        cubeice = cubefull[4]
        
        cube = cubeice.collapsed(['surface','t'],iris.analysis.MEAN)
        
    if fieldreq == 'ice' and MODEL == 'HadGEM2':
        cubefull = iris.load_cube(filename, 'FRACTIONS OF SURFACE TYPES')
        cube_ice = cubefull.extract(iris.Constraint(pseudo=9))
        cube = cube_ice.collapsed(['t'],iris.analysis.MEAN)
    
        #qplt.contourf(cube)
        #plt.show(0)
        #sys.exit(0)
    
    
    return cube

###########################################################################

FIELD = 'ice'   # field can be orog or ice
MODEL = 'HadCM3' # model can be HadGEM2 or HadCM3
OUTSTART = '/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/ancil/'



if MODEL == 'HadGEM2':
    FILENAME_PLIO = ('/nfs/hera1/earjcti/um/HadGEM_data' +  
                   '/xkvjg/netcdf/xkvjga@dan10c1.nc')
    FILENAME_PI = ('/nfs/hera1/earjcti/um/HadGEM_data' +  
                   '/xkvje/netcdf/xkvjea@dao00c1.nc')
    
if MODEL == 'HadCM3':
    FILENAME_PLIO = ('/nfs/hera1/earjcti/um/' +  
                   '/xibol/netcdf/xibola@daz00c1.nc')
    FILENAME_PI = ('/nfs/hera1/earjcti/um/' +  
                   '/xiboi/netcdf/xiboia@dat45c1.nc')

plio_cube = get_data(FILENAME_PLIO, FIELD)
pi_cube = get_data(FILENAME_PI, FIELD)

cube_anom = (plio_cube - pi_cube) 
if FIELD == 'orog':
    V = np.arange(-2000, 2100, 100)
if FIELD == 'ice':
    V = np.arange(-1.0, 1.1, 0.1)
cs = iplt.contourf(cube_anom, V, cmap='bwr', extend='both')
plt.gca().coastlines()
cbar = plt.colorbar(cs, orientation='horizontal')

if FIELD == 'orog' and MODEL == 'HadGEM2':
    plt.title('a) HadGEM2: mPWP-PI orography', loc='left')
    cbar.set_label('m')
    plt.savefig(OUTSTART + 'HadGEM2_orog_anom.pdf')
if FIELD == 'orog' and MODEL == 'HadCM3':
    plt.title('b) HadCM3: mPWP-PI orography', loc='left')
    cbar.set_label('m')
    plt.savefig(OUTSTART + 'HadCM3_orog_anom.pdf')
if FIELD == 'ice' and MODEL == 'HadGEM2':
    plt.title('c) HadGEM2: mPWP-PI ice fraction', loc='left')
    cbar.set_label('fraction')
    plt.savefig(OUTSTART + 'HadGEM2_ice_anom.pdf')
if FIELD == 'ice' and MODEL == 'HadCM3':
    plt.title('d) HadCM3: mPWP-PI ice fraction', loc='left')
    cbar.set_label('fraction')
    plt.savefig(OUTSTART + 'HadCM3_ice_anom.pdf')

plt.close()
#HadGEM2
#xkvjg pliocene (ancil_jct/HadGEM/3ma)
#                   HadGEM_Pliocene_orog
#                   HadGEM_pliocene_vegfrac
#xkvje preindustrial
#                   
#HG2 filename
#/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/xkvjga@dan10c1.nc
#'OROGRAPHY (/STRAT LOWER BC)'
#'FRACTIONS OF SURFACE TYPES' level 9
#...xkfje/netcdf/xkvjea@dao00c1