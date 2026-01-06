#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2023 by Julia

This program will look at the pliomip2 standard and enhanced 
files and compare them

"""

import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import netCDF4
import sys



###############################################
# main program


# get Pliocene  lsm

lsm_std = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_std/Plio_std/Plio_std_LSM_v1.0.nc'

lsmstd_cube = iris.load_cube(lsm_std,'etopo1_lsm')


lsm_enh = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'

lsmenh_cube = iris.load_cube(lsm_enh,'p4_lsm')


lsm_diff_cube = lsmenh_cube - lsmstd_cube
lsm_diff_cube.long_name = 'pliomip2 enh - standard'


# plot new land sea mask
plt.subplot(1,1,1)
qplt.contourf(lsm_diff_cube, cmap='bwr',levels=(-1.5,-0.5,0.5,1.5))
plt.gca().coastlines()
plt.show()
