#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 31.10.2022 by Julia

This program will plot E280 - NOAA-ERSSTv5 for Harry for the North atlantic region as requested
"""

import os
import sys
import numpy as np
#import matplotlib as mp
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#import netCDF4
#from netCDF4 import Dataset, MFDataset
import iris
import iris.analysis.cartography
import iris.quickplot as qplt
import iris.coord_categorisation
import cartopy
import cartopy.crs as ccrs

#os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
#os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
#from mpl_toolkits.basemap import Basemap, shiftgrid



##########################################################
# main program
# set up variable information

E280_SST_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/SST_multimodelmean.nc','SSTmean_pi')

NOAA_ERSST_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/NOAAERSST5/E280.SST.allmean.nc')

NOAA_ERSSTr_cube = NOAA_ERSST_cube.regrid(E280_SST_cube,iris.analysis.Linear())

anom_cube = E280_SST_cube - NOAA_ERSSTr_cube


ax1 = plt.subplot(1, 1, 1, projection=ccrs.PlateCarree())
ax1.set_extent([-80, 35, 20, 80], ccrs.PlateCarree())
V = np.arange(-4.0,4.2,0.2)
qplt.contourf(anom_cube,cmap='RdBu_r',levels=V,extend='both')
plt.gca().coastlines()
plt.title('E280 (MMM) - NOAA-ERSSTV5')
plt.savefig('E280_mmm_NOAA_ERSSTv5_Natl.eps')
plt.savefig('E280_mmm_NOAA_ERSSTv5_Natl.png')
plt.close()


V = np.arange(-4.0,4.2,0.2)
qplt.contourf(anom_cube,cmap='RdBu_r',levels=V,extend='both')
plt.gca().coastlines()
plt.title('E280 (MMM) - NOAA-ERSSTV5')
plt.savefig('E280_mmm_NOAA_ERSSTv5_globe.eps')
plt.savefig('E280_mmm_NOAA_ERSSTv5_globe.png')
plt.close()
