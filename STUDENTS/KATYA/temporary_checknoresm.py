#!/usr/bin/env python2.7
#NAME
#   quick program to check whether I have regridded Noresm correctly for katya


import os
import numpy as np
#import matplotlib as mp
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from netCDF4 import Dataset, MFDataset
import iris
import iris.analysis.cartography
import iris.quickplot as qplt
#import cartopy.crs as ccrs
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
#import subprocess




################################
def check_julia_processed_data():

  #preind_file = '/nfs/hera1/earjcti/regridded/NorESM-L/timeseries/E280.SST.0-200_#timeseries_Katya_no_ann_cycle.nc'
#plio_file = '/nfs/hera1/earjcti/regridded/NorESM-L/timeseries/EOI400.SST.0-200_#timeseries_Katya_no_ann_cycle.nc'

    preind_file = '/nfs/hera1/earjcti/regridded/NorESM-L/timeseries/E280.SST.timeseries_no_ann_cycle.nc'
    plio_file = '/nfs/hera1/earjcti/regridded/NorESM-L/timeseries/EOI400.SST.timeseries_no_ann_cycle.nc'
    

#preind_file = '/nfs/hera1/earjcti/regridded/IPSLCM6A/timeseries/E280.SST.2850-3049_timeseries_Katya_no_ann_cycle.nc'
#plio_file = '/nfs/hera1/earjcti/regridded/IPSLCM6A/timeseries/EOI400.SST.1850-2049_timeseries_Katya_no_ann_cycle.nc'
    
    preind_cube = iris.load_cube(preind_file)
    plio_cube = iris.load_cube(plio_file)
    
#=========================================================
# plot global mean to look for a drift
    preind_cube.coord('longitude').guess_bounds()
    preind_cube.coord('latitude').guess_bounds()
    preind_gridareas = iris.analysis.cartography.area_weights(preind_cube)#
    
    preind_mean_cube = preind_cube.collapsed(['latitude','longitude'],
                                             iris.analysis.MEAN, 
                                             weights=preind_gridareas)
    plio_mean_cube = plio_cube.collapsed(['latitude','longitude'],
                                         iris.analysis.MEAN, 
                                         weights=preind_gridareas)
    
#    plt.subplot(2,1,1)
#    plt.plot(preind_mean_cube.data)
#    plt.subplot(2,1,2)
#    plt.plot(plio_mean_cube.data)
#    plt.show()
#    sys.exit(0)
    
    ixmax = plio_mean_cube.data.argmax() 
    ixmin = plio_mean_cube.data.argmin() 
    
    print('noresm plio',ixmax,ixmin)
    
    pliomaxcube = iris.util.squeeze(plio_cube[ixmax,:,:,:])
    pliomincube = iris.util.squeeze(plio_cube[ixmin,:,:,:])
    
    ixmax = preind_mean_cube.data.argmax() 
    ixmin = preind_mean_cube.data.argmin() 
    
    print('noresm preind',ixmax,ixmin)
    
    
    preindmaxcube = iris.util.squeeze(preind_cube[ixmax,:,:,:])
    preindmincube = iris.util.squeeze(preind_cube[ixmin,:,:,:])
    
    figure = plt.figure(figsize=[9.0,9.0])
    
    vals = np.arange(-2.0, 2.2, 0.2)
    plt.subplot(2,3,1)
    qplt.contourf(pliomaxcube,levels=vals,cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('warmest pliocene')
    
    plt.subplot(2,3,2)
    qplt.contourf(pliomincube,levels=vals,cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('coolest pliocene')
    
    plt.subplot(2,3,3)
    qplt.contourf(pliomaxcube - pliomincube,levels=vals,cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('warmest - coolest pliocene')

    plt.subplot(2,3,4)
    qplt.contourf(preindmaxcube,levels=vals,cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('warmest pi')
    
    plt.subplot(2,3,5)
    qplt.contourf(preindmincube,levels=vals,cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('coolest pi')
    
    plt.subplot(2,3,6)
    qplt.contourf(preindmaxcube - preindmincube,levels=vals,cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('warmest-coolest pi')
#    plt.show()
    
    plt.savefig('noresm_warmest_and_coldest_arthurs_paper.png')
    

#============================================================
# extract a temperature in the NINO3.4 region
# use longitude =210, latitude = 0.5
#latcube = preind_cube.extract(iris.Constraint(latitude=0.5))
#preind_pointcube = latcube.extract(iris.Constraint(longitude=210.))


#latcube = plio_cube.extract(iris.Constraint(latitude=0.5))
#plio_pointcube = latcube.extract(iris.Constraint(longitude=210.))


#plt.subplot(2,1,1)
#plt.plot(preind_pointcube.data)
#plt.subplot(2,1,2)
#plt.plot(plio_pointcube.data)
#plt.show()


def check_raw_data():
    """
    have a look and see what the raw data is like.
    I want to see whether I have regridded it correctly.  
    """



####

check_julia_processed_data()
sys.exit(0)
