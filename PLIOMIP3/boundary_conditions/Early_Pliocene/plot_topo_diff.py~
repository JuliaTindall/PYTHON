#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will plot the topography in the CAS region.
This is for inclusion in Alans paper.  

"""

import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import netCDF4
from datetime import date
import cartopy.crs as ccrs
import sys
import matplotlib.colors




def get_files(fileend, eoi400_fieldname, ep_fieldname):
    """
    gets the original and new pliocene
    """

    filestart = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'     
    eoi400_filename = 'Plio_enh/Plio_enh/Plio_enh_' + fileend + '.nc'
    eoi400_cube = iris.load_cube(filestart + eoi400_filename,eoi400_fieldname)

    ep_cube = iris.load_cube('EP_topo_v1.0.nc',ep_fieldname)
    return eoi400_cube, ep_cube


###############################################
def plot_difference(cube_ep, cubeEOI400):
    """
    here we plot the topography in the ep experiment
    """

   
    
    # get pliocene land sea mask
    eoi400lsm = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    eoi400_lsmcube=iris.load_cube(eoi400lsm)
  
    e280lsm = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc'
    e280_lsmcube=iris.load_cube(e280lsm)
  
    # set up colormap
    #boundaries = np.arange(-1000., 1100., 100)
    boundaries=[-3000,-2000,-1000,-500,-300,-200,-100,0,
                 100,200,300,500,1000,2000,3000]
    cmap_terrain = plt.cm.get_cmap('terrain',len(boundaries))
    colors=list(cmap_terrain(np.arange(len(boundaries))))
    cmap_terrain.set_over(colors[-1])  # set to last color
    cmap_terrain.set_under(colors[0])  # set to first color

    # save it to a file.  We are going to try and plot it with a good 
    # colorscale on Jasmin
    iris.save([cubeEOI400, cube_ep],'to_plot.nc')
    
    # plot
    reg_EP_cube = cubeEOI400.extract(iris.Constraint(latitude = lambda cell: 0< cell < 20, longitude=lambda cell: -100 < cell < -70))
    print(reg_EP_cube)
   
    ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
    cs=iplt.pcolormesh(reg_EP_cube, cmap=cmap_terrain,
        norm = matplotlib.colors.BoundaryNorm(boundaries, 
                                              ncolors=len(boundaries)-1, 
                                              clip=False))
    cbar=plt.colorbar(cs,orientation='vertical', extend='min')
    cbar.set_ticks(boundaries)
    #qplt.contour(e280_lsmcube,[0.5],colors='red',linewidths=0.75)
    #qplt.contour(eoi400_lsmcube,[0.5],colors='black',linewidths=0.75,linestyle='dashed')
    ax.gridlines()
    ax.coastlines()
    plt.show()
    sys.exit(0)
   
   # plt.gca().coastlines()
    plt.title('LP_MIN_LSM - EOI400')

    ax=plt.subplot(1,2,2,projection=ccrs.PlateCarree())
    cs=iplt.pcolormesh(cubediff, cmap='winter',vmin=-200,vmax=0,extend='both')
    plt.colorbar(cs,orientation='horizontal')
    #qplt.contour(e280_lsmcube,[0.5],colors='red',linewidths=0.75)
    #qplt.contour(eoi400_lsmcube,[0.5],colors='black',linewidths=0.75,linestyle='dashed')
    ax.gridlines()
   # plt.gca().coastlines()
    plt.title('LP_MIN_LSM - EOI400')


    plt.tight_layout()
    plt.show()
    
    sys.exit(0)
    plt.savefig('LSM_anom.png')
   
     
###############################################
# main program


# get Pliocene and preindustrial files

eoi400_topo_cube, EP_cube = get_files('topo_v1.0','p4_topo',
                                             'p4_topo_Open_CAS')

# plot the diagnostics for the Early Pliocene
plot_difference(EP_cube, eoi400_topo_cube)

