#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will produce the boundary conditions for the 
LP_Mod_EAIS experiment.

We will do this as follows:
1. start with late Pliocene topography.
2. Plot differences in topography between this and preindustrial topography
3. plot where one is land and the other is ocean


"""

import numpy as np
import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as mplcol
import netCDF4
from datetime import date
import sys


def get_files(fileend, pi_fieldname, plio_fieldname):
    """
    gets the pliocene and the preindustrial files
    """

    filestart = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'     
    pi_filename = 'Modern_std/Modern_std/Modern_std_' + fileend + '.nc'
    plio_filename = 'Plio_enh/Plio_enh/Plio_enh_' + fileend + '.nc'

    pi_cube = iris.load_cube(filestart + pi_filename,pi_fieldname)
    plio_cube = iris.load_cube(filestart + plio_filename,plio_fieldname)


    newcube = iris.load_cube('LP_ModEAIS_topo_v1.0.nc','p4_topo_Mod_EAIS')

    return plio_cube, pi_cube, newcube

def contourf_new_topo(cube, name):
    """
    plots the topography from the new cube
    """

    V=np.arange(0,4000,500)


    mycmap=matplotlib.cm.get_cmap('viridis')
    #mycmap=matplotlib.cm.get_cmap('cubehelix')
    #mycmap=matplotlib.cm.get_cmap('winter')
    #mycmap.set_under('blue')
    colors=list(mycmap(np.arange(len(V)+1)))
    mycmap.set_under([173.0/256.0, 216.0 / 256.0, 230.0/256.0])
    print(colors)
    print(np.shape(colors))
    #mycmap.set_over(colors[7])

    ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = qplt.contourf(cube, levels=V, cmap=mycmap,extend='both')
    ax.coastlines()
    ax.gridlines()
    plt.title(name)

    plt.savefig(name+ '.eps')


def colormesh_new_topo(cube, name):
    """
    plots the topography from the new cube
    """

    boundaries=np.arange(0,4500,500)


    mycmap=matplotlib.cm.get_cmap('viridis', len(boundaries))
    #mycmap=matplotlib.cm.get_cmap('cubehelix')
    #mycmap=matplotlib.cm.get_cmap('winter')
    #mycmap.set_under('blue')
    colors=list(mycmap(np.arange(len(boundaries))))
    mycmap.set_under([173.0/256.0, 216.0 / 256.0, 230.0/256.0])
    print(colors)
    print(np.shape(colors))
    mycmap.set_over(colors[8])

    ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(cube, cmap=mycmap,
                         norm=mplcol.BoundaryNorm(boundaries,
                                                  ncolors=len(boundaries),
                                                  clip=False))
    ax.coastlines()
    ax.gridlines()
    plt.title(name)


    cbar=plt.colorbar(cs,orientation='horizontal',extend='both')
    #cbar.set_label('depth (m)',size=30)
    #cbar.ax.tick_params(labelsize=30)
   

    plt.savefig(name+ 'colormesh.eps')



def colormesh_anom_topo(cubenew, cubeorig, name):
    """
    plots the topography from the new cube
    """

    boundaries=np.arange(-2000,2100,100)


    mycmap=matplotlib.cm.get_cmap('RdBu_r', len(boundaries))
    #mycmap=matplotlib.cm.get_cmap('cubehelix')
    #mycmap=matplotlib.cm.get_cmap('winter')
    #mycmap.set_under('blue')
    colors=list(mycmap(np.arange(len(boundaries))))
    mycmap.set_under([173.0/256.0, 216.0 / 256.0, 230.0/256.0])
    print(colors)
    print(np.shape(colors))
    mycmap.set_over(colors[8])

    ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(cubenew-cubeorig, cmap=mycmap,
                         norm=mplcol.BoundaryNorm(boundaries,
                                                  ncolors=len(boundaries),
                                                  clip=False))
    ax.coastlines()
    ax.gridlines()
    plt.title(name)


    cbar=plt.colorbar(cs,orientation='horizontal',extend='both')
    #cbar.set_label('depth (m)',size=30)
    #cbar.ax.tick_params(labelsize=30)
   

    plt.savefig(name+ 'colormesh.eps')

def colormesh_all_topo(cube1, name1, cube2, name2, cube3, name3):
    """
    plots the topography from the new cube
    """

    boundaries=np.arange(0,4500,500)

    
    mycmap=matplotlib.cm.get_cmap('viridis', len(boundaries))
    #mycmap=matplotlib.cm.get_cmap('cubehelix')
    #mycmap=matplotlib.cm.get_cmap('winter')
    #mycmap.set_under('blue')
    colors=list(mycmap(np.arange(len(boundaries))))
    mycmap.set_under([173.0/256.0, 216.0 / 256.0, 230.0/256.0])
    print(colors)
    print(np.shape(colors))
    mycmap.set_over(colors[8])

    fig = plt.figure(figsize=[9.0, 4.0])
    # plot first cube
    ax=plt.subplot(1,3,1,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(cube1, cmap=mycmap,
                         norm=mplcol.BoundaryNorm(boundaries,
                                                  ncolors=len(boundaries),
                                                  clip=False))
    ax.coastlines()
    ax.gridlines()
    plt.title(name1,size=15)

    # plot second cube
    ax=plt.subplot(1,3,2,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(cube2, cmap=mycmap,
                         norm=mplcol.BoundaryNorm(boundaries,
                                                  ncolors=len(boundaries),
                                                  clip=False))
    ax.coastlines()
    ax.gridlines()
    plt.title(name2,size=15)


    # plot third cube
    ax=plt.subplot(1,3,3,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(cube3, cmap=mycmap,
                         norm=mplcol.BoundaryNorm(boundaries,
                                                  ncolors=len(boundaries),
                                                  clip=False))
    ax.coastlines()
    ax.gridlines()
    plt.title(name3,size=15)

    # make an axis for the colorbar 
    plt.subplots_adjust(bottom=0.25)
    cax = plt.axes([0.1, 0.15, 0.8, 0.05])
    cbar=plt.colorbar(cs,cax=cax,orientation='horizontal',extend='both')
    cbar.set_label('height (m)',size=13)
    cbar.ax.tick_params(labelsize=13)
   

    plt.savefig('allplots_colormesh.eps')
    plt.savefig('allplots_colormesh.png')


     
###############################################
# main program


# get Pliocene and preindustrial topography and also the topography from 
# the new file

(pliocore_topo_cube, 
 pi_topo_cube,
 new_topo_cube) = get_files('topo_v1.0', 'etopo1_topo','p4_topo')

# plot the topography using contourf
#contourf_new_topo(new_topo_cube,'new_topography')
#contourf_new_topo(pliocore_topo_cube,'EOI400_topography')
#contourf_new_topo(pi_topo_cube,'E280_topography')


# plot the topography using colormesh
#colormesh_new_topo(new_topo_cube,'new_topography')
#colormesh_new_topo(pliocore_topo_cube,'EOI400_topography')
#colormesh_new_topo(pi_topo_cube,'E280_topography')


# plot the anomalies using colormesh
#colormesh_anom_topo(new_topo_cube,pliocore_topo_cube,'new-EOI400')
#colormesh_anom_topo(new_topo_cube,pi_topo_cube,'new-E280')
#colormesh_anom_topo(pliocore_topo_cube, pi_topo_cube,'EOI400-E280')

# plot the topographies using colormesh on one plot
colormesh_all_topo(pi_topo_cube, 'PI', pliocore_topo_cube, 'LP',
                   new_topo_cube, 'LP_pi-EAIS')
