#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created 16 May 2023
#
#@author: earjcti
#"""
#
#  plot the SH (antarctic) salinity averaged over the year
#
#

import cartopy.crs as ccrs
import cartopy as cy
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib as mpl
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import iris.plot as iplt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap



def get_mean_sal(expt,startyear,endyear):
    """
    gets the mean salinity southwards of 60S.
    returns a cube of shape depths,latitudes,longitudes

    also get depths from vertical velocity
    """

    filestart = '/home/earjcti/um/' + expt + '/pg/' + expt + 'o#pg00000'

    allSH_cubes = CubeList([])
    for year in range(startyear,endyear):
        filename = filestart + str(year) + 'c1+.nc'
        cube = iris.load_cube(filename,'SALINITY (OCEAN)       (PSU-35)/1000')

        lat_constraint = iris.Constraint(latitude=lambda lat: lat < -59.5)
        SH_cube = cube.extract(lat_constraint)
        allSH_cubes.append(SH_cube)

    iris.util.equalise_attributes(allSH_cubes)
    SH_cubes = allSH_cubes.concatenate_cube()

    avg_SH_cube = SH_cubes.collapsed('t',iris.analysis.MEAN)
    # convert to psu
    avg_SH_cube.data = (avg_SH_cube.data * 1000.) + 35.

    # get depths
    wcube = iris.load_cube(filename,'VERT.VEL. ON OCEAN HALF LEVELS  CM/S')
    depths=wcube.coord('depth').points
    dz=np.empty_like(depths)
    dz[0]=depths[0]
    for k in range(1,len(depths)):
        dz[k]=depths[k]-depths[k-1]
    
    return avg_SH_cube,dz
    
def polarCentral_set_latlim(lat_lims, ax):
    ax.set_extent([-180, 180, lat_lims[0], lat_lims[1]], ccrs.PlateCarree())
    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)

    ax.set_boundary(circle, transform=ax.transAxes)


def plot_SH_salin(cube,month):
    """
    plots the SH salinity for the experiment
    """
    outstart = ('/nfs/hera1/earjcti/um/' + expt + '/avgplots/SH_salin/' + 
                str(startyear) + '_' + str(endyear) + '_')
    
    levels=np.arange(31.0, 37.00, 1.00)
    levels = [31.0,33.0,34.0,34.2,34.4,34.6,34.8,35.0,36.0,40.0]
    lat_lims=[-60,-90]

    custom_cmap = plt.cm.get_cmap('Reds',len(levels))
    custom_cmap.set_under('grey')
    custom_cmap.set_over('white')

   
    # plot
    ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    polarCentral_set_latlim(lat_lims,ax)
    axplot=iplt.pcolormesh(cube,cmap=custom_cmap,
                 norm=mpl.colors.BoundaryNorm(levels, 
                                             ncolors=len(levels)-1,
                                             clip=False))
    cbar = plt.colorbar(axplot,  orientation= 'vertical')
    cbar.set_label('psu')   
    #cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
    plt.title(month+': ' + name.get(expt))

    gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
    gl.ylocator = plt.FixedLocator(range(-90, -59, 10))  # every 10deg lat
    gl.top_labels = False
    gl.right_labels = False


    plt.show()
    plt.close()
    sys.exit(0)
    plt.savefig(outstart + expt + '_' + month + '_SH_salin.eps')

 


def plot_SH_salin_anom(cubeexpt,cubecntl,month):
    """
    plots the salinity for the SH 
    """
    outstart = ('/nfs/hera1/earjcti/um/' + expt + '/avgplots/SH_salin/'
                + str(startyear) + '-' +str(endyear))
    
    #levels=np.arange(0.0, 1.01, 0.01)
    levels=[-1.0,-0.4,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.4,1.0]
    #levels=levels/10.
    lat_lims=[-60,-90]

    custom_cmap = plt.cm.get_cmap('RdBu_r',len(levels)+1)
    custom_cmap.set_under('grey')
    custom_cmap.set_over('grey')

    cube = cubeexpt-cubecntl
    # plot
    ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    polarCentral_set_latlim(lat_lims,ax)
    axplot=iplt.pcolormesh(cube,cmap=custom_cmap,
                 norm=mpl.colors.BoundaryNorm(levels, 
                                             ncolors=len(levels)+1,
                                             clip=False))
    cbar = plt.colorbar(axplot,  orientation= 'vertical')
    cbar.set_label('psu')   
    #cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
    cbar.set_ticks(levels)  
    plt.title(month+': ' + name.get(expt) +  '-' + name.get(cntl))

    gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray')
    gl.ylocator = plt.FixedLocator(range(-90, -59, 10))  # every 10deg lat
    gl.top_labels = False
    gl.right_labels = False


    plt.show()
    plt.close()
    sys.exit(0)
    plt.savefig(outstart + expt + '-' +  cntl + '_' + month + '_SH_salin.png')
    plt.savefig(outstart + expt + '-' + cntl + '_' + month + '_SH_salin.eps')



#################
# MAIN PROGRAM
################

name = {'xqbwc':'PI','xpsic':'PI',
        'xqbwd':'Late Pliocene', 'xqbwe':'Early Pliocene - 400ppmv',
        'xqbwr':'Late Pliocene (min_LSM)', 'xqbwg':'Early Pliocene',
        'xpsig':'Early Pliocene',
        'xpsie':'Early Pliocene 400ppmv',
        'xpsid':'Late Pliocene'}

# read in multimodel mean monthly salinity data (EOI400-E280)
expt='xpsig'
cntl='xpsid'
startyear=1400
endyear=1500

# get average cube all depths
mean_sal_expt_cube,dz = get_mean_sal(expt,startyear,endyear)
mean_sal_cntl_cube,dz = get_mean_sal(cntl,startyear,endyear)

# extract top 1000m and average 
depcons=iris.Constraint(depth_1=lambda d: d <= 1000)

sal_expt_1000m_cube = mean_sal_expt_cube.extract(depcons)
sal_cntl_1000m_cube = mean_sal_cntl_cube.extract(depcons)

# get depths from vertical velocity and find weighted mean

zpoints = len(sal_expt_1000m_cube.coord('depth_1').points)
dzuse = dz[:zpoints]

sal_expt_avg_cube = sal_expt_1000m_cube.collapsed('depth_1',
                                                  iris.analysis.MEAN,
                                                  weights=dzuse)
sal_cntl_avg_cube = sal_cntl_1000m_cube.collapsed('depth_1',
                                                  iris.analysis.MEAN,
                                                  weights=dzuse)


# plot SH salinity anomaly
#plot_SH_salin(sal_expt_avg_cube,'Annual')
plot_SH_salin_anom(sal_expt_avg_cube,sal_cntl_avg_cube,'Annual')

