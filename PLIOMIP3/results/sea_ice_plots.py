#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created 16 May 2023
#
#@author: earjcti
#"""
#
#  plot the sea ice for summer and winter for NH and SH
#
#

import cartopy.crs as ccrs
import cartopy as cy
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib as mpl
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap



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

def get_lsm(expt,cntl):
    """
    land sea mask is where the point is ocean in both pliocene and pi
    """
    filestart = '/nfs/hera1/earjcti/um/'
    lsm_plio_file = (filestart + expt + '/database_averages/' + expt + 
                     '_Annual_Average_#pf_SST_3900_4000.nc')
    lsm_pi_file = (filestart + cntl + '/database_averages/' + cntl + 
                     '_Annual_Average_#pf_SST_3900_4000.nc')

    sst_pi_cube = iris.load_cube(lsm_pi_file,)
    sst_plio_cube = iris.load_cube(lsm_plio_file)
    
    lsm_pi_cube = sst_pi_cube.copy(data=sst_pi_cube.data.mask)
    lsm_plio_cube = sst_plio_cube.copy(data=sst_plio_cube.data.mask)
  
    lsm_pi_cube = iris.util.squeeze(lsm_pi_cube)
    lsm_plio_cube = iris.util.squeeze(lsm_plio_cube)

    # find where both cubes are land
    lsm_both_data = np.minimum(lsm_pi_cube.data,lsm_plio_cube.data)
    lsm_both_cube = lsm_plio_cube.copy(lsm_both_data)
    
    return lsm_pi_cube, lsm_plio_cube, lsm_both_cube


def plot_seaice(hemisphere,cube,month,lsm_cube):
    """
    plots the sea ice for the experiment for either the NH or the southern
    hemisphere
    """
    outstart = '/nfs/hera1/earjcti/um/' + expt + '/avgplots/seaice/'
    
    levels=np.arange(0.0, 1.01, 0.01)
    if hemisphere == 'NH':
        lat_lims=[60,90]
    if hemisphere == 'SH':
        lat_lims=[-60,-90]

    custom_cmap = plt.cm.get_cmap('Blues_r',len(levels))
    custom_cmap.set_under('grey')
    custom_cmap.set_over('white')


    # plot
    if hemisphere == 'NH':   
         ax=plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
    if hemisphere == 'SH':   
         ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    polarCentral_set_latlim(lat_lims,ax)
    axplot=iplt.pcolormesh(cube,cmap=custom_cmap,
                 norm=mpl.colors.BoundaryNorm(levels, 
                                             ncolors=len(levels)-1,
                                             clip=False))
    cbar = plt.colorbar(axplot,  orientation= 'vertical')
    cbar.set_label('Ice Fraction')   
    cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
    plt.title(month+': ' + name.get(expt))
    iplt.contour(lsm_cube, levels=[-2,0.5,2],
                  colors='black', linewidths=0.1)
    plt.savefig(outstart + expt + '_' + month + hemisphere + '_seaice.png')
    plt.savefig(outstart + expt + '_' + month + hemisphere + '_seaice.eps')

 


def plot_seaice_anom(hemisphere,cubeexpt,cubecntl,month,lsm_cube):
    """
    plots the sea ice for the experiment for either the NH or the southern
    hemisphere
    """
    outstart = '/nfs/hera1/earjcti/um/' + expt + '/avgplots/seaice/'
    
    #levels=np.arange(0.0, 1.01, 0.01)
    levels=[-1.0,-0.4,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.4,1.0]
    #levels=levels/10.
    if hemisphere == 'NH':
        lat_lims=[60,90]
    if hemisphere == 'SH':
        lat_lims=[-60,-90]

    custom_cmap = plt.cm.get_cmap('RdBu',len(levels)+1)
    custom_cmap.set_under('grey')
    custom_cmap.set_over('white')

    cube = cubeexpt-cubecntl
    # plot
    if hemisphere == 'NH':   
         ax=plt.subplot(1,1,1,projection=ccrs.NorthPolarStereo())
    if hemisphere == 'SH':   
         ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    polarCentral_set_latlim(lat_lims,ax)
    axplot=iplt.pcolormesh(cube,cmap=custom_cmap,
                 norm=mpl.colors.BoundaryNorm(levels, 
                                             ncolors=len(levels)+1,
                                             clip=False))
    cbar = plt.colorbar(axplot,  orientation= 'vertical')
    cbar.set_label('Ice Fraction')   
    #cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
    cbar.set_ticks(levels)  
    plt.title(month+': ' + name.get(expt) +  '-' + name.get(cntl))
    iplt.contour(lsm_cube, levels=[-2,0.5,2],
                  colors='black', linewidths=0.1)
    plt.savefig(outstart + expt + '-' +  cntl + '_' + month + hemisphere + '_seaice.png')
    plt.savefig(outstart + expt + '-' + cntl + '_' + month + hemisphere + '_seaice.eps')

 

#################
# MAIN PROGRAM
################

name = {'xqbwc':'PI',
        'xqbwd':'Late Pliocene', 'xqbwe':'Early Pliocene - 400ppmv',
        'xqbwr':'Late Pliocene (min_LSM)', 'xqbwg':'Early Pliocene',
        'xqbwe':'Early Pliocene 400ppmv'}

# read in multimodel mean monthly SST data (EOI400-E280)
expt='xqbwr'
cntl='xqbwd'

MP_cube = iris.load_cube('/nfs/hera1/earjcti/um/'+expt+'/database_averages/'+expt+'_Monthly_Average_#pd_SeaIceConc_3900_4000.nc', 'SeaIceConc')
PI_cube = iris.load_cube('/nfs/hera1/earjcti/um/'+cntl+'/database_averages/'+cntl+'_Monthly_Average_#pd_SeaIceConc_3900_4000.nc', 'SeaIceConc')
#anom_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/NearSurfaceTemperature_multimodelmean_month.nc',
#                           'NearSurfaceTemperatureplio - pi')

lsm_cube, lsm_plio_cube, lsm_both_cube = get_lsm(expt,cntl)
# set land to -100
for i in range(0,12):
    MP_cube.data[i,0,:,:] = np.where(lsm_plio_cube.data == 1.0, -100.,
                                     MP_cube.data[i,0,:,:])




# day 45.5 represents february, 255.5 represents july
febMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=45.5)))
sepMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=255.5)))
annMP_cube = iris.util.squeeze(MP_cube.collapsed('time',iris.analysis.MEAN))
febPI_cube = iris.util.squeeze(PI_cube.extract(iris.Constraint(time=45.5)))
sepPI_cube = iris.util.squeeze(PI_cube.extract(iris.Constraint(time=255.5)))
annPI_cube = iris.util.squeeze(PI_cube.collapsed('time',iris.analysis.MEAN))


# plot NH sea ice for all months
#plot_seaice('NH',febMP_cube,'February',lsm_plio_cube)
#plot_seaice('NH',sepMP_cube,'September',lsm_plio_cube)
#plot_seaice('NH',annMP_cube,'Annual',lsm_plio_cube)

# plot SH sea ice for all months
#plot_seaice('SH',febMP_cube,'February',lsm_plio_cube)
#plot_seaice('SH',sepMP_cube,'September',lsm_plio_cube)
#plot_seaice('SH',annMP_cube,'Annual',lsm_plio_cube)


# plot NH sea ice anomaly for all months
plot_seaice_anom('NH',febMP_cube,febPI_cube,'February',lsm_both_cube)
plot_seaice_anom('NH',sepMP_cube,sepPI_cube,'September',lsm_both_cube)
plot_seaice_anom('NH',annMP_cube,annPI_cube,'Annual',lsm_both_cube)

plot_seaice_anom('SH',febMP_cube,febPI_cube,'February',lsm_both_cube)
plot_seaice_anom('SH',sepMP_cube,sepPI_cube,'September',lsm_both_cube)
plot_seaice_anom('SH',annMP_cube,annPI_cube,'Annual',lsm_both_cube)
