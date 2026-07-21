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
from matplotlib import gridspec
import matplotlib.path as mpath
import matplotlib as mpl
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import sys



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

def get_lsm(expt):
    """
    land sea mask is where the point is ocean in both pliocene and pi
    """
    lsm_plio_file = (filestart + '/um/' + expt + '/database_averages/' + expt + 
                     '_Annual_Average_#pf_SST_3900_4000.nc')
  
    sst_plio_cube = iris.load_cube(lsm_plio_file)
    
    lsm_plio_cube = sst_plio_cube.copy(data=sst_plio_cube.data.mask)
  
    lsm_plio_cube = iris.util.squeeze(lsm_plio_cube)

    
    return lsm_plio_cube


def plot_seaice(hemisphere,cube,month,lsm_cube,exptno):
    """
    plots the sea ice for the experiment for either the NH or the southern
    hemisphere
    """
    outstart = filestart + 'um/' + expt + '/avgplots/seaice/'
    
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
         ax=fig.add_subplot(spec[exptno],projection=ccrs.NorthPolarStereo())
    if hemisphere == 'SH':   
        # ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
         ax=fig.add_subplot(spec[exptno],projection=ccrs.SouthPolarStereo())
  
    polarCentral_set_latlim(lat_lims,ax)
    axplot=iplt.pcolormesh(cube,cmap=custom_cmap,
                 norm=mpl.colors.BoundaryNorm(levels, 
                                             ncolors=len(levels)-1,
                                             clip=False))

    plt.title(name.get(expt),fontsize=14)
    iplt.contour(lsm_cube, levels=[-2,0.5,2],
                  colors='black', linewidths=0.1)

    if exptno == 2:
        cax = fig.add_subplot(spec[2:4])
        pos = cax.get_position()
        cax.set_position([pos.x0,pos.y0+0.05,pos.width,pos.height*0.1])
        cbar = fig.colorbar(axplot, cax=cax,orientation='horizontal')
        cbar.ax.tick_params(labelsize=14)
        cbar.set_label('Ice Fraction',fontsize=14)
        cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
   
        # add title
        fig.text(
            0.5, 0.98, # x, y in figure coordinates
            month + ' Sea Ice',
            ha='center',
            va='top',
            fontsize=16
        )

    #cbar = plt.colorbar(axplot,  orientation= 'vertical')
    
   
    
   


def calc_SH_cice_area(title,cube):
    """
    calculates the area of SH sea ice in 'cube'
    """

    # extract SH only and replace land with NAN
    SH_constraint = iris.Constraint(latitude = lambda cell: cell  < 0)
    SH_cube = cube.extract(SH_constraint)
    
    # get area weights (if land put grid areas to zero)
    SH_cube.coord('longitude').guess_bounds()
    SH_cube.coord('latitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(SH_cube)
    for j in range(0,72):
        for i in range(0,288):
            if SH_cube.data[j,i] < 0.0:
                SH_cube.data[j,i]=0.0

    # calculate total sea ice area and print out
    seaicearea = SH_cube.collapsed(['latitude','longitude'],
                                   iris.analysis.SUM,weights=grid_areas)
    # divide by 1E6 to convert from m2 to km2 and another 1E6 to convert to
    # million km2
    print(title,seaicearea.data / 1E12)
 
    return seaicearea.data / 1E12

#################
# MAIN PROGRAM
################

filestart = '/home/earjcti/'

name = {'xqbwc':'PI',
        'xqbwd':'LP', 'xqbwe':'EP$_{400}$',
        'xqbwr':'Late Pliocene (min_LSM)', 'xqbwg':'EP',
        'xpsid':'LP',
        'xpsie':'EP - 400ppmv','xpsig':'EP'}

expts=['xqbwd','xqbwg','xqbwe']


# setup for figure
fig = plt.figure(figsize=[10.5,5.0])

spec = gridspec.GridSpec(ncols=3, nrows=2,hspace=0.1,
                         wspace=0.1, width_ratios=[1,1,1],
                         height_ratios=[1,0.3])


for exptno,expt in enumerate(expts):
    MP_cube = iris.load_cube(filestart + 'um/'+expt+'/database_averages/'+expt+'_Monthly_Average_#pf_SeaIceConc_3900_4000.nc', 'SeaIceConc')

    lsm_plio_cube = get_lsm('xqbwd')
    
    # set land to -100
    for i in range(0,12):
        MP_cube.data[i,0,:,:] = np.where(lsm_plio_cube.data == 1.0, -100.,
                                     MP_cube.data[i,0,:,:])


    # day 45.5 represents february, 255.5 represents july
    #print(MP_cube)
    #print(MP_cube.coord('time').points)
    #print(MP_cube.extract(iris.Constraint(time=45.)))
    febMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=45.)))
    sepMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=255.)))
    annMP_cube = iris.util.squeeze(MP_cube.collapsed('time',iris.analysis.MEAN))
  

    # plot NH sea ice for all months
    #plot_seaice('NH',febMP_cube,'February',lsm_plio_cube)
    #plot_seaice('NH',sepMP_cube,'September',lsm_plio_cube)
    #plot_seaice('NH',annMP_cube,'Annual',lsm_plio_cube)

    # plot SH sea ice for all months
    #plot_seaice('SH',febMP_cube,'February',lsm_plio_cube)
    plot_seaice('SH',sepMP_cube,'September',lsm_plio_cube,exptno)
    #plot_seaice('SH',annMP_cube,'Annual',lsm_plio_cube)

plt.savefig('sea_ice_area_supp.png')
    
plt.close()

# calculate area of SH sea ice
#ciceareamp = calc_SH_cice_area('Feb ' + expt,febMP_cube)
#ciceareapi = calc_SH_cice_area('Feb ' + cntl,febPI_cube)
#print('anomaly feb',ciceareamp - ciceareapi)
#print(' ')
    
#ciceareamp=calc_SH_cice_area('Sep ' + expt,sepMP_cube)
#ciceareapi=calc_SH_cice_area('Sep ' + cntl,sepPI_cube)
#print('anomaly sep',ciceareamp - ciceareapi)
#print(' ')
    

#ciceareamp=calc_SH_cice_area('Ann ' + expt,annMP_cube)
#ciceareapi=calc_SH_cice_area('Ann ' + cntl,annPI_cube)
#print('anomaly ann',ciceareamp - ciceareapi)
#print(' ')
