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

def get_lsm(expt,cntl):
    """
    land sea mask is where the point is ocean in both pliocene and pi
    """
    lsm_plio_file = (filestart + '/um/' + expt + '/database_averages/' + expt + 
                     '_Annual_Average_#pf_SST_3900_4000.nc')
    lsm_pi_file = (filestart + '/um/' + cntl + '/database_averages/' + cntl + 
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
    outstart = filestart + 'um/' + expt + '/avgplots/seaice/depth_'
    
    levels=np.arange(0.0, 6.5, 0.5)
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
    cbar.set_label('Depth')
    tickvals = np.arange(0,100,10)
    #cbar.set_ticks(tickvals)
    #cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]*100.)  
    plt.title(month+': ' + name.get(expt))
    iplt.contour(lsm_cube, levels=[-2,0.5,2],
                  colors='black', linewidths=0.1)
    plt.savefig(outstart + expt + '_' + month + hemisphere + '2_22_seaice.png')
    plt.savefig(outstart + expt + '_' + month + hemisphere + '2_22_seaice.eps')
    plt.close()

    
   


def plot_seaice_anom(hemisphere,cubeexpt,cubecntl,month,lsm_cube):
    """
    plots the sea ice for the experiment for either the NH or the southern
    hemisphere
    """
    outstart = filestart + 'um/' + expt + '/avgplots/seaice/depth_'

    print('in anomaly',outstart)
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
    cbar.set_label('metres')   
    #cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
    cbar.set_ticks(levels)  
    plt.title(month+': ' + name.get(expt) +  '-' + name.get(cntl))
    iplt.contour(lsm_cube, levels=[-2,0.5,2],
                  colors='black', linewidths=0.1)
    plt.savefig(outstart + expt + '-' +  cntl + '_' + month + hemisphere + '2_22_seaice.png')
    plt.savefig(outstart + expt + '-' + cntl + '_' + month + hemisphere + '2_22_seaice.eps')
    plt.close()


def calc_SH_cice_volume(title,cube):
    """
    calculates the volume of SH sea ice in 'cube'
    """

    # extract SH only and replace land with NAN
    SH_constraint = iris.Constraint(latitude = lambda cell: cell  < 0)
    SH_cube = cube.extract(SH_constraint)
    
    # get volume weights (if land put grid volumes to zero)
    SH_cube.coord('longitude').guess_bounds()
    SH_cube.coord('latitude').guess_bounds()
    grid_volumes = iris.analysis.cartography.area_weights(SH_cube)
    for j in range(0,72):
        for i in range(0,288):
            if SH_cube.data[j,i] < 0.0:
                SH_cube.data[j,i]=0.0

    # calculate total sea ice volume and print out
    seaicevolume = SH_cube.collapsed(['latitude','longitude'],
                                   iris.analysis.SUM,weights=grid_volumes)
    # divide by 1E9 to convert from m3 to km3
    
    print(title,seaicevolume.data)
 
    return seaicevolume.data 


def calc_volume_ocean(time_ind,latitude,depth):
    """
    calculates the volume of the ocean southeward of latitude and to a 
    depth of depth
    """

    if time_ind == 'pliocene':
        file = '/home/earjcti/um/xpsid/pg/xpsido#pg000001899c1+.nc'
        cube = iris.load_cube(file,'POTENTIAL TEMPERATURE (OCEAN)  DEG.C')


    depth_constraint = iris.Constraint(depth_1 = lambda cell: cell  < depth)
    cube_surf = cube.extract(depth_constraint)

    lat_constraint = iris.Constraint(latitude = lambda cell: cell  < latitude)
    cube_subset = cube_surf.extract(lat_constraint)

    # calculate volume
    #1. use mask
    data_use = np.where(cube_subset.data.mask == True, 0.0, 1.0)
    cube_use = cube_subset.copy(data=data_use)
    cube_use = iris.util.squeeze(cube_use)
    print(cube_use)
    #qplt.contourf(cube_use[0,0,:,:])

    # 2. get weights
    cube_use.coord('longitude').guess_bounds()
    cube_use.coord('latitude').guess_bounds()
    cube_use.coord('depth_1').guess_bounds()
    grid_areas = np.squeeze(iris.analysis.cartography.area_weights(cube_use))
    depth_bounds = cube_use.coord('depth_1').bounds
    thickness = np.diff(depth_bounds, axis=1).squeeze()  # shape: (depth,)

    
    volume_weights = np.ones(cube_use.shape)

    volume_weights *= thickness[:, None, None]
    volume_weights *= grid_areas[0, :, :]

  
    #3: calculate volume of each level
   
    volume = cube_use.collapsed(['depth_1','latitude','longitude'],
                                   iris.analysis.SUM,weights=volume_weights)

    print(f"{volume.data:2e}")
    print('volume of ocean km3',f"{volume.data/1E9}")
    
    sys.exit(0)
    

#################
# MAIN PROGRAM
################

filestart = '/home/earjcti/'

name = {'xqbwc':'PI',
        'xqbwd':'Late Pliocene', 'xqbwe':'Early Pliocene - 400ppmv',
        'xqbwr':'Late Pliocene (min_LSM)', 'xqbwg':'Early Pliocene',
        'xqbwe':'Early Pliocene 400ppmv','xpsid':'LP',
        'xpsie':'EP - 400ppmv','xpsig':'EP'}



# read in multimodel mean monthly SST data (EOI400-E280)
expt='xpsig'
cntl='xpsid'

MP_cube = iris.load_cube(filestart + 'um/'+expt+'/database_averages/'+expt+'_Monthly_Average_#pf_SeaIceDepth_2_22.nc', 'SeaIceDepth')
PI_cube = iris.load_cube(filestart + 'um/'+cntl+'/database_averages/'+cntl+'_Monthly_Average_#pf_SeaIceDepth_2_22.nc', 'SeaIceDepth')
#anom_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/NearSurfaceTemperature_multimodelmean_month.nc',
#                           'NearSurfaceTemperatureplio - pi')



lsm_cube, lsm_plio_cube, lsm_both_cube = get_lsm('xqbwd','xqbwd')
# set land to -100
for i in range(0,12):
    MP_cube.data[i,0,:,:] = np.where(lsm_plio_cube.data == 1.0, -100.,
                                     MP_cube.data[i,0,:,:])




# day 45.5 represents february, 255.5 represents july
#print(MP_cube)
#print(MP_cube.coord('time').points)
#print(MP_cube.extract(iris.Constraint(time=45.)))
#febMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=45.5)))
#sepMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=255.5)))
febMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=45.)))
sepMP_cube = iris.util.squeeze(MP_cube.extract(iris.Constraint(time=255.)))
annMP_cube = iris.util.squeeze(MP_cube.collapsed('time',iris.analysis.MEAN))
#febPI_cube = iris.util.squeeze(PI_cube.extract(iris.Constraint(time=45.5)))
#sepPI_cube = iris.util.squeeze(PI_cube.extract(iris.Constraint(time=255.5)))
febPI_cube = iris.util.squeeze(PI_cube.extract(iris.Constraint(time=45.)))
sepPI_cube = iris.util.squeeze(PI_cube.extract(iris.Constraint(time=255.)))
annPI_cube = iris.util.squeeze(PI_cube.collapsed('time',iris.analysis.MEAN))


# plot NH sea ice for all months
#plot_seaice('NH',febMP_cube,'February',lsm_plio_cube)
#plot_seaice('NH',sepMP_cube,'September',lsm_plio_cube)
#plot_seaice('NH',annMP_cube,'Annual',lsm_plio_cube)

# plot SH sea ice for all months
plot_seaice('SH',febMP_cube,'February',lsm_plio_cube)
plot_seaice('SH',sepMP_cube,'September',lsm_plio_cube)
plot_seaice('SH',annMP_cube,'Annual',lsm_plio_cube)


# plot NH sea ice anomaly for all months
#plot_seaice_anom('NH',febMP_cube,febPI_cube,'February',lsm_both_cube)
#plot_seaice_anom('NH',sepMP_cube,sepPI_cube,'September',lsm_both_cube)
#plot_seaice_anom('NH',annMP_cube,annPI_cube,'Annual',lsm_both_cube)

plot_seaice_anom('SH',febMP_cube,febPI_cube,'February',lsm_both_cube)
plot_seaice_anom('SH',sepMP_cube,sepPI_cube,'September',lsm_both_cube)
plot_seaice_anom('SH',annMP_cube,annPI_cube,'Annual',lsm_both_cube)


# calculate volume of SH sea ice
cicevolumemp = calc_SH_cice_volume('Feb ' + expt,febMP_cube)
cicevolumepi = calc_SH_cice_volume('Feb ' + cntl,febPI_cube)
print('volume anomaly feb',cicevolumemp - cicevolumepi)
print(' ')
    
cicevolumemp=calc_SH_cice_volume('Sep ' + expt,sepMP_cube)
cicevolumepi=calc_SH_cice_volume('Sep ' + cntl,sepPI_cube)
print('volume anomaly sep',cicevolumemp - cicevolumepi)
print(' ')
    

cicevolumemp=calc_SH_cice_volume('Ann ' + expt,annMP_cube)
cicevolumepi=calc_SH_cice_volume('Ann ' + cntl,annPI_cube)
print('cicevolumemp',f'{cicevolumemp:2e}')
print('volume anomaly ann (m3)',f"{cicevolumemp - cicevolumepi:2e}")
print('volume anomaly ann (km3)',f"{(cicevolumemp - cicevolumepi)/1E9}")
print(' ')
sys.exit(0)

volume_ocean = calc_volume_ocean('pliocene',-70,100)
print('volume of ocean polewards of 70S to a depth of 100m is',volume_ocean)
sys.exit(0)
