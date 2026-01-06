"""
#NAME
#    salinity_depth.py
#PURPOSE 
#
#  This program will plot the salinity at a given depth
#  it will also provide a difference from a control
"""

# Import necessary libraries

import os
import numpy as np
import math
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import iris
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.quickplot as qplt
import iris.plot as iplt
from netCDF4 import Dataset, MFDataset
import sys
#from netCDF4 import Dataset, MFDataset

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


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
    lsm_plio_file = (filestart + expt + '/pg/' + expt + 
                     'o#pg000003900c1+.nc')
    lsm_pi_file = (filestart + cntl + '/pg/' + cntl + 
                     'o#pg000003900c1+.nc')

    sst_pi_cube = iris.load_cube(lsm_pi_file,
                                 'SALINITY (OCEAN)       (PSU-35)/1000')
    sst_plio_cube = iris.load_cube(lsm_plio_file,
                                   'SALINITY (OCEAN)       (PSU-35)/1000')
    
    lsm_piall_cube = sst_pi_cube.copy(data=sst_pi_cube.data.mask)
    lsm_plioall_cube = sst_plio_cube.copy(data=sst_plio_cube.data.mask)
  
    lsm_piall_cube = iris.util.squeeze(lsm_piall_cube)
    lsm_plioall_cube = iris.util.squeeze(lsm_plioall_cube)

    lsm_pi_cube = lsm_piall_cube[LEVEL,:,:]
    lsm_plio_cube = lsm_plioall_cube[LEVEL,:,:]


    # find where both cubes are land
    lsm_both_data = np.minimum(lsm_pi_cube.data,lsm_plio_cube.data)
    lsm_both_cube = lsm_plio_cube.copy(lsm_both_data)

    
    return lsm_pi_cube, lsm_plio_cube, lsm_both_cube


def plot_SH_salin(cube,lsm_cube,expt,startyear,endyear):
    """
    plots the SH salinity for the experiment
    """
    outstart = ('/nfs/hera1/earjcti/um/' + expt + '/avgplots/salinity/' + 
                str(startyear) + '_' + str(endyear) + '_')
    
    levels=np.arange(31.0, 35.5, 0.5)
    lat_lims=[-60,-90]

    custom_cmap = plt.cm.get_cmap('viridis',len(levels))
    custom_cmap.set_under('grey')
    custom_cmap.set_over('white')


    # plot
    ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    polarCentral_set_latlim(lat_lims,ax)
    axplot=iplt.pcolormesh(cube,cmap=custom_cmap,
                 norm=mp.colors.BoundaryNorm(levels, 
                                             ncolors=len(levels)-1,
                                             clip=False))
    cbar = plt.colorbar(axplot,  orientation= 'vertical')
    cbar.set_label('psu')   
    cbar.set_ticks(np.arange(31.0,35.5,0.5))  
    plt.title(name.get(expt))

    iplt.contour(lsm_cube, levels=[-2,0.5,2],
                  colors='black', linewidths=0.1)
    print(outstart)
    plt.savefig(outstart + expt + '_'  + '_SH_salin.eps')
    plt.savefig(outstart + expt + '_'  + '_SH_salin.png')

 


def plot_SH_salin_anom(cubeexpt,cubecntl,lsm_cube,expt,cntl,startyear,endyear):
    """
    plots the salinity for the SH 
    """
    outstart = ('/nfs/hera1/earjcti/um/' + expt + '/avgplots/salinity/'
                + str(startyear) + '-' +str(endyear))
    
    #levels=np.arange(0.0, 1.01, 0.01)
    levels=[-2.0,-1.0,-0.4,-0.2,-0.1,-0.05,0.05,0.1,0.2,0.4,1.0,2.0]
    #levels=levels/10.
    lat_lims=[-60,-90]

    custom_cmap = plt.cm.get_cmap('RdBu_r',len(levels)+1)
    custom_cmap.set_under('grey')
    custom_cmap.set_over('pink')

    cube = cubeexpt-cubecntl
    # plot
    ax=plt.subplot(1,1,1,projection=ccrs.SouthPolarStereo())
    polarCentral_set_latlim(lat_lims,ax)
    axplot=iplt.pcolormesh(cube,cmap=custom_cmap,
                 norm=mp.colors.BoundaryNorm(levels, 
                                             ncolors=len(levels)+1,
                                             clip=False))
    cbar = plt.colorbar(axplot,  orientation= 'vertical')
    cbar.set_label('psu')   
    #cbar.set_ticks([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])  
    cbar.set_ticks(levels)  
    plt.title(name.get(expt) +  '-' + name.get(cntl))
    iplt.contour(lsm_cube, levels=[-2,0.5,2],
                  colors='black', linewidths=0.1)


    plt.savefig(outstart + expt + '-' +  cntl + '_'  + '_SH_salin.png')
    plt.savefig(outstart + expt + '-' + cntl + '_'  + '_SH_salin.eps')

 




def get_avg(jobid, startyear):
    """
    gets the average data fpr the field
    """  
    
    allsalinity = np.zeros((144,288))
    count=0.
    field = 'SALINITY (OCEAN)       (PSU-35)/1000'

    # get template map
    cube = iris.load_cube('/nfs/hera1/earjcti/um/xqbwc/pg/' + 
                          'xqbwco#pg' + str(3900).zfill(9) + 'c1+.nc',
                          field)

    cube = iris.util.squeeze(cube)[LEVEL,:,:]
      
    for year in range(startyear, startyear+NYEARS):
        
        files = ('/nfs/hera1/earjcti/um/' + jobid + '/pg/' + 
                 jobid + 'o#pg' + str(year).zfill(9) + 'c1+.nc')
        print(files)
        f=MFDataset(files)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        depths = f.variables['depth_1'][:]
        sal=np.squeeze(f.variables['salinity'][:])
        sallev = sal[LEVEL,:,:]
        levreq=depths[LEVEL]

        allsalinity = allsalinity + sallev
        count=count+1.

    avgsaldata = ((allsalinity / count) * 1000.) + 35.0
       
    avgsalcube=cube.copy(data=avgsaldata)

    
    # plot average
    vals = np.arange(30.0, 38.5,0.5)

    qplt.contourf(avgsalcube,levels=vals,extend='both')
    titlename = (jobid + '. Years:' + str(startyear) + '-' + str(startyear + NYEARS) +  ' depth=' + np.str(levreq))
    plt.title(titlename, fontsize=10)
    plt.gca().coastlines()
      
    plt.savefig('/nfs/hera1/earjcti/um/' + jobid +  '/avgplots/salinity/salinity_' + jobid  + '_' + np.str(levreq) + '.eps')
    plt.savefig('/nfs/hera1/earjcti/um/' + jobid +  '/avgplots/salinity/salinity_' + jobid  + '_' + np.str(levreq) + '.png')
    
    plt.close()



    return avgsalcube,levreq
                   
    
#=================================================================
# MAIN PROGRAM STARTS HERE

name = {'xqbwc':'PI','xpsic':'PI',
        'xqbwd':'Late Pliocene', 'xqbwe':'Early Pliocene - 400ppmv',
        'xqbwr':'Late Pliocene (min_LSM)', 'xqbwg':'Early Pliocene',
        'xpsig':'Early Pliocene',
        'xpsie':'Early Pliocene 400ppmv'}



LINUX_WIN='l'
NYEARS = 100
SEASON = 'ann'

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPT = 'xqbwr'  # xpsic PI,  xpsij-lp490  xpsik - lp560
EXPT_STARTYEAR = 3900
#EXPT = 'Eoi400_ARC4_2450-2499'

# data from good experiment
CNTL = 'xqbwd'  # xpsic pi, xpsid lp400
CNTL_STARTYEAR = 3900

LEVEL=0

expt_cube, levreq = get_avg(EXPT,EXPT_STARTYEAR)

print('got expt_cube')
cntl_cube, levreq = get_avg(CNTL,CNTL_STARTYEAR)
print('got cntl cube')
diff_cube = expt_cube - cntl_cube

print(np.shape(diff_cube))

# get means
diff_cube.coord('longitude').guess_bounds()
diff_cube.coord('latitude').guess_bounds()
grid_areas = iris.analysis.cartography.area_weights(diff_cube)
meandiff = diff_cube.collapsed(['longitude','latitude'],
                               iris.analysis.MEAN, weights=grid_areas)


# get LSM
lsm_cube, lsm_plio_cube, lsm_both_cube = get_lsm(EXPT,CNTL)
# set land to -100
#MP_cube.data[0,0,:,:] = np.where(lsm_plio_cube.data == 1.0, -100.,
#                                     MP_cube.data[0,0,:,:])

# plot SH salinity cubes and anomaly
plot_SH_salin(expt_cube,lsm_cube,EXPT,EXPT_STARTYEAR,EXPT_STARTYEAR+NYEARS)
plot_SH_salin(cntl_cube,lsm_cube,CNTL,CNTL_STARTYEAR,CNTL_STARTYEAR+NYEARS)
plot_SH_salin_anom(expt_cube,cntl_cube,lsm_both_cube, EXPT,CNTL,
                   EXPT_STARTYEAR,EXPT_STARTYEAR+NYEARS)

  
