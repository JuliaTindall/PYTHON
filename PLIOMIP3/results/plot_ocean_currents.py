"""
#NAME
#    plot_ocean_currents.py
#PURPOSE 
#
#  This program will plot the average ocean currents over the top n metres of the
#  ocean
"""

# Import necessary libraries

import os
import numpy as np
import math
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from matplotlib.colors import LogNorm,SymLogNorm,TwoSlopeNorm

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.quickplot as qplt
import iris.plot as iplt
from netCDF4 import Dataset, MFDataset
import sys
#from netCDF4 import Dataset, MFDataset

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



def get_avg(jobid, startyear,field):
    """
    gets the average data fpr the field
    """  

    filestart ='/home/earjcti/um/' + jobid + '/pg/' + jobid + 'o#pg'

    allcubes = CubeList([])
    
    for year in range(startyear, startyear+NYEARS):
        
        fieldcube = iris.load_cube(filestart + str(year).zfill(9) + 'c1+.nc',
                                   field)
        allcubes.append(fieldcube)
        

    iris.util.equalise_attributes(allcubes)
    onecube = allcubes.concatenate_cube()
    
    avgcube=onecube.collapsed(['t'],iris.analysis.MEAN)


    # average over top n metresm

    depth_coord = avgcube.coord('depth_1')
    nearest_depth = depth_coord.points[np.argmin(np.abs(depth_coord.points - DEPREQ))]

    top_constraint = iris.Constraint(depth_1=lambda d: d <= nearest_depth)
    shallow_cube = avgcube.extract(top_constraint)
    ndepths = len(shallow_cube.coord('depth_1').points)
   
    topdepth = np.zeros(ndepths)
    botdepth = np.zeros(ndepths)
    topdepth[0]=0.0
    for i in range(1,ndepths):
        topdepth[i]=((depth_coord.points[i] -depth_coord.points[i-1])
                     + topdepth[i-1])
        botdepth[i-1]=((depth_coord.points[i] -depth_coord.points[i-1])
                     + topdepth[i-1])
    botdepth[ndepths-1]=((depth_coord.points[ndepths-1]
                          -depth_coord.points[ndepths-2])
                         + topdepth[ndepths-1])
    
    layer_thickness = botdepth-topdepth

   
    weights = iris.util.broadcast_to_shape(layer_thickness, shallow_cube.shape, (0,))

    depthavg_cube = shallow_cube.collapsed('depth_1', iris.analysis.MEAN, weights=weights)

    
    return depthavg_cube


def plot_ocean_currents(u_cube,v_cube,jobid):


    lats = u_cube.coord('latitude').points
    lons = u_cube.coord('longitude').points

    # Extract data arrays (e.g., surface layer or collapsed depth)
    u = u_cube.data / 100.
    v = v_cube.data / 100.

    # i am going to convert all vectors to length 1 so
    #  we can see the direction of the currents only
    u_unit = np.zeros_like(u)
    v_unit = np.zeros_like(v)
    mag = np.hypot(u, v)
    mask = mag > 0
    u_unit[mask] = u[mask] / (mag[mask] * 5.)
    v_unit[mask] = v[mask] / (mag[mask] * 5.)
   
    # Compute speed for color scale
    speed = np.sqrt(u**2 + v**2)

    # Create meshgrid for plotting
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Plot quiver
    plt.figure(figsize=(10, 6))
    #Q = plt.quiver(lon_grid, lat_grid, u, v, speed, cmap='viridis', scale=5)
    #contour = plt.contourf(lons, lats, speed, cmap='viridis')


    levels = np.logspace(-3, 0, num=10)  # 10 divisions between 10^-3 and 10^0
    contour = plt.contourf(lons, lats, speed, cmap='viridis',levels=levels,
                       norm=LogNorm(vmin=speed.min()+1e-3, vmax=speed.max()))

    step=8
    plt.quiver(lon_grid[::step, ::step], lat_grid[::step, ::step],
           u_unit[::step, ::step], v_unit[::step, ::step],
           color='black', scale=10)


    mag = np.hypot(u, v)
    mask = mag > 0.1  # or your anomaly threshold
    ii, jj = np.where(mask)
    plt.quiver(lon_grid[ii, jj], lat_grid[ii, jj],
           u_unit[ii, jj], v_unit[ii, jj],
           color='black', scale=10)

    
    cbar=plt.colorbar(contour, label='Speed (m/s)',ticks=levels)
    cbar.ax.yaxis.set_major_formatter(mticker.LogFormatter(base=10, labelOnlyBase=False))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Ocean Currents: ' + jobid)
        
    plt.show()

def basin_currents(u_cube,v_cube,jobid,basin):
    """
    plots currents only for a basin so we can see them better
    """
    
    lats = u_cube.coord('latitude').points
    lons = u_cube.coord('longitude').points

    # Extract data arrays (e.g., surface layer or collapsed depth)
    u = u_cube.data / 100.
    v = v_cube.data / 100.

    
    if basin == 'Atlantic':
        lons = np.where(lons > 180, lons - 360, lons)
        xmin=-100
        xmax=20
        ymin=-60
        ymax=60


    # set up basin

    proj_map = ccrs.PlateCarree()   # map projection
    proj_data = ccrs.PlateCarree()  # your data is in lon/lat degrees

    fig = plt.figure(figsize=(11, 7))
    ax = plt.axes(projection=proj_map)

    ax.set_extent([xmin, xmax, ymin, ymax], crs=proj_data)

    # --- Coastlines and land ---
    ax.coastlines(resolution='110m', linewidth=1.0)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=0)
    ax.add_feature(cfeature.OCEAN, facecolor='white', zorder=0)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, edgecolor='gray')

    # --- Gridlines ---
    gl = ax.gridlines(crs=proj_data, draw_labels=True,
                  linewidth=0.6, color='gray', alpha=0.5, linestyle='--')
    gl.right_labels = False
    gl.top_labels = False


  
    # i am going to convert all vectors to length 1 so
    #  we can see the direction of the currents only
    u_unit = np.zeros_like(u)
    v_unit = np.zeros_like(v)
    mag = np.hypot(u, v)
    mask = mag > 0
    u_unit[mask] = u[mask] / (mag[mask] * 5.)
    v_unit[mask] = v[mask] / (mag[mask] * 5.)
   
    # Compute speed for color scale
    speed = np.sqrt(u**2 + v**2)

    # Create meshgrid for plotting
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Plot quiver
    #Q = plt.quiver(lon_grid, lat_grid, u, v, speed, cmap='viridis', scale=5)
    #contour = plt.contourf(lons, lats, speed, cmap='viridis')


    levels = np.logspace(-3, 0, num=10)  # 10 divisions between 10^-3 and 10^0
    contour = ax.pcolormesh(lon_grid, lat_grid,
                            speed, cmap='viridis',
                       norm=LogNorm(vmin=speed.min()+1e-3, vmax=speed.max()),
                                    transform=proj_data,zorder=1)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=3)

    step=2
    ax.quiver(lon_grid[::step, ::step], lat_grid[::step, ::step],
           u_unit[::step, ::step], v_unit[::step, ::step],
           color='black', scale=10)


    #mag = np.hypot(u, v)
    #mask = mag > 0.1  # or your anomaly threshold
    #ii, jj = np.where(mask)
    #ax.quiver(lon_grid[ii, jj], lat_grid[ii, jj],
    #       u_unit[ii, jj], v_unit[ii, jj],
    #       color='black', scale=10)

    
    cbar=plt.colorbar(contour, label='Speed (m/s)',ticks=levels)
    cbar.ax.yaxis.set_major_formatter(mticker.LogFormatter(base=10, labelOnlyBase=False))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Ocean Currents: ' + jobid)
        
    plt.show()


def plot_ocean_currents_anom(u_cube,u_cntl,v_cube,v_cntl):
    """
    plots the anomaly between u and v
    """

    lats = u_cube.coord('latitude').points
    lons = u_cube.coord('longitude').points

    # Extract data arrays (e.g., surface layer or collapsed depth)
    u = u_cube.data / 100.
    v = v_cube.data / 100.
    speed_expt = np.sqrt(u**2 + v**2)

    uctl = u_cntl.data / 100.
    vctl = v_cntl.data / 100.
    speed_cntl = np.sqrt(uctl**2 + v**2)

    speed_anom = speed_expt - speed_cntl

    # we are going to plot quivers (from the experiment not the anomaly)
    # only on points where the magnitude of the anomaly is > 0.01
    
    # threshold for plotting arrows
    thresh = 0.01
    mask = np.abs(speed_anom) > thresh   # True where we want arrows

    # Normalize u, v to unit vectors (direction only)
    mag = np.hypot(u, v)
    u_unit = np.zeros_like(u)
    v_unit = np.zeros_like(v)
    valid = mag > 0
    u_unit[valid] = u[valid] / mag[valid]
    v_unit[valid] = v[valid] / mag[valid]

    # Choose a constant arrow length in data units (tune as needed)
    dx = np.median(np.diff(lons, axis=1)) if lons.ndim == 2 else np.median(np.diff(lons))
    dy = np.median(np.diff(lats, axis=0)) if lats.ndim == 2 else np.median(np.diff(lats))
    L = 0.8 * min(dx, dy)

    u_plot = u_unit * L * 8
    v_plot = v_unit * L * 8

    # Create meshgrid for plotting
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    
    # Mask everything except where 'mask' is True
    u_plot_m = np.ma.masked_where(~mask, u_plot)
    v_plot_m = np.ma.masked_where(~mask, v_plot)
    lons_m    = np.ma.masked_where(~mask, lon_grid)
    lats_m    = np.ma.masked_where(~mask, lat_grid)


   
    
    
    # do a log anomaly plot and add quivers
    plt.figure(figsize=(10, 6))
    vmin, vmax = -0.1, 0.1
    #norm = SymLogNorm(linthresh=1e-2, vmin=vmin, vmax=vmax)
    #levels = np.concatenate([
    #        -np.logspace(-3, np.log10(vmax), 10),
    #        np.logspace(-3, np.log10(vmax), 10)
    #    ])
    #levels=np.sort(levels)

    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0,vmax=vmax)
    levels = np.linspace(vmin,vmax,21)
 
    contour = plt.contourf(
        lons, lats, speed_anom,
        cmap='RdBu_r', 
        norm=norm,
        levels=levels,extend='both')

    cbar = plt.colorbar(contour)
    cbar.set_label("Anomaly")

    step = 5
    sl = (slice(None, None, 2), slice(None, None, 10))

    
    Q = plt.quiver(
    lons_m[sl], lats_m[sl], u_plot_m[sl], v_plot_m[sl],
    angles='xy',
    scale_units='xy',
    scale=1,
    pivot='mid',
    color='k',
    width=0.002
)


    plt.show()
    sys.exit(0)

    #Q = plt.quiver(lon_grid, lat_grid, u, v, speed, cmap='viridis', scale=5)
    #contour = plt.contourf(lons, lats, speed, cmap='viridis')


    #levels = np.logspace(-3, 0, num=10)  # 10 divisions between 10^-3 and 10^0
    contour = plt.contourf(lons, lats, speed_anom, cmap='RdBu_r')

    #step=8
    #plt.quiver(lon_grid[::step, ::step], lat_grid[::step, ::step],
    #       u_unit[::step, ::step], v_unit[::step, ::step],
    #      color='black', scale=10)

    cbar=plt.colorbar(contour, label='Speed (m/s)')
    #cbar.ax.yaxis.set_major_formatter(mticker.LogFormatter(base=10, labelOnlyBase=False))
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title('Ocean Currents: ' + EXPT + '-' + CNTL)
        
    plt.show()

#=================================================================
# MAIN PROGRAM STARTS HERE



LINUX_WIN='l'
NYEARS = 30
SEASON = 'ann'
DEPREQ=200

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPT = 'xqbwg'  # xpsic PI,  xpsij-lp490  xpsik - lp560
CNTL = 'xqbwc'
EXPT_STARTYEAR = 3970
#EXPT = 'Eoi400_ARC4_2450-2499'


v_cube = get_avg(EXPT,EXPT_STARTYEAR,'TOTAL OCEAN V-VELOCITY      CM S**-1')
u_cube = get_avg(EXPT,EXPT_STARTYEAR,'TOTAL OCEAN U-VELOCITY      CM S**-1')

#plot_ocean_currents(u_cube,v_cube,EXPT)
basin_currents(u_cube,v_cube,EXPT,'Atlantic')
sys.exit(0)
v_cntl = get_avg(CNTL,EXPT_STARTYEAR,'TOTAL OCEAN V-VELOCITY      CM S**-1')
u_cntl = get_avg(CNTL,EXPT_STARTYEAR,'TOTAL OCEAN U-VELOCITY      CM S**-1')


plot_ocean_currents_anom(u_cube,u_cntl,v_cube,v_cntl)



print(u_cube)
sys.exit(0)

print('got cntl cube')
 
for EXPT in EXPTS:
    expt_cube = get_avg(EXPT,EXPT_STARTYEAR)
    print('got expt_cube')
    diff_cube = expt_cube - cntl_cube

    # get means
    diff_cube.coord('longitude').guess_bounds()
    diff_cube.coord('latitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(diff_cube)
    meandiff = diff_cube.collapsed(['longitude','latitude'],
                                   iris.analysis.MEAN, weights=grid_areas)
    diffchar = str(np.around(meandiff.data,2))
   
    print('about to plot')
    vals = np.arange(-4.0, 4.5,0.5)
    #vals = np.arange(-8.0, 9.0, 1.0)

    if CNTL == 'xpsid' or CNTL == 'xpsig' or CNTL == 'xpsie':
        vals = np.arange(-4.0, 4.5,0.5)
#    if CNTL == 'xpsic' or CNTL == 'xqbwc' or CNTL == 'xqbwi':
#        vals = np.arange(-8.0, 9.0, 1.0)
    if EXPT == 'xpsin' or EXPT=='xpsio' or EXPT == 'xpsip' or EXPT == 'xpsiq':
        vals = np.arange(-6.0, 7.0, 1.0)
    qplt.contourf(diff_cube,levels=vals,extend='both',cmap='RdBu_r')
    titlename = EXPT + '-' +  CNTL + '. Years:' + str(EXPT_STARTYEAR) + '-' + str(EXPT_STARTYEAR + NYEARS) + '.Meandiff =' +  diffchar
    plt.title(titlename, fontsize=10)
    plt.gca().coastlines()
      
    print('about to write to file')
    plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja_' + EXPT + '-' + CNTL + '_' + field + '.eps')
    plt.savefig('/nfs/hera1/earjcti/um/' + EXPT +  '/avgplots/jja_' + EXPT + '-' + CNTL + '_' + field + '.png')
    plt.close()
