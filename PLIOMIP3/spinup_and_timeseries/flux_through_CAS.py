#!/usr/bin/env python2.7
#NAME
#    flux_through_CAS
#PURPOSE
#
#   This program will find
#   a) the waterflux through the CAS both northwards and eastwards
#   b) the ocean heat transport through the CAS relative to a reference
#
# Julia 22/11/2016



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess
import matplotlib.ticker as mticker

exptname='xpsig'
startyear=12
endyear=2999
filestart = '/uolstore/Research/a/hera1/earjcti/um/' + exptname 

# we will sum v * area from lonmin to lonmax
# and u* area from latmin to matmax.
# latmin,latmax,lonmin,lonmax are all given on the v-grid
lonmin=276.0
lonmax=283.125
vacross=13.75

latmax=13.75
latmin=11.0
uacross=283.125



#paper looks like 9N-12N inclusive
#79W-84W inclusive 281E-276E


def where_is_cas():
    """
    this program will compare the LSM with and without the CAS open
    """

    CAS_cube = iris.load_cube('/home/earjcti/um/xqbwe/pg/xqbweo#pg000003970c1+.nc','TOTAL OCEAN V-VELOCITY      CM S**-1')
    noCAS_cube = iris.load_cube('/home/earjcti/um/xqbwd/pg/xqbwdo#pg000003970c1+.nc','TOTAL OCEAN V-VELOCITY      CM S**-1')

    lons = CAS_cube.coord('longitude').points - 360.
    lats = CAS_cube.coord('latitude').points
    
    CAS_mask_cube = CAS_cube.copy(data=CAS_cube.data.mask)
    noCAS_mask_cube = noCAS_cube.copy(data=noCAS_cube.data.mask)
    print(np.shape(noCAS_cube.data.mask))
    print(np.shape(CAS_cube.data.mask))
    diff_mask_cube = noCAS_cube.copy(data=noCAS_cube.data.mask & (~CAS_cube.data.mask))  # this is the way to difference two boolean arrays

    fig, [ax0,ax1] = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    plt.sca(ax0)
    qplt.pcolormesh(CAS_mask_cube[0,0,:,:])
    ax0.set_extent([-100, -60, 0, 25], crs=ccrs.PlateCarree())
    gl = ax0.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(lons)
    gl.ylocator = mticker.FixedLocator(lats)
    gl.bottom_labels = True
    gl.right_labels = False
    ax0.xaxis.set_major_formatter(LongitudeFormatter(zero_direction_label=True))

    # draw line where we are averaging
    
    print(lons)
   
    ax0.plot(
    [-84.375, -76.875],   # Longitudes
    [13.75, 13.75],      # Latitudes (constant)
    color='red',   # Line color
    linewidth=2,   # Line width
    transform=ccrs.PlateCarree() ) # Specify the coordinate system

    ax0.plot(
    [-76.875, -76.875],   # Longitudes
    [11.25, 13.75],      # Latitudes (constant)
    color='red',   # Line color
    linewidth=2,   # Line width
    transform=ccrs.PlateCarree() ) # Specify the coordinate system

    
    #ax0.xaxis.set_visible(True)
   


    plt.sca(ax1)
    qplt.pcolormesh(diff_mask_cube[0,0,:,:])
    ax1.set_extent([-90, -68, 5, 20], crs=ccrs.PlateCarree())
    gl = ax1.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlocator = mticker.FixedLocator(lons)
    gl.ylocator = mticker.FixedLocator(lats)
    gl.bottom_labels=True
    ax1.xaxis.set_visible(True)

    # notes the CAS opennig runs across 13.75N from 

    plt.show()
    print(CAS_cube.data.mask)
    sys.exit(0)

    

def get_mean_levels(cube, weights, startlev, endlev):
    """
    get volume weighted mean between startlev and endlev
    """
    cubelevs = cube[:,startlev:endlev+1,:,:]
    weightlev = weights[:,startlev:endlev+1,: :]

    meancube = cubelevs.collapsed(['longitude','latitude','depth_1'],
                              iris.analysis.MEAN,
                              weights=weightlev)
    return meancube


def get_avg(year,weights):
    """
    gets the northwards and eastwards 
    """
 
    yearuse = str(year).zfill(9)
    filename=(filestart+'/pg/'+exptname+'o#pg'+ yearuse + 'c1+.nc')
    print(filename)
    cubes = iris.load(filename)
    Ucube = cubes.extract('TOTAL OCEAN U-VELOCITY      CM S**-1')[0]
    Vcube = cubes.extract('TOTAL OCEAN V-VELOCITY      CM S**-1')[0]
    # convert to m/s
    Vcube.data = Vcube.data / 100.0
    Ucube.data = Ucube.data / 100.0
    # this is in degC so we don't need a reference temperature
    Tcube_origgrid = iris.load_cube(filename,
                    'POTENTIAL TEMPERATURE (OCEAN)  DEG.C')
    Tcube = Tcube_origgrid.regrid(Vcube,iris.analysis.Nearest())
   


    # find northward flux from V cube across latitude=vacross and
    # between longitudes=lonmin and lonmax
    depths=Vcube.coord('depth_1').points
    lats=Vcube.coord('latitude').points
    lons=Vcube.coord('longitude').points
    # calculate the depthlayers from the midpoints
    depdiff = np.zeros(20)
    for k in range(1,19):
        depdiff[k]=(depths[k+1] - depths[k-1]) / 2.0
    depdiff[0]=depths[1]-depths[0]
    depdiff[19]=depths[19]-depths[18]
    coslat = np.cos(vacross * 2.0 * np.pi / 360.) # cos latitude of vacross
    onedeg=111100.   # 111.3 km
    
    northflux=np.zeros(20) # this is V x area
    northheatflux = np.zeros(20) 
    
    for j, lat in enumerate(lats):
        if lat == vacross:
            for i,lon in enumerate(lons):
                if lonmin < lon < lonmax:
                    #print('foundlat',j,lat,i,lon)
                    for k,dep in enumerate(depths):
                        if not Vcube.data.mask[0,k,j,i]:
                            # add on the northward flux for all depths
                            flux_add = (Vcube.data[0,k,j,i]
                                        * depdiff[k]
                                        * coslat * onedeg * (lons[1]-lons[0]))
                            northflux[k] = northflux[k] + (flux_add / 1.0E6)
                            northheatflux[k] = (northheatflux[k] +
                                      (flux_add * Tcube.data[0,k,j,i]))
                            #if k==0:
                            #    print(Vcube.data[0,k,j,i],lon,dep,flux_add / 1.0E6,northflux[k])
                            
    eastflux=np.zeros(20) # this is U x area
    eastheatflux = np.zeros(20)
    for i, lon in enumerate(lons):
        if lon == uacross:
            for j,lat in enumerate(lats):
                if latmin < lat < latmax:
                    for k,dep in enumerate(depths):
                        if not Ucube.data.mask[0,k,j,i]:
                            # add on the eastward flux for all depths
                            flux_add = (Ucube.data[0,k,j,i]
                                        * depdiff[k]
                                        * onedeg * (lats[1]-lats[0]))
                            eastflux[k] = eastflux[k] + (flux_add / 1.0E6)
                            eastheatflux[k] = (eastheatflux[k] +
                                      (flux_add * Tcube.data[0,k,j,i]))



    #print('total northward flux',np.sum(northflux),' levels ',northflux)
    #print('total eastwardflux',np.sum(eastflux),' levels ',eastflux )
    
    heatflux_multiplier = 1025 * 4000 / 1.0E15  # density * cp / petawatt conv
    tot_north_heat_flux = np.sum(northheatflux)  * heatflux_multiplier
    tot_east_heat_flux = np.sum(eastheatflux)  * heatflux_multiplier
    #print('total northward heat flux',tot_north_heat_flux)
    #print('total eastward heat flux',tot_east_heat_flux)
    #sys.exit(0)

    return (np.sum(northflux),np.sum(eastflux),tot_north_heat_flux,tot_east_heat_flux)


#######################################################
def gettransport(HadCM3,exptname,startyear,endyear):
    """
    gets the northwards and eastwards transport across the CAS between
    years startyear and endyear
    """

    # arrays for storing mean temperatures
    northflux_arr = np.zeros(endyear-startyear+1)
    eastflux_arr = np.zeros(endyear - startyear+1)
    northheat_arr = np.zeros(endyear-startyear+1)
    eastheat_arr = np.zeros(endyear - startyear+1)
    year_arr = np.zeros(endyear-startyear+1)

    # obtain means for each year and store in the arrays
    weights = [0.0]
    for year in range(startyear,endyear+1):
        (northflux,eastflux,northheat,eastheat) = get_avg(year,weights)

        northflux_arr[year-startyear]  = northflux
        eastflux_arr[year-startyear]  = eastflux
        northheat_arr[year-startyear]  = northheat
        eastheat_arr[year-startyear]  = eastheat
        year_arr[year-startyear] = year

    # plot and save
    #plt.plot(year_arr,northflux_arr + eastflux_arr,label='total')
    #plt.plot(year_arr,northflux_arr,label='north')
    #plt.legend()
    #plt.show()

    # save to file
    fileout = filestart + '/spinup/CASflow_' + exptname + '_' + str(startyear) + '_' + str(endyear) + '.tex'
    
    with open(fileout, 'w') as f:
        f.write("year,northflux,eastflux, northheatflux,eastheatflux\n")
        for x, y, z,a,b in zip(year_arr, northflux_arr, eastflux_arr,
                               northheat_arr,eastheat_arr):
            f.write(f"{x},{y},{z},{a},{b}\n")
    f.close()

    


################################
# main program


#where_is_cas()



# annual mean
figureno=0


HadCM3='y'
plt.figure(figureno)
gettransport(HadCM3,exptname,startyear,endyear)
figureno=figureno+1






sys.exit(0)

####

