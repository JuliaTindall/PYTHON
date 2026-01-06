#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created October 2024 by Julia

This program will setup up the SST and sea ice for HadGEM2
"""

import numpy as np
import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib
import matplotlib.pyplot as plt
from iris.cube import CubeList
#import matplotlib.ticker as mticker
#import netCDF4
#from datetime import date
#import cmocean to use topo colormap
import sys

def plot_cube_diff(cube1, cube2, title1, title2, titleanom):
    """
    plots two cube1 and the anomaly between them over the globe
    """

    plt.figure(figsize=[12.0, 8.0])
    diff_cube = cube1.copy(cube2.data - cube1.data)
    print('number of non zero points',np.sum(np.where((diff_cube.data) != 0.0, 1.0, 0.0)))
    print('number of positive points',np.sum(np.where((diff_cube.data) > 0.5, 1.0, 0.0)))
    print('number of negative points',np.sum(np.where((diff_cube.data) < -0.5, 1.0, 0.0)))
    vals=np.arange(-1.1,1.3,0.2)  
    plt.subplot(2,2,1)
    qplt.contourf(cube1,cmap='bwr',levels=vals)
    plt.gca().coastlines()
    plt.title(title1)
    plt.subplot(2,2,2)
    qplt.contourf(cube2,cmap='bwr',levels=vals)
    plt.gca().coastlines()
    plt.title(title2)
    plt.subplot(2,2,3)
    qplt.contourf(diff_cube,cmap='bwr',levels=vals)
    plt.gca().coastlines()
    plt.title(titleanom)
    plt.show()
    sys.exit(0)


def plot_cmesh(cube1, cube2, cube3, title1, title2, title3):
    """
    plots two cube1 and the anomaly between them over the globe
    """
    boundaries=np.arange(-1.1,1.3,0.2)  
  
    plt.figure(figsize=[12.0, 8.0])
    vals=np.arange(-1.1,1.3,0.2)  
  
    cmap_bwr = plt.cm.get_cmap('bwr',len(boundaries))
  

    ax=plt.subplot(2,2,1,projection=ccrs.PlateCarree())
    cs=iplt.pcolormesh(cube1, cmap=cmap_bwr,
                       norm = matplotlib.colors.BoundaryNorm(boundaries, 
                                              ncolors=len(boundaries)-1, 
                                              clip=False))
  #  cbar=plt.colorbar(cs,orientation='vertical', extend='min')
    cbar=plt.colorbar(cs,orientation='vertical')
    cbar.set_ticks(boundaries)
    ax.gridlines(xlocs=cube1.coord('longitude').points-358.125, ylocs=cube1.coord('latitude').points-1.25)
    ax.coastlines()
    plt.title(title1)
   
  
    ax=plt.subplot(2,2,2,projection=ccrs.PlateCarree())
    cs=iplt.pcolormesh(cube2, cmap=cmap_bwr,
                       norm = matplotlib.colors.BoundaryNorm(boundaries, 
                                              ncolors=len(boundaries)-1, 
                                              clip=False))
    cbar=plt.colorbar(cs,orientation='vertical')
    cbar.set_ticks(boundaries)
    ax.gridlines(xlocs=cube2.coord('longitude').points-0.5, ylocs=cube2.coord('latitude').points-0.5)
    ax.coastlines()
    plt.title(title2)
   
    ax=plt.subplot(2,2,3,projection=ccrs.PlateCarree())
    cs=iplt.pcolormesh(cube3, cmap=cmap_bwr,
                       norm = matplotlib.colors.BoundaryNorm(boundaries, 
                                              ncolors=len(boundaries)-1, 
                                              clip=False))
  #  cbar=plt.colorbar(cs,orientation='vertical', extend='min')
    cbar=plt.colorbar(cs,orientation='vertical')
    cbar.set_ticks(boundaries)
    ax.gridlines(xlocs=cube3.coord('longitude').points-358.125, ylocs=cube3.coord('latitude').points-1.25)
    ax.coastlines()
    plt.title(title1)
 
   
    plt.show()


def make_lsmfrac_LP_1():
    """
    this program will 
    1. read in the preindustrial ancil file from HadGEM2
        (/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/mask_frac_l32161_n96_ref.nc)
    2. read in the preindustrial and Late Pliocene mask files
        (/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/
        Modern_std_exp_data/Modern_std_LSM_v1.0.nc
        LP_exp_data/LP_LSM_v1.0.nc
    3.  Put both the pliomip3 boundary conditions on the HadGEM grid and
        find out which gridboxes have changed (Use linear interpolation)
    4.  For the gridboxes that have changed, change them in the mask_frac file
        cube. 
    5.  Partially correct.  Over Hudson Bay, Australia, CAA 
    6.  write out to a file.    
    """

    # 1. 2. 3 read in files and regrid as necessary  
    # also read in the preindustrial ancil files because we want to see
    # what has changed between the pliocene and the preindustrial

    # PI HadGEM2 put on a 0-360 deg grid
    filename = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/mask_frac_l32161_n96_ref.nc'
    PI_cube = iris.load_cube(filename)
    PI_cube.coord('longitude').guess_bounds()
    PI_cube.coord('latitude').guess_bounds()
    
  
    
    # PLIOMIP3 boundary conditions for mask and regrid
    startPL2 = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'
    mod = 'Modern_std/Modern_std/Modern_std_'
    P2 = 'Plio_enh/Plio_enh/Plio_enh_'
    mask_ungridded_pi_cube = iris.load_cube(startPL2 + mod + 'LSM_v1.0.nc')
    mask_ungridded_p2_cube = iris.load_cube(startPL2 + P2 + 'LSM_v1.0.nc')
    mask_ungridded_pi_cube.coord('longitude').circular = True
    mask_ungridded_p2_cube.coord('longitude').circular = True
  
    mask_ungridded_pi_cube.coord('longitude').guess_bounds()
    mask_ungridded_pi_cube.coord('latitude').guess_bounds()
    mask_ungridded_p2_cube.coord('longitude').guess_bounds()
    mask_ungridded_p2_cube.coord('latitude').guess_bounds()

    mask_pliomip2_pi_cube = mask_ungridded_pi_cube.regrid(PI_cube,iris.analysis.AreaWeighted(mdtol=0.5))
    mask_pliomip2_p2_cube = mask_ungridded_p2_cube.regrid(PI_cube,iris.analysis.AreaWeighted(mdtol=0.5))
    masklin_pliomip2_pi_cube = mask_ungridded_pi_cube.regrid(PI_cube,iris.analysis.Linear())
    masklin_pliomip2_p2_cube = mask_ungridded_p2_cube.regrid(PI_cube,iris.analysis.Linear())
    masklin_pliomip2_p2_cube.long_name = 'P2 regridded linear interpolation'
    mask_diff_p2_and_hg2 = PI_cube - mask_pliomip2_pi_cube
    outfile = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/diffmaskfrac_pi_lsm.nc'
   
    # which gridboxes changed
    diff_cube = mask_pliomip2_p2_cube - mask_pliomip2_pi_cube
    diff_cube.long_name = 'gridded_p2-pi from PlioMIP'
    ungridded_diff_cube = mask_ungridded_p2_cube - mask_ungridded_pi_cube
    ungridded_diff_cube.long_name = 'ungridded p2-pi from PlioMIP'
    iris.save([mask_ungridded_pi_cube, mask_ungridded_p2_cube,
               ungridded_diff_cube,
                mask_pliomip2_pi_cube,mask_pliomip2_p2_cube, 
                PI_cube,mask_diff_p2_and_hg2,diff_cube,
               masklin_pliomip2_p2_cube],
              outfile, netcdf_format='NETCDF3_CLASSIC')
  

    #4. add difference onto the PI  
    diff_data = diff_cube.data
    hg_data = np.copy(PI_cube.data + diff_data)

    # corrections

    # fractions > 1.0 or < 0.0fractions
    hg_data = np.where(hg_data <  0.0, 0.0, hg_data)
    hg_data = np.where(hg_data > 1.0, 1.0, hg_data)

    lats = PI_cube.coord('latitude').points
    lons = PI_cube.coord('longitude').points
   
    cube_int1_boundscorr = PI_cube.copy(data=np.copy(hg_data))
    outfile = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/P4_maskfrac_unfixed.nc'
    iris.save(cube_int1_boundscorr,outfile)


    linear_pliomip2_data = masklin_pliomip2_p2_cube.data
    # manual corrections
  

  # correct canadian arctic archipelago and other values in hudson bay
    for j,lat in enumerate(lats):    
   
        if lat>75:   
            for i,lon in enumerate(lons):
                if 260<lon<300:
                    hg_data[0,0,j,i] = mask_pliomip2_p2_cube.data[j,i]
        if 72 <lat < 77:
            for i,lon in enumerate(lons):
                if 268<lon<273:        
                    if hg_data[0,0,j,i] < 1.0: # if frac lsm or zero
                        hg_data[0,0,j,i]=1.0
   
    # if gridbox is more than 50% land and is not next to an ocean point
    # set to land but try to move to nearby ocean
    for j, lat in enumerate(lats):
        jmin1=j-1
        if jmin1 < 0: jmin1=jmin1+145
        jpl1=j+1
        if jpl1 > 144: jmin1=jmin1-145

        for i,lon in enumerate(lons):
            if 0.5 < hg_data[0,0,j,i]:
                imin1=i-1
                if imin1 < 0: imin1=imin1+192
                ipl1=i+1
                if ipl1 > 191: ipl1=ipl1-192
                if (hg_data[0,0,jmin1,i] > 0.0 and
                    hg_data[0,0,jpl1,i] > 0.0 and
                    hg_data[0,0,j,imin1] > 0.0 and
                    hg_data[0,0,j,ipl1] > 0.0):
                    to_remove = 1.0 - hg_data[0,0,j,i]
                    hg_data[0,0,j,i]=1.0

                    
                    # try and move to nearby coast if possible
                    # 1. find out which coast
                    coast_points=np.zeros(4)  # north, south, east, west
                    coast_ind=np.zeros(4,dtype=bool) 
                    for dist in range(6,1,-1):
                        if hg_data[0,0,j+dist,i] == 0.0:
                            coast_points[0] = coast_points[0] + 1
                            coast_ind[0]=True
                        if hg_data[0,0,j-dist,i] == 0.0:
                            coast_points[1] = coast_points[1] + 1
                            coast_ind[1]=True
                        if hg_data[0,0,j,(i-dist)%192] == 0.0:
                            coast_points[2] = coast_points[2] + 1
                            coast_ind[2]=True
                        if hg_data[0,0,j,(i+dist)%192] == 0.0:
                            coast_points[3] = coast_points[3] + 1
                            coast_ind[3]=True
                    print(coast_points)
                    print(coast_ind)
                    # if multiple coasts have been found drop down to 1
                    if np.sum(coast_ind) > 1.0:
                        ind=np.argmax(coast_points)
                        coast_ind[:]=False
                        coast_ind[ind]=True
                        print('multiple coasts found to use',coast_ind)
                      
                    # if north coast nearest move the stuff northwards
                    index=1
                    if coast_ind[0] == True:
                        while to_remove > 0 and index < 7:
                            if hg_data[0,0,j+index,i] < 1.0:
                                can_decrease = hg_data[0,0,j+index,i]  
                                hg_data[0,0,j+index,i] = (
                                    hg_data[0,0,j+index,i] -
                                         min(to_remove,can_decrease))
                                to_remove = (to_remove -
                                             min(to_remove, can_decrease))

                            index=index+1

                    # if south coast nearest move the stuff southwards
                    index=1
                    if coast_ind[1] == True:
                        while to_remove > 0 and index < 7:
                            if hg_data[0,0,j-index,i] < 1.0:
                                can_decrease = hg_data[0,0,j-index,i]  
                                hg_data[0,0,j-index,i] = (
                                    hg_data[0,0,j-index,i] -
                                         min(to_remove,can_decrease))
                                to_remove = (to_remove -
                                             min(to_remove, can_decrease))
                            index=index+1

                    # if west coast nearest move the stuff west
                    index=1
                    if coast_ind[2] == True:
                        while to_remove > 0 and index < 7:
                            if hg_data[0,0,j,i-index] < 1.0:
                                can_decrease = hg_data[0,0,j,i-index]  
                                hg_data[0,0,j,i-index] = (
                                    hg_data[0,0,j,i-index] -
                                         min(to_remove,can_decrease))
                                to_remove = (to_remove -
                                             min(to_remove, can_decrease))

                            index=index+1

                    # if east coast nearest move the stuff east
                    index=1
                    if coast_ind[3] == True:
                        while to_remove > 0 and index <7:
                            ixuse=(i+index)%192
                            if hg_data[0,0,j,ixuse] < 1.0:
                                can_decrease = hg_data[0,0,j,ixuse]  
                                hg_data[0,0,j,ixuse] = (
                                    hg_data[0,0,j,ixuse] -
                                         min(to_remove,can_decrease))
                                to_remove = (to_remove -
                                             min(to_remove, can_decrease))

                            index=index+1
                            
                            
                    
    # don't allow isolated fractional gridboxes in high northern latitudes
    for j,lat in enumerate(lats):    
        # set points to land if they are fractional and mid continent
        if 40<lat<85:   
            for i,lon in enumerate(lons):
              if 5 <lon < 350:
                if hg_data[0,0,j,i] < 1.0: # if frac lsm
             
                    # if land below and above set to land - move the land
                    # fraction if possible
                    if hg_data[0,0,j+1,i] == 1.0 and hg_data[0,0,j-1,i] == 1.0:
                        to_remove = 1.0- hg_data[0,0,j,i]
                        hg_data[0,0,j,i]=1.0
                        if hg_data[0,0,j,i-1] < 1.0:
                            hg_data[0,0,j,i-1] = max(0.0, 
                                          hg_data[0,0,j,i-1] - to_remove)
                        elif hg_data[0,0,j,i+1] < 1.0:
                            hg_data[0,0,j,i+1] = max(0.0, 
                                          hg_data[0,0,j,i+1] - to_remove)
                            
                    # if to right and left set to land
                    if hg_data[0,0,j,i+1] == 1.0 and hg_data[0,0,j,i-1] == 1.0:
                        to_remove = 1.0- hg_data[0,0,j,i]
                        hg_data[0,0,j,i]=1.0
                        if hg_data[0,0,j-1,i-1] < 1.0:
                            hg_data[0,0,j-1,i] = max(0.0, 
                                          hg_data[0,0,j-1,i] - to_remove)
                        elif hg_data[0,0,j+1,i] < 1.0:
                            hg_data[0,0,j+1,i] = max(0.0, 
                                          hg_data[0,0,j+1,i] - to_remove)


    
    # correction for bering strait not needed
    # correction near california not needed
    
    # manual correction near australia
        if -16<lat<-13:
            for i,lon in enumerate(lons):
                if 133 < lon < 138:
                    hg_data[0,0,j,i]=1.0
      
    

    P2_maskfrac_cube = PI_cube.copy(data=hg_data)
    diffcube = P2_maskfrac_cube - cube_int1_boundscorr
    print(diffcube)
    diffcube.long_name='partially corrected minus uncorrected'
   
    outfile = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/P4_maskfrac_partfixed.nc'
    iris.save([P2_maskfrac_cube,diffcube],outfile,netcdf_format='NETCDF3_CLASSIC')
    
    


   
def make_lsmfrac_LP_2():
    """
    this program will 
    1. read in the partially corrected land sea mask produced in make_lsmfrac_LP_1
    2. read in the regridded pliocene mask
    2. Correct over antarctica
    """

    lsm_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/P4_maskfrac_partfixed.nc','Land fraction in grid box')
    lsm_p2_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/diffmaskfrac_pi_lsm.nc','p4_lsm_0')
    lsm_p2linint_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/diffmaskfrac_pi_lsm.nc','P2 regridded linear interpolation')
    filename = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/mask_frac_l32161_n96_ref.nc'
    PI_cube = iris.load_cube(filename)


  
    # overwrite Antarctica with the P2 data
    lons = lsm_cube.coord('longitude').points
    lats = lsm_cube.coord('latitude').points
    lsm_cube_data = lsm_cube.data
    lsm_p2_cube_data = lsm_p2_cube.data

    for j,lat in enumerate(lats):
        if lat < -60.:
            for i,lon in enumerate(lons):
                if 150. < lon < 350.:
                    lsm_cube_data[0,0,j,i] = lsm_p2_cube_data[j,i]
    # corect some fractional gridboxes in the middle of Antarctica
        if -73<lat<-67:
            for i,lon in enumerate(lons):
                if 104<lon<113:
                    lsm_cube_data[0,0,j,i]=1.0


    # connect some isolated sea points south east asia with the ocean

    for j,lat in enumerate(lats):
        if 28.5 < lat < 35.:
            for i,lon in enumerate(lons):
                if 123. < lon < 124.:
                    lsm_cube_data[0,0,j,i] = 0.0
   
    # correct some jagged edges in east antarctica

    for j,lat in enumerate(lats):
        if -75 < lat < -66.:
            for i,lon in enumerate(lons):
                if 127. < lon < 142.0:
                    lsm_cube_data[0,0,j,i] = lsm_cube_data[0,0,j,i-1]
                if 143. < lon < 147.0:
                    if 0.2 < lsm_cube_data[0,0,j,i] <  0.9:
                        lsm_cube_data[0,0,j,i] = lsm_cube_data[0,0,j,i+1]
                if 142. < lon < 143.0:
                    if ((0.1 < lsm_cube_data[0,0,j-1,i] <  0.5) and
                        (lsm_cube_data[0,0,j,i]==1.0)):
                        lsm_cube_data[0,0,j,i] = lsm_cube_data[0,0,j-1,i]

      
    #alter a couple of gridpoints to be same as preindustrial in CAS regino.
    for j,lat in enumerate(lats):
        if 7 < lat < 16.:
            for i,lon in enumerate(lons):
                if 275. < lon < 276.0:
                    lsm_cube_data[0,0,j,i] = PI_cube.data[0,0,j,i]

    # alter red sea region to be same as in the Pliocene standard
    # but try and reduce number of fractional points
    # persian gulf closed in original pliocene so keep as it is
    ixlats = []
    for j,lat in enumerate(lats):
        if 10 < lat < 22.:
            ixlats.append(j)
            for i,lon in enumerate(lons):
                if 35. < lon < 46.0:
                    lsm_cube_data[0,0,j,i] = lsm_p2_cube_data[j,i]
    print(ixlats)
    for i,lon in enumerate(lons):
        if 36. < lon < 46.0:
          for sweep in range(0,5):  #5 sweeps  
            vals=[]
            valsix=[]
            for j in ixlats:
                if 0 < lsm_cube_data[0,0,j,i] < 1.0:
                    vals.append(lsm_cube_data[0,0,j,i])
                    valsix.append(j)
            if len(vals) > 1:  # check there is enough in the array
                maxval=max(vals)
                maxix=vals.index(max(vals))
                minval=min(vals)
                minix=vals.index(min(vals))
                maxtomove=1.0-maxval
                maxtotake=minval
                if maxtomove > maxtotake:
                    lsm_cube_data[0,0,valsix[minix],i]=0.0
                    lsm_cube_data[0,0,valsix[maxix],i]=(
                        lsm_cube_data[0,0,valsix[maxix],i] +
                        lsm_cube_data[0,0,valsix[minix],i])
                else:
                    lsm_cube_data[0,0,valsix[maxix],i]=1.0
                    lsm_cube_data[0,0,valsix[minix],i]= (
                        lsm_cube_data[0,0,valsix[minix],i] - maxtomove)
                                                              

   

    
    outfile = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/temporary_P4_maskfrac.nc'
    subsetcubes = CubeList([])
    subsetcubes.append(lsm_cube)
    for i in range(0,6):
        lonmin = (i * 60)-10.
        lonmax = ((i+1)*60) + 10.
        #print(lsm_cube)
        loncon = iris.Constraint(longitude=lambda cell: lonmin < cell < lonmax)
        loncube = lsm_cube.extract(loncon)
        #print(lonmin, lonmax)
        #print(lsm_cube.coord('longitude').points)
        #print(loncube)
        for j in range(0,3):
            latmin=((j*60) -90.) - 10.
            latmax = (((j+1) * 60.)-90) + 10.0
            latcon = iris.Constraint(latitude=lambda cell: latmin < cell < latmax)
            cube = loncube.extract(latcon)
            cube.long_name = 'region' + str(i) + '_' + str(j)
            subsetcubes.append(cube)
        
    # extract a preindustrial for comparison
    lonmin=290.
    lonmax=360.
    latmin=20.
    latmax=90
    loncon=iris.Constraint(longitude=lambda cell: lonmin < cell < lonmax)
    latcon = iris.Constraint(latitude=lambda cell: latmin < cell < latmax)
    loncube = PI_cube.extract(loncon)
    cube=loncube.extract(latcon)
    cube.long_name='preindustrial region 1_1'
    print('julia',cube)
    subsetcubes.append(cube)

    iris.save(subsetcubes,outfile,netcdf_format='NETCDF3_CLASSIC')
        
    
    outfile = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/P4_maskfrac.nc'
    iris.save(lsm_cube,outfile,netcdf_format='NETCDF3_CLASSIC')


def checks_on_mask():
    """
    check to see whether the pliocene - preindustrial anomaly for the Pliomip3
    boundary conditions looks like the pliocene-preindustrial anomaly that
    I made for HadGEM2
    """

    # read in HadGEM2 Pliocene and preindustrial masks

    HadGEM_P2_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/P4_maskfrac.nc','Land fraction in grid box')
    HadGEM_PI_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/mask_frac_l32161_n96_ref.nc')

    startPL2 = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'
    mod = 'Modern_std/Modern_std/Modern_std_'
    P2 = 'Plio_enh/Plio_enh/Plio_enh_'
    pliomip_pi_cube = iris.load_cube(startPL2 + mod + 'LSM_v1.0.nc')
    pliomip_p2_cube = iris.load_cube(startPL2 + P2 + 'LSM_v1.0.nc')
    pliomip_pi_cube.coord('longitude').circular = True
    pliomip_p2_cube.coord('longitude').circular = True

    # get anomalies on raw grid
    HadGEM2_anom = HadGEM_P2_cube - HadGEM_PI_cube
    HadGEM2_anom.long_name = 'HadGEM2 anomaly'
    pliomip_anom = pliomip_p2_cube - pliomip_pi_cube
    pliomip_anom.long_name= 'Pliomip2 anomaly (original grid)'

    # regrid the pliomip2 stuff and calculate anomalies
  
    pliomip_pi_cube.coord('longitude').guess_bounds()
    pliomip_pi_cube.coord('latitude').guess_bounds()
    pliomip_p2_cube.coord('longitude').guess_bounds()
    pliomip_p2_cube.coord('latitude').guess_bounds()
    regrid_pliomip2_pi_cube = pliomip_pi_cube.regrid(HadGEM_P2_cube,iris.analysis.AreaWeighted(mdtol=0.5))
    regrid_pliomip2_p2_cube = pliomip_p2_cube.regrid(HadGEM_P2_cube,iris.analysis.AreaWeighted(mdtol=0.5))
    anom_p2_grid = regrid_pliomip2_p2_cube - regrid_pliomip2_pi_cube
    anom_p2_grid.long_name = 'PlioMIP2 anomaly (regridded)'

    # anomaly of anomalies
    anom_anom_cube = HadGEM2_anom - anom_p2_grid
    anom_anom_cube.long_name = 'HadGEM anomalies - pliomip2 anomalies'

    
    outfile = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/checks_on_mask.nc'
    iris.save([HadGEM2_anom, pliomip_anom, anom_p2_grid,anom_anom_cube],outfile,netcdf_format='NETCDF3_CLASSIC')


    # what is average value of everything
    print('HadGEM2 changes=',HadGEM2_anom.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)
    print('PlioMIP2 changes (orig grid)=',pliomip_anom.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)
    print('PlioMIP2 changes regridded=',anom_p2_grid.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)
    print('difference=',anom_anom_cube.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)
    print(' ')
    print('PlioMIP2 PI land=',pliomip_pi_cube.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)
    print('PlioMIP2 LP land=',pliomip_p2_cube.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)
    print('HadGEM PI land=',HadGEM_PI_cube.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)
    print('HadGEM LP land=',HadGEM_P2_cube.collapsed(['longitude','latitude'],iris.analysis.MEAN).data)


    




def make_other_ancil_EP():
    """
    make all the other ancils.  (Done LSM and rivers in earlier program)
    """
    filestart = '/nfs/hera1/earjcti/ancil/P4_enh/'
    pft_file = filestart + 'P4_enh_mb_qrparm.pft.nc'
    soil_file = filestart + 'P4_enh_mb_qrparm.soil.nc'
    slt_file = filestart + 'P4_enh_mb_qrclim.slt.nc'
    smow_file = filestart + 'P4_enh_mb_qrclim.smow.nc'
    disturb_file = filestart + 'P4_enh_qrfrac.disturb.nc'
    fractype_file = filestart + 'P4_enh_mb_qrfrac.type.nc'
    orogtype_file = filestart + 'P4_enh_qrparm.orog.nc'
    filename_dict = {'pft':pft_file,   'soil':soil_file, 'slt':slt_file,
                     'smow':smow_file,  'disturb':disturb_file,
                     'frac':fractype_file, 'orog':orogtype_file}
    fileout_dict = {'pft':'EP_pft.nc',   'soil':'EP_soil.nc', 
                    'slt': 'EP_slt.nc', 'smow':'EP_smow.nc', 
                    'disturb':'EP_disturb.nc', 'frac':'EP_frac.nc', 
                    'orog':'EP_orog.nc'}
    missingdata_dict = {'pft':2.0E20,   'soil':2.0E20, 'slt':2.0E20,
                     'smow':2.0E20,  'disturb':2.0E20,
                     'frac':2.0E20, 'orog':0.0}

    #filestochg = ['pft','soil','slt','smow','disturb','frac','orog']
    filestochg = ['frac','orog']
    
    for file in filestochg:
        cubelist_out = CubeList([])
        cubes = iris.load(filename_dict.get(file))
        missing = missingdata_dict.get(file)

        for cube in cubes:
            # change the CAS
            latitudes = cube.coord('latitude').points
            longitudes = cube.coord('longitude').points
            data = np.copy(cube.data)
            
            for j,lat in enumerate(latitudes):
                if 7. < lat < 13.:
                    for i,lon in enumerate(longitudes):
                        if 270. < lon < 285.:
                            latlon = '['+np.str(lat) + ',' + np.str(lon) + ']'
                            # change if required otherwise don't
                            changereq = gridbox_missing.get(latlon)
                            if changereq != None:
                                data[:,:,j,i] = missing

            EPcube = cube.copy(data)
            cubelist_out.append(EPcube)

            if file == 'orog':
                iris.save(cubelist_out,fileout_dict.get(file))
            else:
                iris.save(cubelist_out,fileout_dict.get(file),
                          fill_value=missing)


            
    sys.exit(0)

def merge_SST_SAT(cubeSST,cubeSAT):
    """
    we are going to use SAT in regions where there is not an SST
    returns a surface temperature field on the best grid possible
    """  
    # we want the temperature in kelvin so add 273.15 onto SST
    cubeSST.data=cubeSST.data + 273.15
    del cubeSST.attributes["valid_min"]
    del cubeSST.attributes["valid_max"]
    
    # put SAT onto same grid as SST (use linear interpolation because it 
    # will rarely be used

    cubeSATr = cubeSAT.regrid(cubeSST,iris.analysis.Linear())
 
    # if there is some missing data in SST replace with SAT
    newSST=np.where(cubeSST.data > 1.0E10, cubeSATr.data, cubeSST.data)
  
    # if any values are < 271.40 then replace with 271.40
    newSST2 = np.where(newSST < 271.40, 271.40, newSST)

    cubeSST_upd = cubeSST.copy(data=newSST2)
    cubeSST_upd.long_name='combined SATSST'

    
    return cubeSST_upd


def interpolate_SST(cubeSST):
    """
    if we have a missing value fill it in by interpolating SST across continents
    """
    SST_data = np.copy(cubeSST.data)

    # where missing values are in the middle of the grid
    for j in range(0,len(cubeSST.coord('latitude').points)):
    #for j in range(138,139):
        print(j)
        prevval=-999.
        nextval=-999.
        for i in range(0,len(cubeSST.coord('longitude').points)):
            print('i is',i, SST_data[0,0,j,i])
            if SST_data[0,0,j,i] < 1.0E10:
                prevval = SST_data[0,0,j,i]
            elif prevval > -900:  # need to linear interpolate
                count=1
                # find how many you need to linearly interpolate over
                for i2 in range(i+1,len(cubeSST.coord('longitude').points)):
                    print('i2',i2,SST_data[0,0,j,i2])
                    if SST_data[0,0,j,i2] > 1.0E10:
                        count=count+1
                    else:
                        nextval=SST_data[0,0,j,i2]
                        break
                # linearly interpolate over these values
                if prevval > -900 and nextval > -900:
                    for i2 in range(i,i+count):
                        if SST_data[0,0,j,np.mod(i2,288)] > 1.0E10:
                            print('changing',i2,count,SST_data[0,0,j,i2])
                            SST_data[0,0,j,i2] = (prevval + 
                              ((i2-i+1) *(nextval - prevval)/(count+1)))
                i=i+count
                count=0
                prevval=-999
                nextval=-999
                       
            else:
               pass


    # where missing values are on the edge of the grid
    for j in range(0,len(cubeSST.coord('latitude').points)):
    #for j in range(138,139):
        print(j)
        prevval=-999.
        nextval=-999.
        for i in range(140,len(cubeSST.coord('longitude').points+140)):
            print('i is',i, SST_data[0,0,j,i])
            if SST_data[0,0,j,np.mod(i,288)] < 1.0E10:
                prevval = SST_data[0,0,j,i]
            elif prevval > -900:  # need to linear interpolate
                count=1
                # find how many you need to linearly interpolate over
                for i2 in range(i+141,len(cubeSST.coord('longitude').points)+140):
                    print('i2',i2,SST_data[0,0,j,np.mod(i2,288)])
                    if SST_data[0,0,j,np.mod(i2,288)] > 1.0E10:
                        count=count+1
                    else:
                        nextval=SST_data[0,0,j,np.mod(i2,288)]
                        break
                # linearly interpolate over these values
                if prevval > -900 and nextval > -900:
                    for i2 in range(i+140,i+count+140):
                        if SST_data[0,0,j,np.mod(i2,288)] > 1.0E10:
                            SST_data[0,0,j,np.mod(i2,288)] = (prevval + 
                              ((i2-i+1) *(nextval - prevval)/(count+1)))
                i=i+count
                count=0
                prevval=-999
                nextval=-999
                       
            else:
               pass


    #sys.exit(0)
    newcube = cubeSST.copy(data=SST_data)
    print('interpolated')
                     
    return newcube

def extract_HCM3(field):
    """
    extracts the HadCM3 SST and puts into a cube
    """
    cubes = CubeList()
    cubesSST = CubeList()
    cubesSAT = CubeList()
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    filestart = '/nfs/hera1/earjcti/um/' + expt + '/netcdf/' + expt

    for year in startyear_HCM3,endyear_HCM3:
        for month in months:
            if field == 'SST':
                file = filestart + 'o#pf00000' + str(year) + month + '+.nc'
                cubeSST = iris.load_cube(file,'temp')
                cubesSST.append(cubeSST)

                file = filestart + 'a#pd00000' + str(year) + month + '+.nc'
                cubeSAT = iris.load_cube(file,
                                         'SURFACE TEMPERATURE AFTER TIMESTEP')
                cubesSAT.append(cubeSAT)
                
                # one way of getting the temperatures is to merge
                #temperature_cube = merge_SST_SAT(cubeSST,cubeSAT)
                #cubes.append(temperature_cube)
                
                # the other way of getting the temperatures is to linearly
                # interpolate across continents

                alt_temperature_cube = interpolate_SST(cubeSST)
                cubes.append(alt_temperature_cube)

                
    iris.util.equalise_attributes(cubes)
    iris.util.equalise_attributes(cubesSST)
    iris.util.equalise_attributes(cubesSAT)
    for cube in cubes:
        print(cube)
    allcubes = cubes.concatenate_cube()
    SATcubes = cubesSAT.concatenate_cube()
    SSTcubes = cubesSST.concatenate_cube()
    iris.save([allcubes,SATcubes,SSTcubes],'temporary.nc',fill_value=2.0E20)
    sys.exit(0)

    return cube
    
##########################################################################
# setup what to change

expt = 'xpsic'
startyear_HCM3=2899
endyear_HCM3=2900
startyear_HGEM=1850

HCM3_sst_cube = extract_HCM3('SST')
HCM3r_sst_cube = regrid_HCM3(HCM3_sst_cube,'SST')
HGEM_sst_cube = correct_HCM3_cube(HCM3r_sst_cube,'SST')
