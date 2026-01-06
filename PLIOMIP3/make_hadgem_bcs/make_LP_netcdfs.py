#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created July 2023 by Julia

This program will do the LP boundary conditions for HadGAM
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
#import sys

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


    

def make_lsm_and_river_route():
    """
    make land sea mask and river reouting
    we need to have previously made the Pliocene fractional lsm
    """

    fraclsm_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/' + 
                                  'PLIOMIP3/P4_maskfrac.nc')

    fraclsm_data = fraclsm_cube.data
    lsm_data = np.where(fraclsm_data > 0.0, 1.0, 0.0)

    pi_mask_cubes = iris.load('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/' + 
                                  '00.0k/qrparm.mask.nc')
    for cube in pi_mask_cubes:
        if cube.var_name == 'lsm':
            lsm_cube = cube.copy(data=lsm_data)
            
        if cube.var_name == 'field700':
            river_cube=cube

    iris.save([lsm_cube,river_cube],'/nfs/hera1/earjcti/UM_ANCIL/HadGEM/' + 
                                  'PLIOMIP3/P4_qrparm.mask.nc')

    

def make_veg():
    """
    because veg is very difficult to translate we will
    1. get from HadCM3
    1.5 regrid for HadGEM and put on HadGEM LSM.
    2. overwrite ice (from PlioMIP2)
    3. find gridboxes where veg fraction does not add up to 1
    4. Any gridboxes where vegfrac > 1.0 scale down (conserving ice)
    4. Any gridboxes where 0.5 < vegfrac < 1.0 scale up (conserving ice)
    4. If under 0.5 use vegetation from nearest gridbox
    """
    #get preindustrial veg cube on correct grid. 
    #Also get HadCM3 veg cube on wrong grid
    #Also get ice on the original pliomip2 grid (1X1)
    preveg_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/qrparm.veg.frac_1860.nc')
    origgrid_plioveg_hcm3_cube = iris.load_cube('/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_mb_qrfrac.type.nc')
    origgrid_plio_icecube = iris.load_cube('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_icemask_v1.0.nc')

    # regrid
    origgrid_plio_icecube.coord('longitude').guess_bounds()
    origgrid_plio_icecube.coord('latitude').guess_bounds()
    origgrid_plioveg_hcm3_cube.coord('longitude').guess_bounds()
    origgrid_plioveg_hcm3_cube.coord('latitude').guess_bounds()
    preveg_cube.coord('longitude').guess_bounds()
    preveg_cube.coord('latitude').guess_bounds()


    plioveg_hgem2_cube = origgrid_plioveg_hcm3_cube.regrid(preveg_cube,iris.analysis.AreaWeighted())
    plio_icecube = origgrid_plio_icecube.regrid(preveg_cube[0,0,:,:],iris.analysis.Nearest())

    # overwrite the land sea mask with the correct one
    plio_lsm=iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/P4_qrparm.mask.nc','LAND MASK (No halo) (LAND=TRUE)')
    plioveg_hgem2_data = plioveg_hgem2_cube.data
    plioveg_hgem2_data.mask = (plio_lsm.data * -1.0) + 1.0

    # 2.0 overwrite vegetation with ice
    plio_ice_data = plio_icecube.data
    print(np.shape(plio_ice_data))
    for j in range(0,144):
        for i in range(0,192):
            if plio_ice_data[j,i] == 2.0:
                plioveg_hgem2_data[0,8,j,i]=1.0
                plioveg_hgem2_data[0,0:8,j,i]=0.0
  
    #2.1  anywhere ice is fractional it is 'leakage' from HadCM3 so get rid
    # check that ice amount is conserved
    plioveg_hgem2_data[0,8,:,:] = np.ma.where(plioveg_hgem2_data[0,8,:,:] < 0.5,
                                             0.0,
                                             plioveg_hgem2_data[0,8,:,:])
    plioveg_hgem2_data[0,8,:,:] = np.ma.where(plioveg_hgem2_data[0,8,:,:] > 0.5,
                                             1.0,
                                             plioveg_hgem2_data[0,8,:,:])
    hadgem2_ice_cube = plioveg_hgem2_cube[0,8,:,:]
    hadgem2_area_ice = hadgem2_ice_cube.collapsed(['latitude','longitude'],
                                                  iris.analysis.MEAN)
    origicecube_area = origgrid_plio_icecube.collapsed(['latitude','longitude'],
                                                       iris.analysis.MEAN)
    print('hadgem2_area_ice',hadgem2_area_ice.data)
    print('origcube_area_ice',origicecube_area.data - 1.0)
                
    # 3. find gridboxes where total does not add up to 1 and fix
    sum_pfts_cube = plioveg_hgem2_cube.collapsed('pseudo',iris.analysis.SUM)
    #a.) if > 1.0 scale everything downwards
    sum_pfts_data=sum_pfts_cube.data
    for j in range(0,144):
        for i in range(0,192):
            if sum_pfts_data[0,j,i] > 1.0:
               # print('gt1 at',j,i,sum_pfts_data[0,j,i],1/sum_pfts_data[0,j,i])
                plioveg_hgem2_data[0,0:8,j,i] = (plioveg_hgem2_data[0,0:8,j,i]
                          / sum_pfts_data[0,j,i])
                sum_pfts_data[0,j,i]=np.sum(plioveg_hgem2_data[0,:,j,i])
                #print('after',np.sum(plioveg_hgem2_data[0,:,j,i]))
  
    #b.) if > 0.5 scale everything upwards
    for j in range(0,144):
        for i in range(0,192):
            if 0.5 < sum_pfts_data[0,j,i] < 1.0:
                #print('0.5<1.0',j,i,sum_pfts_data[0,j,i],1/sum_pfts_data[0,j,i])
                plioveg_hgem2_data[0,0:8,j,i] = (plioveg_hgem2_data[0,0:8,j,i]
                          / sum_pfts_data[0,j,i])
                sum_pfts_data[0,j,i]=np.sum(plioveg_hgem2_data[0,:,j,i])
                #print('after',np.sum(plioveg_hgem2_data[0,:,j,i]))
  
    #c.) data points > 0.0 and < 0.99
    for j in range(0,144):
        for i in range(0,192):
            if 0. < sum_pfts_data[0,j,i] < 0.99:
                print('not 0/1',j,i,sum_pfts_data[0,j,i],
                      plioveg_hgem2_data[0,:,j,i])
                sys.exit(0)
                # doesn't look like there are any.  If there are you need to
                # process
             
    #d) land points which are zero.  Use nearest neighbour
    for j in range(0,144):
        for i in range(0,192):
            if sum_pfts_data[0,j,i] ==0.0:
                foundnear='y'
               # see if we can overwrite with point to left/right/upper/lower
                if 0.99 < sum_pfts_data[0,j,np.mod(i-1,192)] < 1.01:
                    plioveg_hgem2_data[0,:,j,i]= (
                               plioveg_hgem2_data[0,:,j,np.mod(i-1,192)])
                    sum_pfts_data[0,j,i]=sum_pfts_data[0,j,np.mod(i-1,192)]
                elif 0.99 < sum_pfts_data[0,j,np.mod(i+1,192)] < 1.01:
                    plioveg_hgem2_data[0,:,j,i]=(
                            plioveg_hgem2_data[0,:,j,np.mod(i+1,192)])
                    sum_pfts_data[0,j,i]=sum_pfts_data[0,j,np.mod(i+1,192)]
                elif 0.99 < sum_pfts_data[0,j-1,i] < 1.01:
                    plioveg_hgem2_data[0,:,j,i]=plioveg_hgem2_data[0,:,j-1,i]
                    sum_pfts_data[0,j,i]=sum_pfts_data[0,j-1,i]
                elif 0.99 < sum_pfts_data[0,j+1,i] < 1.01:
                    plioveg_hgem2_data[0,:,j,i]=plioveg_hgem2_data[0,:,j+1,i]
                    sum_pfts_data[0,j,i]=sum_pfts_data[0,j+1,i]
                else:
                    foundnear='n'
                    print('none nearby :{',i,j)
                # if we haven't found one nearby use nearest in horizontal 
                # direction
                if foundnear == 'n':
                    i2=1
                    while i2 < 100 and foundnear == 'n':
                        #print(i2)
                        if 0.99 < sum_pfts_data[0,j,np.mod(i+i2,192)] < 1.01:
                             plioveg_hgem2_data[0,:,j,i]=(
                                 plioveg_hgem2_data[0,:,j,np.mod(i+i2,192)])
                             sum_pfts_data[0,j,i]=(
                                 sum_pfts_data[0,j,np.mod(i+i2,192)])
                             foundnear='y'
                        if 0.99 < sum_pfts_data[0,j,np.mod(i-i2,192)] < 1.01:
                             plioveg_hgem2_data[0,:,j,i]=(
                                 plioveg_hgem2_data[0,:,j,np.mod(i-i2,192)])
                             sum_pfts_data[0,j,i]=(
                                 sum_pfts_data[0,j,np.mod(i-i2,192)])
                             foundnear='y'
                        i2=i2+1
                if foundnear == 'n':
                    print('definitely none nearby :{',i,j)
  
  

                      
    # looks good, but just check to see if there is anywhere that is not 1.0
    sum_pfts_cube = plioveg_hgem2_cube.collapsed('pseudo',iris.analysis.SUM)
    sum_pfts_data=sum_pfts_cube.data
    for j in range(0,144):
        for i in range(0,192):
            if sum_pfts_data[0,j,i] < 0.99:
                print('error at end',i,j,sum_pfts_data[0,j,i])
                print('need to write code to fix')
                sys.exit(0)
            if sum_pfts_data[0,j,i] > 1.01:
                if plioveg_hgem2_data[0,8,j,i]==1.0: # land ice
                    plioveg_hgem2_data[0,0:8,j,i]=0.0
                    sum_pfts_data[0,j,i] = np.sum(plioveg_hgem2_data[0,:,j,i])
                else:
                    print('codefix needed')
                    print('error at end',i,j,sum_pfts_data[0,j,i])
                    print(plioveg_hgem2_data[0,:,j,i])
  
    iris.save([plioveg_hgem2_cube,plio_icecube,
               origgrid_plioveg_hcm3_cube,origgrid_plio_icecube,
               sum_pfts_cube],'temporary.nc',
               fill_value = -99999.)
    
    iris.save(plioveg_hgem2_cube,'PlioMIP3.veg.frac.nc',
                   fill_value = -99999.)


    #print(preveg_cube)
    #print(plioveg_hgem2_cube)
    #print(plio_icecube)



def make_orog():
    """
    This is how we will make the orog.
    1. Find the topography anomaly from PlioMIP2
    2. Read in the HadGEM3 preindustrial orography (qrparm.orog.new.nc)
    3. Put the PlioMIP3 topography onto the HadGEM grid
    4. Read in the PlioMIP3 HadGEM land sea mask.  If it is an ocean point set
       the topography to 0 (but Nan for testing)
    5. Look for things that might need correcting
       a) find out where topography looks weird
       b) find out where the topography is massively different from PI
    """
    #1. read in topography from PlioMIP3
    filestart = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'
    fileplio = filestart + '/Plio_enh/Plio_enh/Plio_enh_topo_v1.0.nc'
    filepi = filestart + '/Modern_std/Modern_std/Modern_std_topo_v1.0.nc'

    PlioMIP3_P3_orog_cube = iris.load_cube(fileplio)
    PlioMIP3_PI_orog_cube = iris.load_cube(filepi)
    PlioMIP3_P3_orog_cube.coord('longitude').circular = True
    PlioMIP3_PI_orog_cube.coord('longitude').circular = True
    PlioMIP3_anom_cube = PlioMIP3_P3_orog_cube - PlioMIP3_PI_orog_cube
    

    #2. read in the HadGEM PI orography and put PlioMIP orog anomaly on 
    #   this grid.
    #   we can then create the first version of the HadGEM pliomip3 orog
    #   by adding our 'pliomip3 anomaly onto the preindustrial
    filename = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/qrparm.orog_new.nc'
    HGEM2_PI_allorog_cubes = iris.load(filename)
    for cube in HGEM2_PI_allorog_cubes:
        if cube.var_name == 'ht':
            HGEM2_PI_orog_cube = cube
            
    HGEM2_PI_orog_cube.coord('longitude').guess_bounds()
    HGEM2_PI_orog_cube.coord('latitude').guess_bounds()
    PlioMIP3_anom_cube.coord('longitude').guess_bounds()
    PlioMIP3_anom_cube.coord('latitude').guess_bounds()

    regrid_P3_anom_cube = PlioMIP3_anom_cube.regrid(HGEM2_PI_orog_cube,
                                          iris.analysis.AreaWeighted())
    P3_first_orog_cube = HGEM2_PI_orog_cube + regrid_P3_anom_cube
    P3_first_orog_cube.long_name='first attempt PlioMIP3'
            

    #3.  Add in the land sea mask from PlioMIP HadGEM (previously calculated)
    lsm_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/PLIOMIP3/P4_qrparm.mask.nc','LAND MASK (No halo) (LAND=TRUE)')
    print(lsm_cube.data)
    P3_first_orog_data = P3_first_orog_cube.data
    P3_sec_orog_data = np.ma.where(lsm_cube.data == 0.0, -99999.,
                                  P3_first_orog_data)
    P3_sec_orog_data.mask=(lsm_cube.data -1.0) * 1.0

    print(P3_sec_orog_data)
    print(P3_first_orog_cube.data)
    P3_sec_orog_cube = P3_first_orog_cube.copy(data=P3_sec_orog_data)
    P3_sec_orog_cube.long_name='second attempt PlioMIP3'

    #4.  There are a number or land points below sea level 
    #    if they are on the coast set the orography height to zero

    P3_third_orog_data=np.ma.copy(P3_sec_orog_data)
    for j in range(0,144):
        for i in range(0,192):
            if (P3_third_orog_data[0,0,j,i] < 0.0 and 
                P3_third_orog_data.mask[0,0,j,i]==0.0):
                if (P3_third_orog_data.mask[0,0,j-1,i]==0.0 or
                    P3_third_orog_data.mask[0,0,j+1,i]==0.0 or
                    P3_third_orog_data.mask[0,0,j,np.mod(i,192)]==0.0 or
                    P3_third_orog_data.mask[0,0,j,np.mod(i,192)]==0.0):
                    #print('fixing coast',i,j,P3_third_orog_data[0,0,j,i])
                    P3_third_orog_data[0,0,j,i]=0.0
                else:
                    print('unfixed',i,j,P3_third_orog_data[0,0,j,i])
    P3_third_orog_cube = P3_sec_orog_cube.copy(data=P3_third_orog_data)
    P3_third_orog_cube.long_name='third attempt PlioMIP3'


    #5.Find difference between PlioMIP and PI orography on 
    diff_PI_cube=P3_third_orog_cube - HGEM2_PI_orog_cube 
    diff_PI_cube.long_name='third attempt - pi'

   

    iris.save([PlioMIP3_P3_orog_cube,PlioMIP3_PI_orog_cube,
               PlioMIP3_anom_cube,regrid_P3_anom_cube,
               P3_first_orog_cube,P3_sec_orog_cube,P3_third_orog_cube,
               diff_PI_cube],
              'temporary.nc',fill_value=-99999.)
    P3_third_orog_cube.data = np.where(P3_third_orog_cube.data < -900, 0,
                                       P3_third_orog_cube.data)
    P3_third_orog_cube.long_name='Pliocene orography'


    iris.save(P3_third_orog_cube,'orog_only.nc',fill_value=-99999.)

 

def make_orog_gradients():
    """
    previously need to have run make_orog which produces the file orog_only.nc
    1. read in pliocene and preindustrial orography
    2. Find out if the Pliocene orography is much different to the preindustrial
       anomaly as follows
       If preindustrial orography > 3m find percentage difference between
       pliocene anomaly and preinddustrial anomaly
       If preindustrial orography < 3m set percentage difference to be the 
       absolute difference
       percentage difference to 'unknown'. Write out the percentage difference
       to the tempeorary file
    """
    #read in pliocene orog cube and preindustrial all cubes  
    HGEM2_P3_orog_cube = iris.load_cube('orog_only.nc')
    
    filename = '/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/qrparm.orog_new.nc'
    HGEM2_PI_allorog_cubes = iris.load(filename)
    for cube in HGEM2_PI_allorog_cubes:
        if cube.var_name == 'ht':
            HGEM2_PI_orog_cube = cube
        if cube.var_name == 'field150':
            HGEM2_PI_stdev_orog_cube = cube
      
    # get percentage change of orography (as detailed in description)
    PI_data = HGEM2_PI_orog_cube.data
    Plio_data = HGEM2_P3_orog_cube.data
    anom_data = Plio_data - PI_data
    percent_change_data = np.ma.zeros(np.shape(PI_data))
    for j in range(0,145): # loop over latitudes
        for i in range(0,192): # loop over longitudes
            if PI_data[0,0,j,i] > 3:
                percent_change_data[0,0,j,i]=(anom_data[0,0,j,i] * 100.0 / 
                                              PI_data[0,0,j,i])
            else:
                percent_change_data[0,0,j,i]=anom_data[0,0,j,i]

    percent_change_cube=HGEM2_PI_orog_cube.copy(data=percent_change_data)
    percent_change_cube.long_name='percentage change'

    # find the gridboxes we need to change.  This will be where the 
    # percentage change is < -20 or > 20.  If > 100 change to 100 so we
    # can see detail
    to_change = np.ma.where(np.abs(percent_change_data) > 20.0,
                            percent_change_data, -99999.)
    to_change = np.ma.where(percent_change_data > 100.0, 100.0,
                            percent_change_data)
    to_change.mask = np.ma.where(np.abs(percent_change_data) > 20.0,
                            0.0, 1.0)
    to_change_cube = percent_change_cube.copy(data=to_change)
    to_change_cube.long_name = 'gridboxes where orography is different'


  
    # now set up standard deviation of orography
    #a) if it is within 20% of predindustrial use preindustrial
    PI_stdev_data=HGEM2_PI_stdev_orog_cube.data
    P3_stdev_orog_data = np.ma.where(np.abs(percent_change_data) < 20.0,
                                     PI_stdev_data,-99999)
    P3_stdev_cube=HGEM2_PI_stdev_orog_cube.copy(data=P3_stdev_orog_data)

    iris.save([percent_change_cube,to_change_cube,HGEM2_PI_orog_cube,
               HGEM2_P3_orog_cube,P3_stdev_cube],
              'temporary.nc',fill_value=-99999.)
  
  
    
##########################################################################
# setup what to change

# DO FRACTIONAL LAND SEA MASK
#  THIS WILL NEED DOING BEFORE WE CAN DO ANYTHING ELSE

#make_lsmfrac_LP_1()  # first set of changes to mask
#make_lsmfrac_LP_2()  # second set of changes to mask
#checks_on_mask()  # check the mask looks nice

# IF FRACTIONAL LAND SEA MASK HAS BEEN DONE


#make_lsm_and_river_route()
#make_veg()
#make_orog()
make_orog_gradients() # need to have previously run make_orog
