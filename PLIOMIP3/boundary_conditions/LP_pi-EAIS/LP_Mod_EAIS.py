#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will produce the boundary conditions for the 
LP_Mod_EAIS experiment.

We will do this as follows:
1. start with late Pliocene topography.nu
2. Plot differences in topography between this and preindustrial topography
3. plot where one is land and the other is ocean


"""

import numpy as np
import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
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

    return plio_cube, pi_cube

def get_EAIS_ice_mask(plio_cube, pi_cube):
    """
    inputs: pi and pliocene_core topography
    output: plot showing the differences in topography and ice sheets

    """

    lats = plio_cube.coord('latitude').points
    lons = plio_cube.coord('longitude').points
    newpliodata = np.zeros((len(lats),len(lons)))
   
    plio_data = plio_cube.data
    pi_data = pi_cube.data

    difftopo = plio_cube - pi_cube
    plio_mask_data = np.where(plio_data > 0.0, 1.0, 0.0)
    preind_mask_data = np.where(pi_data > 0.0, 1.0, 0.0)
    diff_mask = plio_mask_data - preind_mask_data
    diff_mask_cube = plio_cube.copy(data=diff_mask)


    # plots to check how the topogrpahy looks different
#    fig = plt.figure(figsize=[12,12])
#    plt.subplot(2,2,1)
#    V = np.arange(-300, 300, 50)
#    cs = qplt.contourf(difftopo, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('difftopo')

#    plt.subplot(2,2,2)
#    cs = iplt.pcolormesh(diff_mask_cube, vmin=-2, vmax=2)
#    plt.colorbar(cs,orientation='horizontal')
#    plt.gca().coastlines()
#    plt.title('difftopo')

#    ax3=plt.subplot(2,2,3,projection=ccrs.SouthPolarStereo())
#    ax3.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
#    cs = qplt.contourf(difftopo, levels=V, cmap='RdBu_r',extend='both')
#    ax3.coastlines()
#    ax3.gridlines()
#    plt.title('difftopo')

#    ax4=plt.subplot(2,2,4,projection=ccrs.SouthPolarStereo())
#    ax4.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
#    cs = iplt.pcolormesh(diff_mask_cube, vmin=-2, vmax=2)
#    plt.colorbar(cs,orientation='horizontal')
#    ax4.coastlines()
#    ax4.gridlines()
#    plt.title('diff lsm')
#    plt.savefig('difftopo.eps')
#    plt.close()

    # read in the two ice masks
    filestart = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'     
    pi_filename = 'Modern_std/Modern_std/Modern_std_mbiome_v1.0.nc'
    plio_filename = 'Plio_enh/Plio_enh/Plio_enh_icemask_v1.0.nc'

    pi_cube = iris.load_cube(filestart + pi_filename,'Modveg_ModLsm')
    pi_data = pi_cube.data
    pi_ice_data = np.ma.where(pi_data == 28, 2.0, 1.0)
    pi_ice_data = np.ma.where(pi_data.mask == True, 0.0, pi_ice_data)
    pi_ice_data.mask = np.where(pi_data.mask == True, 1.0, 0.0)
    pi_ice_cube = pi_cube.copy(data=pi_ice_data)

    fig = plt.figure(figsize=[12,12])
    ax=plt.subplot(2,2,1,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(pi_ice_cube, vmin=-2, vmax=2)
    #plt.colorbar(cs,orientation='horizontal')
    ax.coastlines()
    ax.gridlines()
    plt.title('piice')
   
    plio_ice_cube = iris.load_cube(filestart + plio_filename,'P4_icemask')
    ax=plt.subplot(2,2,2,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(plio_ice_cube, vmin=-2, vmax=2)
    #plt.colorbar(cs,orientation='horizontal')
    ax.coastlines()
    ax.gridlines()
    plt.title('pliocene')

    # now we are going to find the EAIS
    # 1. South of 60degS and land in Pliocene and preindustrial
    EAIS_mask_data = np.copy(pi_ice_data)
    EAIS_mask_data = np.where(plio_mask_data == 1.0, 1.0, 0.0)
    EAIS_mask_data = np.where(preind_mask_data == 1.0, EAIS_mask_data, 0.0)
    EAIS_mask_cube = diff_mask_cube.copy(data=EAIS_mask_data)
    for j, lat in enumerate(EAIS_mask_cube.coord('latitude').points):
        if lat > -60:
            EAIS_mask_cube.data[j,:]=0.0
 
    # everything between -180 and -30 which is not ice in Pliocene is removed
    plio_ice_data = np.copy(plio_ice_cube.data)
    for j, lat in enumerate(lats):
        if lat < -60:
            for i, lon in enumerate(lons):
                if -180 < lon < -30:
                    if 0.5 < plio_ice_data[j,i] < 1.5:
                        EAIS_mask_data[j,i]=-2.0
 
   # everything between -30 and 160 which is ice in Preindustrial is added
    plio_ice_data = np.copy(plio_ice_cube.data)
    for j, lat in enumerate(EAIS_mask_cube.coord('latitude').points):
        if lat < -60:
            for i, lon in enumerate(EAIS_mask_cube.coord('longitude').points):
                if -30 < lon < 160:
                    if pi_ice_data[j,i] > 1.5 and EAIS_mask_data[j,i] !=1.0:
                        EAIS_mask_data[j,i]=2.0


    ax=plt.subplot(2,2,3,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(EAIS_mask_cube, vmin=-2, vmax=2)
    ax.coastlines()
    ax.gridlines()
    plt.title('EAIS')


    # plot EAIS final mask and write to a file
    EAIS_final_mask = np.ma.where(EAIS_mask_data >0.5, 1.0, 0.0)


    # we are just going to highlight the gridboxes that are
    # ice in the Pliocene and not ice in the preindustrial
    # we will set these to 2.0

    for j,lat in enumerate(lats):
        if lat < -60.0:
            for i, lon in enumerate(lons):
                #if lon == 50.5:
                #    print(lon,lat,plio_ice_cube.data[j,i],pi_ice_data[j,i],
                #          pi_ice_data.mask[j,i])
                if (((plio_ice_cube.data[j,i] > 1.5) and
                    (pi_ice_data[j,i] < 1.5)) or 
                    ((plio_ice_cube.data[j,i] > 1.5) and
                     (pi_ice_data.mask[j,i] == True))):
                    print(lat,lon,plio_ice_cube.data[j,i],pi_ice_data[j,i])
                    EAIS_final_mask[j,i] = 2.0
  
    EAIS_final_mask_cube = EAIS_mask_cube.copy(data=EAIS_final_mask)

    
    EAIS_final_mask_cube.long_name = 'EAIS ice mask'
    ax=plt.subplot(2,2,4,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(EAIS_final_mask_cube, vmin=-2, vmax=2)
    ax.coastlines()
    ax.gridlines()
    plt.title('EAIS mask final (2-Pliocene ice, pi nonice')

    plt.savefig('EAIS_ice_mask.eps')

    iris.save(EAIS_final_mask_cube,'EAIS_ice_mask.nc',
              netcdf_format='NETCDF3_CLASSIC')

    
  
def change_EAIS_topography(plio_cube, pi_cube):
    """
    inputs: pi and pliocene_core topography
    output: pliocene topography with EAIS topography back to pi

    """

    EAIS_mask_cube = iris.load_cube('EAIS_ice_mask.nc','EAIS ice mask')

    # new topography is Pliocene if ice mask=0 or pi if ice mask =1 or 2

    new_topo_data = np.where(EAIS_mask_cube.data > 0.5, 
                             pi_cube.data, plio_cube.data)
    new_topo_cube = plio_cube.copy(data=new_topo_data)
    

    # plot this to check
    fig = plt.figure(figsize=[12,12])

    
    # new topography near antarctica
    ax=plt.subplot(3,3,1,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(new_topo_cube, vmin=-2000, vmax=2000)
    plt.colorbar(cs,orientation='horizontal')
    ax.coastlines()
    ax.gridlines()
    plt.title('new topography')

    # difference between new topography and preindustrial
    ax=plt.subplot(3,3,2,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(new_topo_cube - pi_cube, vmin=-1000, vmax=1000,
                         cmap='RdBu_r')
    plt.colorbar(cs,orientation='horizontal')
    ax.coastlines()
    ax.gridlines()
    plt.title('Topo:  newplio - orig preindu')
   
    # difference between new topography and pliocene
    ax=plt.subplot(3,3,3,projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180,180,-90,-60],ccrs.PlateCarree())
    cs = iplt.pcolormesh(new_topo_cube - plio_cube, vmin=-1000, vmax=1000,
                         cmap='RdBu_r')
    plt.colorbar(cs,orientation='horizontal')
    ax.coastlines()
    ax.gridlines()
    plt.title('Topo:  newplio - orig plio')

    # plot the part between west antarctica and east antarctica

    ax = plt.subplot(3,3,4,projection=ccrs.PlateCarree())
    ax.set_extent([-120,-60,-90,-80])
    cs = iplt.pcolormesh(new_topo_cube - plio_cube, vmin=-1000,vmax=1000,
                         cmap='RdBu_r')
    plt.colorbar(cs,orientation='horizontal')
    plt.title('newplio-plio zoomed in')

    ax = plt.subplot(3,3,5,projection=ccrs.PlateCarree())
    ax.set_extent([-120,-60,-90,-80])
    cs = iplt.pcolormesh(new_topo_cube, vmin=0,vmax=3000,cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    plt.title('New Pliocene Topo zoomed in')

    ax = plt.subplot(3,3,6,projection=ccrs.PlateCarree())
    ax.set_extent([-120,-60,-90,-80])
    cs = iplt.pcolormesh(plio_cube, vmin=0,vmax=3000,cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    plt.title('Original Pliocene Topo zoomed in')

    ax = plt.subplot(3,3,7,projection=ccrs.PlateCarree())
    ax.set_extent([-120,-60,-90,-80])
    cs = iplt.pcolormesh(pi_cube, vmin=0,vmax=3000,cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    plt.title('preindustrial Topo zoomed in')

    # adjust topography between west antarctica and East Antarctica to 
    # avoid sharper gradients than in EOI400
    
    lons=new_topo_cube.coord('longitude').points
    lats = new_topo_cube.coord('latitude').points
    topo_data =new_topo_cube.data
    topo_data_adj = np.copy(topo_data)
    plio_data = plio_cube.data
    for j,lat in enumerate(lats):
        if -90 < lat < -80:
            for i, lon in enumerate(lons):
                if -120 < lon < -60:
                    if topo_data[j,i] - plio_data[j,i] > 0.0:
                        topo_data_adj[j,i] = plio_data[j,i]
    topo_cube_adj = new_topo_cube.copy(data=topo_data_adj)

    ax = plt.subplot(3,3,8,projection=ccrs.PlateCarree())
    ax.set_extent([-120,-60,-90,-80])
    cs = iplt.pcolormesh(topo_cube_adj, vmin=0,vmax=3000,cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    plt.title('Adjusted new plio topo zoomed in')

    ax = plt.subplot(3,3,9,projection=ccrs.PlateCarree())
    ax.set_extent([-120,-60,-90,-80])
    cs = iplt.pcolormesh(topo_cube_adj - plio_cube, vmin=-1000,vmax=1000,cmap='RdBu_r')
    plt.colorbar(cs,orientation='horizontal')
    plt.title('Adj_newplio - orig_plio zoomed in')

    plt.savefig('different_topographies.eps')
    plt.savefig('different_topographies.pdf')

    newpliominpi = topo_cube_adj - pi_cube
    newpliominpi.long_name = 'newplio-pi'
    newpliominplio = topo_cube_adj - plio_cube
    newpliominplio.long_name = 'newplio-origplio'
    unadjminnewplio = new_topo_cube - plio_cube
    unadjminnewplio.long_name = 'unadjusted minus newplio'
    iris.save([topo_cube_adj, newpliominpi, newpliominplio,
               unadjminnewplio],'topograpies.nc',
              netcdf_format='NETCDF3_CLASSIC')
    return 


def topography_and_lsm():
    """
    uses the adjusted topography created in change_EAIS_topography
    writes the topography out to a file
    obtains the LSM and writes that out to a file too
    """

    topography_cube=iris.load_cube('topograpies.nc','p4_topo')
    change_attributes(topography_cube,'Topography')
    iris.save(topography_cube,'LP_pi-EAIS_topo_v1.0.nc')

    # convert to LSM and save

    EOI400_enh = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    EOI400_lsm_cube = iris.load_cube(EOI400_enh,'p4_lsm')

    newdata = topography_cube.data
    lsmdata = np.where(newdata > 0.0, 1.0, 0.0)
    lsmcube = EOI400_lsm_cube.copy(data=lsmdata)
    lsmcube = change_attributes(lsmcube,'LSM') 
    del lsmcube.attributes["invalid_units"]

    #compare_with_EOI400(lsmcube)
    print(lsmcube)
    iris.save(lsmcube,'LP_pi-EAIS_LSM_v1.0.nc')

    return topography_cube, lsmcube
  
###############################################
def compare_with_EOI400(lsmcube):
    """
    here we get the new lsm and compare with the original EOI400 one
    and also the pi
    """
    EOI400_enh = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'

    EOI400_cube = iris.load_cube(EOI400_enh,'p4_lsm')

    modern = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc'

    modern_cube = iris.load_cube(modern,'etopo1_lsm')

    lsm_min_eoi400_cube = lsmcube - EOI400_cube
    lsm_min_modern_cube = lsmcube - modern_cube
    
    # plot
    plt.subplot(1,2,1)
    cs=qplt.contourf(lsm_min_eoi400_cube, cmap='bwr',
                     levels=[-1.5, -0.5, 0.5, 1.5])
    #plt.colorbar(cs,orientation='horizontal')
    plt.gca().coastlines()
    plt.title('LP_Mod_EAIS - EOI400')
    plt.subplot(1,2,2)
    cs=qplt.contourf(lsm_min_modern_cube, cmap='bwr', 
                     levels=[-1.5, -0.5, 0.5, 1.5])
    #plt.colorbar(cs,orientation='horizontal')
    plt.gca().coastlines()
    plt.title('LP_Mod_EAIS - modern')
    plt.savefig('LSM_anom.eps')
   
###############################################
def change_attributes(cube,filetype):
    """
    changes the following attributes to the cube
    Netcdf_author
    Email
    Institution
    Code Used
    Netcdf_Creation_date
    Information
    Information2
    """ 
    datetoday = date.today()
    datestr = str(datetoday)
    cube.attributes.update({'Netcdf_author' : 'Tindall, J. C.'})  
    cube.attributes.update({'Email' : 'J.C.Tindall@leeds.ac.uk'})  
    cube.attributes.update({'Netcdf_Creation_Date' : datestr}) 
    cube.attributes.update({'Information' : "This is the PRISM4 paleogeography reconstruction as presented in Dowsett et al (2016), edited to change the East Antarctic Ice Sheet to modern"})
    try:
        del cube.attributes["Code_author"]
    except:
        print('no code author')
    try:
        del cube.attributes["Code Used"]
    except:
        print('no code author')
    try:
        del cube.attributes["Information2"]
    except:
        print('Information')
    print('long',cube.long_name)
    print('var',cube.var_name)
    print('name',cube.name)
    cube.var_name = cube.var_name + '_Mod_EAIS'
    if filetype == 'Lake':
        cube.attributes.update({'Other_Key':'0=ocean, 1000-1100=Land indicator (1000) + Lake_percentage of grid cell: 1000=0%, 1100=100%'})
  
    return cube

##############################################
def get_ice_and_veg(lsmcube):
    """
    gets the ice cube and the vegetation cube for this experiment.
    Note The ice mask is pliocene, except for over Antarctica where
    it is a hybrid of Pliocene and preindustrial
    """

    icefile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_icemask_v1.0.nc')
    plio_icecube = iris.load_cube(icefile,'P4_icemask')

    EAIS_icecube = iris.load_cube('EAIS_ice_mask.nc')

    newicedata = np.ma.copy(plio_icecube.data)
    lats = plio_icecube.coord('latitude').points
    lons = plio_icecube.coord('longitude').points
    
    for j, lat in enumerate(lats):
        if lat < -60.0:  # If antarctica replace with new values
            for i in range(0,len(lons)):
                if EAIS_icecube.data[j,i] == 1.0:
                    newicedata[j,i] = 2.0
                else:
                    newicedata[j,i] = lsmcube.data[j,i]
                    
    newicedata = np.where(newicedata > 1000., 0.0, newicedata)
    newicecube = plio_icecube.copy(data=newicedata)
    change_attributes(newicecube,'Icemask')
    newicecube.attributes["Details"] =  "Ocean=0, Land=1, Ice=2" 

  

    print('icecube',newicecube)
    iris.save(newicecube,'LP_pi-EAIS_icemask_v1.0.nc')


    ################################################
    # vegetation.

    # see if there are changes in the LSM between the pliocene and the this one
    lsm_modeais_cube = iris.load_cube('LP_pi-EAIS_LSM_v1.0.nc')
    pliolsm = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc')
    lsm_plio_cube = iris.load_cube(pliolsm)
    cubediff = lsm_modeais_cube - lsm_plio_cube
              


    vegfile_plio = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_mbiome_v1.0.nc')

    vegcube_plio = iris.load_cube(vegfile_plio,'P3veg_P4Lsm')

    vegfile_pi = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
                  'Modern_std/Modern_std/Modern_std_mbiome_v1.0.nc')
    vegcube_pi = iris.load_cube(vegfile_pi,'Modveg_ModLsm')

    # Use Pliocene vegetation
    # but use preindustrial vegetation where the ice mask >0.5

    # default
    vegnewdata = np.ma.where(EAIS_icecube.data > 0.5, vegcube_pi.data, 
                          vegcube_plio.data)
    vegnewdata.mask = np.where(vegnewdata < 0.5, True, False)
   
    # where lsm_modeais_cube = 0(sea) and lsm_plio_cube=1 (land)
    # set the vegetation to sea / masked
    vegnewdata = np.ma.where(cubediff.data == -1.0, 0.0, vegnewdata)
    vegnewdata.mask = np.where(cubediff.data == -1.0, True, vegnewdata.mask)

   
    # where lsmmodeais_cube = 1(land) and lsm_plio_cube=1 (sea)
    # set to preindustrial
    vegnewdata = np.ma.where(cubediff.data == 1.0, vegcube_pi.data, vegnewdata)
    vegnewdata.mask = np.where(cubediff.data == 1.0, False, vegnewdata.mask)

#    plt.subplot(2,2,1)
#    qplt.contourf(EAIS_icecube.copy(data=vegnewdata.mask))

    # look for Nan and correct
    ny,nx=np.shape(vegnewdata)
    for j in range(0,ny):
        for i in range(0, nx):
            if np.isfinite(vegnewdata[j,i]):
                pass
            else:
                if vegnewdata.mask[j,i] == False:
#                    print('unknown veg code',i,j,vegcube_plio.coord('longitude').points[i],vegcube_plio.coord('latitude').points[j],vegcube_plio.data[j,i],vegcube_plio.data.mask[j,i],vegnewdata[j,i],vegnewdata.mask[j,i],cubediff.data[j,i])
                    # need to try and correct
#                    print(i,j,np.shape(vegcube_plio.data),np.shape(vegnewdata.mask))
#                    print('up above',vegcube_plio.data[j-1,i],vegcube_plio.data.mask[j-1,i],vegnewdata[j-1,i],vegnewdata.mask[j-1,i])
#                    print('down below',vegcube_plio.data[j+1,i],vegcube_plio.data.mask[j+1,i],vegnewdata[j+1,i],vegnewdata.mask[j+1,i])
#                    print('to the left',vegcube_plio.data[j,i-1],vegcube_plio.data.mask[j+1,i],vegnewdata[j,i-1],vegnewdata.mask[j,i-1])
#                    if i < 359:
#                        print('to the right',vegcube_plio.data[j,i+1],vegcube_plio.data.mask[j,i+1],vegnewdata[j,i+1],vegnewdata.mask[j,i+1])
#                    print(' ')
                    # we have found that one of the problems is fiji
                    # set to tropical forest (biome0)
                    if i ==359 and j == 73:
                        vegnewdata[j,i]=1
                    # the other problem is near the Bering Strait
                    # set to 7
                    if j==161:
                        vegnewdata[j,i]=7


#    plt.subplot(2,2,2)
#    qplt.contourf(EAIS_icecube.copy(data=vegnewdata.mask))


    # look to see where the LSM doesn't match and correct.
    # where lsm cube doesn't mach newvegcube
    # remember lsm=1: land 0: ocean

    veg_01mask = np.where(vegnewdata.mask==True, 0.0, 1.0)
    for j in range(0,ny):
        for i in range(0, nx):
            if lsmcube.data[j,i] != veg_01mask[j,i]:
                if lsmcube.data[j,i] == 0.0: # ocean
                    vegnewdata[j,i] = 0.0
                    vegnewdata.mask[j,i] = True
                else:
                    # we are going to set it to the biome at j-1
                    # tests show this is appropriate
                    vegnewdata[j,i] = vegnewdata[j-1,i]
                    vegnewdata.mask[j,i] = vegnewdata.mask[j-1,i]
                    #print('mask diff',j,i,
                    #      vegcube_plio.coord('longitude').points[i],
                    #      vegcube_plio.coord('latitude').points[j],
                    #      lsmcube.data[j,i], 
                    #      veg_01mask[j,i],vegnewdata[j,i])
                    #print('surrounddata',vegnewdata[j-1,i], vegnewdata[j+1,i],
                    #      vegnewdata[j,i-1])
                    #if i < 359:
                    #    print('surrounddata2',vegnewdata[j,i+1])
                    #else:
                    #    print('surrounddata2',vegnewdata[j,0])

 
    vegnewdata = np.where(vegnewdata > 1000., 0.0, vegnewdata)
   
    newvegcube = vegcube_plio.copy(data=vegnewdata)
    change_attributes(newvegcube,'megabiomes')

    newvegcube.attributes.update({'Veg_Key_0':'0 = Ocean'})


    print(newvegcube)
    iris.save(newvegcube,'LP_pi-EAIS_mbiome_v1.0.nc')
 
   
    return newicecube, newvegcube

##############################################
def get_soils_and_lakes(lsmcube):
    """
    gets the soils cube and the lakes cube for this experiment.
    Note The ice mask is pliocene, except for over Antarctica where
    it is a hybrid of Pliocene and preindustrial
    """
    # see if there are changes in the LSM between the pliocene and the this one
    lsm_modeais_cube = iris.load_cube('LP_pi-EAIS_LSM_v1.0.nc')
    pliolsm = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc')
    lsm_plio_cube = iris.load_cube(pliolsm)

    cubediff = lsm_modeais_cube - lsm_plio_cube
              

    # lakes (lake file is 0/mask for sea, 
    # and a number between 1000 and 1100 overland

    lakefile_plio = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_lake_v1.0.nc')

    lakecube_plio = iris.load_cube(lakefile_plio,'P4lake_P4lsm')

    # where lsm_modeais_cube and lsm_plio_cube are same set lake to EOI400
    newdata = np.where(cubediff.data == 0.0, lakecube_plio.data, 5000.)

    # where lsm_modeais_cube = 0(sea) and lsm_plio_cube=1 (land)
    # set the lake to sea / masked
    newdata = np.where(cubediff.data == -1.0, 0.0, newdata)

    # where lsmmodeais_cube = 1(land) and lsm_plio_cube=1 (sea)
    # set to land without lake 1000/not masked
    newdata = np.where(cubediff.data == 1.0, 1000.0, newdata)

    lakecube_new = lakecube_plio.copy(data=newdata)
    change_attributes(lakecube_new,'Lake')
    iris.save(lakecube_new,'LP_pi-EAIS_lake_v1.0.nc')

    ################################################
    # soils

    soilfile_plio = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_soil_v1.0.nc')
    pliosoil_cube = iris.load_cube(soilfile_plio)
 
    vegfile_pi = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
                  'Modern_std/Modern_std/Modern_std_mbiome_v1.0.nc')
    vegcube_pi = iris.load_cube(vegfile_pi)

    icecube = iris.load_cube('EAIS_ice_mask.nc')

   
    # Use Pliocene soils
    # but use preindustrial vegetation where the ice mask >0.5

    soilnewdata = np.where(icecube.data > 0.5, vegcube_pi.data, 
                          pliosoil_cube.data)
    newsoilcube = pliosoil_cube.copy(data=soilnewdata)
    change_attributes(newsoilcube,'soils')

    newsoilcube.attributes.update({'Other_Key':'0 = Ocean, 28 = Ice'})

    print(newsoilcube)
    sys.exit(0)

    iris.save(newsoilcube,'LP_pi-EAIS_soil_v1.0.nc')
   
    return 

def check_lsm(file_to_check, details, mask_ind):
    """
    check all the files created in this section have the same land sea mask
    details is what will be written out to show what section we are in
    mask_ind = 'm' if sea is masked out
    mask_ind = 0 if sea is set to zero
    mask_ind = 't' if sea is < 1 land > 1
    """     
    print(' ')
    print(details) 
    cubeLSM = iris.load_cube('LP_pi-EAIS_LSM_v1.0.nc')  # 0 sea 1 land
    cubecheck = iris.load_cube(file_to_check)
    if mask_ind == 'm':
        filelsm = np.where(cubecheck.data.mask == True, 0.0, 1.0)
    elif mask_ind == '0':
        filelsm = np.where(cubecheck.data == 0, 0.0, 1.0)
    elif mask_ind == 't':
        filelsm = np.where(cubecheck.data <= 0.0, 0.0, 1.0)
    else:
        print('wrong mask indicator')
        sys.exit(0)
      
        
    diff = cubeLSM.data - filelsm
    #print(cubecheck.data)
    #plt.subplot(211)
    #qplt.contourf(cubeLSM)
    #plt.subplot(212)
    #qplt.contourf(cubeLSM.copy(data=filelsm))
    #plt.show()
    #sys.exit(0)

    differences = np.sum(np.where(diff != 0.0, 1.0, 0.0))
    print('number of different gridboxes =' + np.str(differences))
    
    if differences > 0.0:
        if details == "veg":
            print('veg')
        else:
            print('julia',details)
        diffcube = cubeLSM.copy(data=diff)
        for j, lat in enumerate(diffcube.coord('latitude').points):
            for i, lon in enumerate(diffcube.coord('longitude').points):
                if diffcube.data[j,i] != 0.0:
                    print(lat, lon, diffcube.data[j,i], cubeLSM.data[j,i], filelsm[j,i],cubecheck.data[j,i])
        
        
###############################################
# main program


# get Pliocene and preindustrial files

#pliocore_topo_cube, pi_topo_cube = get_files('topo_v1.0',
#                                             'etopo1_topo','p4_topo')

#get_EAIS_ice_mask(pliocore_topo_cube, pi_topo_cube)
#sys.exit(0)

# change EAIS topography must have the EAIS_ice_mask.nc file that was produced 
# by get_EAIS_ice_mask (above)

#change_EAIS_topography(pliocore_topo_cube, pi_topo_cube)
#sys.exit(0)

# get topography and land sea mask and write out to a file
# also check that the land sea mask agrees with other files
# you must have proviously obtained topographies.nc from change_EAIS_topography

topography_cube, lsmcube = topography_and_lsm()

# create ice mask file and vegetation file
icecube, vegcube = get_ice_and_veg(lsmcube)

# create soils file and lakes file
get_soils_and_lakes(lsmcube)

print(topography_cube)
print(lsmcube)
print(icecube)
print(vegcube)
sys.exit(0)

# check all the files have the same land sea mask
check_lsm('LP_pi-EAIS_topo_v1.0.nc','topography', 't')
check_lsm('LP_pi-EAIS_icemask_v1.0.nc','ice', 'm')
check_lsm('LP_pi-EAIS_mbiome_v1.0.nc','veg', 'm')
check_lsm('LP_pi-EAIS_soil_v1.0.nc','soil', 'm')
check_lsm('LP_pi-EAIS_lake_v1.0.nc','lake', 'm')
# note veg might be nan at some location
