#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will produce the boundary conditions for the 
PlioMIP3_basic_Eoi400 experiment.

Called LP_Min_LSM_Change

The topo will be as for P4 (EOI400) but will be modern over
Canadian Archeapelego
Bering Straight
Those shelves near Australia

The ice sheet file will need redoing. It is EOI400, we just need a different LSM
The vegetation file will need redoing.  It is EOI400 but some of the gridboxes will now be sea

"""

import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import netCDF4
from datetime import date
import sys
import cartopy.crs as ccrs

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

def change_bering_strait(plio_cube, pi_cube):
    """
    inputs: pi and pliocene_core topography
    output: pliocene topography with bering straight changed back to pi

    """

    fig = plt.figure(figsize=[12,6])
    lats = plio_cube.coord('latitude').points
    lons = plio_cube.coord('longitude').points
    newpliodata = np.zeros((len(lats),len(lons)))
   
    plio_data = plio_cube.data
    pi_data = pi_cube.data

    count=0
    for j, lat in enumerate(lats):
        for i, lon in enumerate(lons):
            if (50 < lat < 80 and lon < -150
                and (plio_data[j,i] > 0 and pi_data[j,i] < 0.0)):
                count=count+1
                newpliodata[j,i] = pi_data[j,i]
            elif (50 < lat < 80 and lon > 175
                and (plio_data[j,i] > 0 and pi_data[j,i] < 0.0)):
                count=count+1
                newpliodata[j,i] = pi_data[j,i]
            else:
                newpliodata[j,i] = plio_data[j,i]    
    
    print('bering straits points changed',count)
    new_plio_cube = plio_cube.copy(data = newpliodata)

    ax=plt.subplot(1,3,1,projection=ccrs.NorthPolarStereo())
    ax.set_extent([0,360,55,90],ccrs.PlateCarree())
    #plio_reg_cube = plio_cube.extract(iris.Constraint(latitude = lambda cell: 50 < cell < 80, longitude=lambda cell: cell < -150))
    V = np.arange(-5, 6, 1)
#    #cs = iplt.pcolormesh(plio_reg_cube, vmin=-5, vmax=5, cmap='RdBu_r')
    cs = qplt.contourf(plio_cube, levels=V, cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('plio1')
    

    ax=plt.subplot(1,3,2,projection=ccrs.NorthPolarStereo())
    ax.set_extent([0,360,55,90],ccrs.PlateCarree())
#    V = np.arange(-5, 6, 1)  
#    #cs = iplt.pcolormesh(pi_reg_cube, vmin=-5, vmax=5, cmap='RdBu_r')
    cs = qplt.contourf(pi_cube, levels=V, cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('pi1')
 
    ax=plt.subplot(1,3,3,projection=ccrs.NorthPolarStereo())
    ax.set_extent([0,360,55,90],ccrs.PlateCarree())
    V = np.arange(-5, 6, 1)  
#    #cs = iplt.pcolormesh(newplio_reg_cube, vmin=-5, vmax=5, cmap='RdBu_r')
    cs = qplt.contourf(new_plio_cube, levels=V, cmap='RdBu_r',extend='both')
    plt.gca().coastlines()
    plt.title('newplio')
 
    plt.tight_layout()
    #plt.show()
    #sys.exit(0)
   
#    plt.subplot(3,3,4)
#    V = np.arange(-500, 600, 100)
#    #cs = iplt.pcolormesh(plio_reg_cube, vmin=-500, vmax=500, cmap='RdBu_r')
#    cs = qplt.contourf(plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')
    

#    plt.subplot(3,3,5)
#    V = np.arange(-500, 600, 100)
#    #cs = iplt.pcolormesh(pi_reg_cube, vmin=-500, vmax=500, cmap='RdBu_r')
#    cs = qplt.contourf(pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('pi2')
   
#    plt.subplot(3,3,6)
#    V = np.arange(-500, 600, 100)
#    #cs = iplt.pcolormesh(newplio_reg_cube, vmin=-500, vmax=500,cmap='RdBu_r')
#    cs = qplt.contourf(newplio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('newplio')
   
#    plt.subplot(3,3,7)
#    V = np.arange(-50, 60, 10)
#    #cs = iplt.pcolormesh(plio_reg_cube - pi_reg_cube, vmin=-500, vmax=500,cmap#='RdBu_r')
#    cs = qplt.contourf(plio_reg_cube - pi_reg_cube, levels=V, cmap='RdBu_r',ext#end='both')
#    plt.gca().coastlines()
#    plt.title('Plio-PI')
   
   
#    plt.subplot(3,3,8)
#    V = np.arange(-50, 60, 10)
#    #cs = iplt.pcolormesh(newplio_reg_cube - plio_reg_cube,  vmin=-500, vmax=50#0,cmap='RdBu_r')
#    cs = qplt.contourf(newplio_reg_cube - plio_reg_cube, levels=V, cmap='RdBu_r#',extend='both')

#    plt.gca().coastlines()
#    plt.title('newplio - old plio')
  

#    plt.subplot(3,3,9)
#    V = np.arange(-50, 60, 10)
#    #cs = iplt.pcolormesh(newplio_reg_cube - pi_reg_cube,  vmin=-500, vmax=500,#cmap='RdBu_r')
#    cs = qplt.contourf(newplio_reg_cube - pi_reg_cube, levels=V, cmap='RdBu_r',#extend='both')

#    plt.gca().coastlines()
#    plt.title('newplio - pi')
    
#    plt.show()
#    plt.savefig('Bering_Strait_Changes.png')
#    sys.exit(0)
    return new_plio_cube


def change_canadian_arctic(plio_cube, pi_cube):
    """
    inputs: pi and pliocene topography
    output: pliocene topography with CAA changed back to pi

    """

    lats = plio_cube.coord('latitude').points
    lons = plio_cube.coord('longitude').points
    newpliodata = np.zeros((len(lats),len(lons)))
   
    plio_data = plio_cube.data
    pi_data = pi_cube.data

    count=0
    for j, lat in enumerate(lats):
        for i, lon in enumerate(lons):
            if (lat > 75 and  -75 < lon < -60 and plio_data[j,i] > 0.0):
                newpliodata[j,i] = pi_data[j,i]
            else:
                newpliodata[j,i] = plio_data[j,i]  
            
    print('CAA points changed',count)
  
    new_plio_cube = plio_cube.copy(data = newpliodata)

#    plt.subplot(3,3,1)
#    plio_reg_cube = plio_cube.extract(iris.Constraint(latitude = lambda cell: 60 < cell < 90, longitude=lambda cell: -130 < cell < -50))
#    V = np.arange(-5, 6, 1)
#    #cs = iplt.pcolormesh(plio_reg_cube, vmin=-5, vmax=5, cmap='RdBu_r')
#    cs = qplt.contourf(plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio1')
    

#    plt.subplot(3,3,2)
#    pi_reg_cube = pi_cube.extract(iris.Constraint(latitude = lambda cell: 60 < cell < 90, longitude=lambda cell: -130 < cell < -50))
#    V = np.arange(-5, 6, 1)  
#    #cs = iplt.pcolormesh(pi_reg_cube,  vmin=-5, vmax=5, cmap='RdBu_r')
#    cs = qplt.contourf(pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('pi1')
 

#    plt.subplot(3,3,3)
#    newplio_reg_cube = new_plio_cube.extract(iris.Constraint(latitude = lambda cell: 60 < cell < 90, longitude=lambda cell: -130 < cell < -50))
#    V = np.arange(-5, 6, 1)  
#    #cs = iplt.pcolormesh(newplio_reg_cube,  vmin=-5, vmax=5, cmap='RdBu_r')
#    cs = qplt.contourf(newplio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('newplio')
 

#    plt.subplot(3,3,4)
#    V = np.arange(-500, 600, 100)
#    #cs = iplt.pcolormesh(plio_reg_cube,  vmin=-500, vmax=500, cmap='RdBu_r') 
#    cs = qplt.contourf(plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')
    

#    plt.subplot(3,3,5)
#    V = np.arange(-500, 600, 100)
#    #cs = iplt.pcolormesh(pi_reg_cube,  vmin=-500, vmax=500, cmap='RdBu_r')
#    cs = qplt.contourf(pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')

#    plt.gca().coastlines()
#    plt.title('pi2')
   
#    plt.subplot(3,3,6)
#    V = np.arange(-500, 600, 100)
#    #cs = iplt.pcolormesh(newplio_reg_cube,  vmin=-500, vmax=500, cmap='RdBu_r'#)
#    cs = qplt.contourf(newplio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')

#    plt.gca().coastlines()
#    plt.title('newplio')
   
#    plt.subplot(3,3,7)
#    V = np.arange(-50, 60, 10)
#    #cs = iplt.pcolormesh(plio_reg_cube - pi_reg_cube,  vmin=-500, vmax=500, cm#ap='RdBu_r')
#    cs = qplt.contourf(plio_reg_cube-pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')

#    plt.gca().coastlines()
#    plt.title('orig diff')
   
   
#    plt.subplot(3,3,8)
#    V = np.arange(-50, 60, 10)
#    #cs = iplt.pcolormesh(newplio_reg_cube - plio_reg_cube,  vmin=-5, vmax=5, cmap='RdBu_r')
#    cs = qplt.contourf(newplio_reg_cube - plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')

#    plt.gca().coastlines()
#    plt.title('newplio - old plio')
  

#    plt.subplot(3,3,9)
#    V = np.arange(-50, 60, 10)
#    #cs = iplt.pcolormesh(newplio_reg_cube - pi_reg_cube,  vmin=-500, vmax=500, cmap='RdBu_r')
#    cs = qplt.contourf(newplio_reg_cube - pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')

#    plt.gca().coastlines()
#    plt.title('newplio - pi')
    
#    plt.savefig('CAA.png')
    return new_plio_cube


def change_australian_shelves(plio_cube, pi_cube):
    """
    inputs: pi and pliocene topography
    output: pliocene topography with shelves near Australia changed back to pi

    """

    lats = plio_cube.coord('latitude').points
    lons = plio_cube.coord('longitude').points
    newpliodata = np.zeros((len(lats),len(lons)))
   
    plio_data = plio_cube.data
    pi_data = pi_cube.data

    count=0
    for j, lat in enumerate(lats):
        for i, lon in enumerate(lons):
         #   if (-30 < lat < 30 and  90 < lon < 180
         #       and (plio_data[j,i] > 0 and pi_data[j,i] < 0.0)):
            if (-25 < lat < 13 and  90 < lon < 150
                 and (plio_data[j,i] > 0 and pi_data[j,i] <= 0.0)):
                count=count+1
                newpliodata[j,i] = pi_data[j,i]
            else:
                newpliodata[j,i] = plio_data[j,i]    
    

    new_plio_cube = plio_cube.copy(data = newpliodata)

#    plt.subplot(3,3,1)
#    plio_reg_cube = plio_cube.extract(iris.Constraint(latitude = lambda cell: -#30 < cell < 30, longitude=lambda cell: 90 < cell < 180))
#    V = np.arange(-5, 6, 1)
#    cs = qplt.contourf(plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio1')
    

#    plt.subplot(3,3,2)
#    pi_reg_cube = pi_cube.extract(iris.Constraint(latitude = lambda cell: -30 <# cell < 30, longitude=lambda cell: 90 < cell < 180))
#    V = np.arange(-5, 6, 1)  
#    cs = qplt.contourf(pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#   
#    plt.title('pi1')
 

#    plt.subplot(3,3,3)
#    newplio_reg_cube = new_plio_cube.extract(iris.Constraint(latitude = lambda #cell: -30 < cell < 30, longitude=lambda cell: 90 < cell < 180))
#    V = np.arange(-5, 6, 1)  
#    cs = qplt.contourf(newplio_reg_cube, levels=V, cmap='RdBu_r',extend='both')##

#    plt.gca().coastlines()
#    plt.title('newplio')
 

#    plt.subplot(3,3,4)
#    V = np.arange(-500, 600, 100)
#    cs = qplt.contourf(plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio2')
    

#    plt.subplot(3,3,5)
#    V = np.arange(-500, 600, 100)
#    cs = qplt.contourf(pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('pi2')
   
#    plt.subplot(3,3,6)
#    V = np.arange(-500, 600, 100)
#    cs = qplt.contourf(newplio_reg_cube, levels=V, cmap='RdBu_r',extend='both')#   plt.gca().coastlines()
#    plt.title('newplio')
   
#    plt.subplot(3,3,7)
#    V = np.arange(-500, 600, 100)
#    cs = qplt.contourf(plio_reg_cube - pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('orig diff')
   
   
#    plt.subplot(3,3,8)
#    V = np.arange(-50, 60, 10)
#    cs = qplt.contourf(newplio_reg_cube - plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('newplio - old plio')
  

#    plt.subplot(3,3,9)
#    V = np.arange(-500, 600, 100)
#    cs = qplt.contourf(newplio_reg_cube - pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('newplio - pi')
#    plt.savefig('Australian_shelves.png')
    return new_plio_cube

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
    plt.subplot(2,1,1)
    qplt.contourf(lsm_min_eoi400_cube, cmap='bwr', levels=[-1.5,-0.5,0.5,1.5])
    plt.gca().coastlines()
    plt.title('LP_MIN_LSM - EOI400')
    plt.subplot(2,1,2)
    qplt.contourf(lsm_min_modern_cube, cmap='bwr', levels=[-1.5,-0.5,0.5,1.5])
    plt.gca().coastlines()
    plt.title('LP_MIN_LSM - modern')
    plt.savefig('LSM_anom.png')
   
###############################################
def change_attributes(cube):
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
    cube.attributes.update({'Information' : "This is the PRISM4 paleogeography reconstruction as presented in Dowsett et al (2016), with some changes to topography"})
    cube.attributes.update({'Information2' : "The changes are that the Bering Strait, the Canadian Arctic Archipelago and the Sahul and Sunda Shelf all use topogrphy from the Modern_std_topo_v1.0.nc file"})
   
    
    try:
        del cube.attributes["Code Used"]
    except:
        print('no code used')

    try:
        del cube.attributes["Code_author"]
    except:
        print('no code author')

    print('long',cube.long_name)
    print('var',cube.var_name)
    print('name',cube.name)
    cube.var_name = cube.var_name + '_Min_LSM_Chg'
  
    return cube

##############################################
def get_ice_and_veg(lsmcube):
    """
    gets the ice cube and the vegetation cube for this experiment.
    Note this is exactly the same as PlioMIP2 EOI400.  
    However some of the land gridboxes will now be ocean.  
    Therefore we just need to remove some gridboxes
    """

    icefile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_icemask_v1.0.nc')
    icecube = iris.load_cube(icefile,'P4_icemask')
    icedata = icecube.data
    icedata = np.where(lsmcube.data == 0.0, 0.0, icedata)
   
    newicecube = icecube.copy(data=icedata)
    change_attributes(newicecube)
    newicecube.attributes["Details"] =  "Ocean=0, Land=1, Ice=2 - ice is as PRISM4 (Dowsett et al 2016)" 


    #########################
    # vegetation

    vegfile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_mbiome_v1.0.nc')
    vegcube = iris.load_cube(vegfile,'P3veg_P4Lsm')
    vegdata=vegcube.data
    vegdata = np.ma.where(lsmcube.data == 0.0, 0.0, vegdata)
    vegdata.mask = np.where(vegdata == 0.0, 1.0, 0.0)
   
       # look for Nan and correct (copied from LP_Mod_EAIS experiments
    ny,nx=np.shape(vegdata)
    for j in range(0,ny):
        for i in range(0, nx):
            if np.isfinite(vegdata[j,i]):
                pass
            else:
                if vegdata.mask[j,i] == False:
                    #print('unknown veg code',i,j,vegcube.coord('longitude').points[i],vegcube.coord('latitude').points[j],vegcube.data[j,i],vegcube.data.mask[j,i],vegdata[j,i],vegdata.mask[j,i])
                    # need to try and correct
                    #print(i,j,np.shape(vegcube.data),np.shape(vegdata.mask))
                    #print('up above',vegcube.data[j-1,i],vegcube.data.mask[j-1,i],vegdata[j-1,i],vegdata.mask[j-1,i])
                    #print('down below',vegcube.data[j+1,i],vegcube.data.mask[j+1,i],vegdata[j+1,i],vegdata.mask[j+1,i])
                    #print('to the left',vegcube.data[j,i-1],vegcube.data.mask[j+1,i],vegdata[j,i-1],vegdata.mask[j,i-1])
                    if i < 359:
                        print('to the right',vegcube.data[j,i+1],vegcube.data.mask[j,i+1],vegdata[j,i+1],vegdata.mask[j,i+1])
                    #print(' ')
                    # we have found that one of the problems is fiji
                    # set to tropical forest (biome0)
                    if i ==359 and j == 73:
                        vegdata[j,i]=1
                    # the other problem is near the Bering Strait
                    # set to 7 - not a problem for this exp
                    #if j==161:
                    #    vegdata[j,i]=7

    #sys.exit(0)
     # look to see where the LSM doesn't match and correct.
    # where lsm cube doesn't mach newvegcube
    # remember lsm=1: land 0: ocean

    veg_01mask = np.where(vegdata.mask==True, 0.0, 1.0)
    for j in range(0,ny):
        for i in range(0, nx):
            if j==161 and i==359:
                print(lsmcube.data[j,i], veg_01mask[j,i],vegdata[j,i])
            if ((lsmcube.data[j,i] != veg_01mask[j,i]) or
                (lsmcube.data[j,i] == 1.0 and vegdata[j,i] == 0.0)):
            
                print('found one')
                if lsmcube.data[j,i] == 0.0: # ocean
                    print('setting land to ocean')
                    vegndata[j,i] = 0.0
                    vegdata.mask[j,i] = True
                else:
                    if j==161 and i==359: # set to 7 because its nearest
                        vegdata[j,i] = 7
                        vegdata.mask[j,i] = False
                    else:
                        print('need to decide what this biome should be')
                        print(i,j)
                        sys.exit(0)
                    print('mask diff',j,i,
                          vegcube.coord('longitude').points[i],
                          vegcube.coord('latitude').points[j],
                          lsmcube.data[j,i], 
                          veg_01mask[j,i],vegdata[j,i])
                    print('surrounddata',vegdata[j-1,i], vegdata[j+1,i],
                          vegdata[j,i-1])
                    if i < 359:
                        print('surrounddata2',vegdata[j,i+1])
                    else:
                        print('surrounddata2',vegdata[j,0])

   
    vegdata = np.where(vegdata > 1000., 0.0, vegdata)
    newvegcube = vegcube.copy(data=vegdata)
    change_attributes(newvegcube)
    newvegcube.attributes.update({'Veg_Key_0':'0 = Ocean'})

   
    return newicecube, newvegcube

##############################################
def get_soil_and_lake(lsmcube):
    """
    gets the soil cube and the lake cube for this experiment.
    Note this is exactly the same as PlioMIP2 EOI400.  
    However some of the land gridboxes will now be ocean.  
    Therefore we just need to remove some gridboxes
    """

    soilfile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_soil_v1.0.nc')
    soilcube = iris.load_cube(soilfile,'P4soil_P4lsm')
    soildata = soilcube.data
    soildata = np.where(lsmcube.data == 0.0, 0.0, soildata)
    soildata = np.where(soildata > 1000, 0.0, soildata)
   
    newsoilcube = soilcube.copy(data=soildata)
    change_attributes(newsoilcube)
    newsoilcube.attributes.update({'Other_Key':'0 = Ocean, 28 = Ice'})
    
  

    lakefile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_lake_v1.0.nc')
    lakecube = iris.load_cube(lakefile,'P4lake_P4lsm')
    lakedata=lakecube.data
    lakedata = np.where(lsmcube.data == 0.0, 0.0, lakedata)
   
    newlakecube = lakecube.copy(data=lakedata)
    change_attributes(newlakecube)
    newlakecube.attributes.update({'Other_Key':'0=ocean, 1000-1100=Land indicator (1000) + Lake_percentage of grid cell: 1000=0%, 1100=100%'})
   
    return newsoilcube, newlakecube
     
###############################################
# main program


# get Pliocene and preindustrial files

pliocore_topo_cube, pi_topo_cube = get_files('topo_v1.0',
                                             'etopo1_topo','p4_topo')


# change bering strait 
new_plio_topo = change_bering_strait(pliocore_topo_cube, pi_topo_cube)

# change Canadian arctic archipelego
new_plio_topo2 = change_canadian_arctic(new_plio_topo, pi_topo_cube)

# change Australasian shelves
new_plio_topo = change_australian_shelves(new_plio_topo2, pi_topo_cube)
new_plio_topo = change_attributes(new_plio_topo)

# save new topography

iris.save(new_plio_topo,'LP_MinLSM_topo_v1.0.nc')

# convert to LSM and save

EOI400_enh = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
EOI400_lsm_cube = iris.load_cube(EOI400_enh,'p4_lsm')

newdata = new_plio_topo.data
lsmdata = np.where(newdata > 0.0, 1.0, 0.0)
lsmcube = EOI400_lsm_cube.copy(data=lsmdata)
lsmcube = change_attributes(lsmcube) 
del lsmcube.attributes["invalid_units"]
iris.save(lsmcube,'LP_MinLSM_LSM_v1.0.nc')

# check LSM has changed properly by comparing with PlioMIP2_enh and PI
#compare_with_EOI400(lsmcube)

# create ice mask file and vegetation file
icecube, vegcube = get_ice_and_veg(lsmcube)
iris.save(icecube,'LP_MinLSM_icemask_v1.0.nc')
iris.save(vegcube,'LP_MinLSM_mbiome_v1.0.nc')

# create soil file and lake file
soilcube, lakecube = get_soil_and_lake(lsmcube)
iris.save(soilcube,'LP_MinLSM_soil_v1.0.nc')
iris.save(lakecube,'LP_MinLSM_lake_v1.0.nc')

#print(new_plio_topo)
print(lsmcube)
#print(icecube)
#print(vegcube)
#print(soilcube)
#print(lakecube)
