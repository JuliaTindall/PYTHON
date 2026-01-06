#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will create the boundary conditions for the 
PIi560 experiment.

This is a preindustrial experiment with Pliocene ice

We are going to change West Antarctica to Ocean
If ice in the preindustrial or pliocene
a) change topography to Pliocene
"""


THIS PROGRAM IS UNFINISHED

import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import netCDF4
from datetime import date
import sys


def get_files(fileendtopo, pi_toponame, plio_toponame,
              fileendveg, pi_vegname, plio_vegname):
    """
    gets the pliocene and the preindustrial topography and ice
    """

    # get topography
    filestart = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'     
    pi_filename = 'Modern_std/Modern_std/Modern_std_' + fileendtopo + '.nc'
    plio_filename = 'Plio_enh/Plio_enh/Plio_enh_' + fileendtopo + '.nc'

    pi_topocube = iris.load_cube(filestart + pi_filename,pi_toponame)
    plio_topocube = iris.load_cube(filestart + plio_filename,plio_toponame)

    # get vegetation
    pi_filename = 'Modern_std/Modern_std/Modern_std_' + fileendveg + '.nc'
    plio_filename = 'Plio_enh/Plio_enh/Plio_enh_' + fileendveg + '.nc'

    pi_vegcube = iris.load_cube(filestart + pi_filename, pi_vegname)
    plio_vegcube = iris.load_cube(filestart + plio_filename,plio_vegname)
   
    return plio_topocube, pi_topocube,  plio_vegcube, pi_vegcube

def get_ice(vegcube):
    """
    geta an icemask out of a vegetation file
    """
    vegdata = vegcube.data
    ice_data = np.ma.where(vegdata == 28, 1.0, 0.0)
    ice_data = np.ma.where(vegdata.mask == True, 0.0, ice_data)
    ice_data.mask = np.where(vegdata.mask == True, 1.0, 0.0)
    ice_cube = vegcube.copy(data=ice_data)
    ice_cube.long_name='ice'


    return ice_cube


def change_topo(plio_cube, pi_cube, plio_icecube, pi_icecube):
    """
    inputs: pi and pliocene_core topography and ice mask
    output: newpliocene topography which is
               if pliocene or preindustrial is ice: pliocene topography
               otherwise: preindustrial topography 

    """
    
    icemask = np.where(plio_icecube.data==1.0, 1.0, 0.0)
    icemask = np.where(pi_icecube.data==1.0, 1.0, icemask)
    icecube = plio_icecube.copy(data=icemask)
    fig = plt.figure(figsize=[12.0,12.0])
    plt.subplot(2,2,1)
    qplt.contourf(icecube)
    plt.title('a) ice in pliocene or preindustrial')
    plt.subplot(2,2,2)
    qplt.contourf(icecube -pi_icecube)
    plt.title('b) a-preindustrial ice')
    plt.subplot(2,2,3)
    qplt.contourf(plio_icecube)
    plt.title('c) a-preindustrial ice')
    plt.savefig('iceconfig.pdf')
    plt.show()
    sys.exit(0)


GOT TO HERE  THE REST OF THIS SECTION IS FROM CHANGING THE BERING STRAIT
    fig = plt.figure(figsize=[12,12])
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
            else:
                newpliodata[j,i] = plio_data[j,i]    
    

    new_plio_cube = plio_cube.copy(data = newpliodata)

#    plt.subplot(3,3,1)
#    plio_reg_cube = plio_cube.extract(iris.Constraint(latitude = lambda cell: 50 < cell < 80, longitude=lambda cell: cell < -150))
#    V = np.arange(-5, 6, 1)
#    #cs = iplt.pcolormesh(plio_reg_cube, vmin=-5, vmax=5, cmap='RdBu_r')
#    cs = qplt.contourf(plio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('plio1')
    

#    plt.subplot(3,3,2)
#    pi_reg_cube = pi_cube.extract(iris.Constraint(latitude = lambda cell: 50 < cell < 80, longitude=lambda cell: cell < -150))
#    V = np.arange(-5, 6, 1)  
#    #cs = iplt.pcolormesh(pi_reg_cube, vmin=-5, vmax=5, cmap='RdBu_r')
#    cs = qplt.contourf(pi_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('pi1')
 

#    plt.subplot(3,3,3)
#    newplio_reg_cube = new_plio_cube.extract(iris.Constraint(latitude = lambda cell: 50 < cell < 80, longitude=lambda cell: cell < -150))
#    V = np.arange(-5, 6, 1)  
#    #cs = iplt.pcolormesh(newplio_reg_cube, vmin=-5, vmax=5, cmap='RdBu_r')
#    cs = qplt.contourf(newplio_reg_cube, levels=V, cmap='RdBu_r',extend='both')
#    plt.gca().coastlines()
#    plt.title('newplio')
 

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
    
#    plt.savefig('Bering_Strait_Changes.png')
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
    codeused = cube.attributes.get("Code Used")
    cube.attributes.update({"Code Used" : codeused + " LP_MIN_LSM_Change.py"})
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
    icedata = np.ma.where(lsmcube.data == 0.0, 0.0, icedata)
    icedata.mask = np.where(icedata == 0.0, 1.0, 0.0)
   
    newicecube = icecube.copy(data=icedata)
    change_attributes(newicecube)
    newicecube.attributes["Details"] = icecube.attributes.get("Details")


    vegfile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_mbiome_v1.0.nc')
    vegcube = iris.load_cube(vegfile,'P3veg_P4Lsm')
    vegdata=vegcube.data
    vegdata = np.ma.where(lsmcube.data == 0.0, 0.0, vegdata)
    vegdata.mask = np.where(vegdata == 0.0, 1.0, 0.0)
   
    newvegcube = vegcube.copy(data=vegdata)
    change_attributes(newvegcube)
   
    return newicecube, newvegcube
     
###############################################
# main program


# get Pliocene and preindustrial files

(pliocore_topo_cube, pi_topo_cube,
 pliocore_veg_cube, pi_veg_cube) = get_files('topo_v1.0',
                                             'etopo1_topo','p4_topo',
                                             'mbiome_v1.0','Modveg_ModLsm',
                                             'P3veg_P4Lsm')

plio_ice = get_ice(pliocore_veg_cube)
pi_ice = get_ice(pi_veg_cube)

newtopo = change_topo(pliocore_topo_cube, pi_topo_cube,plio_ice,pi_ice)
#START HERE START HERE!!!!!!

# save new topography

iris.save(new_plio_topo,'LP_MinLSM_topo_v1.0.nc',netcdf_format='NETCDF3_CLASSIC')

# convert to LSM and save

EOI400_enh = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
EOI400_lsm_cube = iris.load_cube(EOI400_enh,'p4_lsm')

newdata = new_plio_topo.data
lsmdata = np.where(newdata > 0.0, 1.0, 0.0)
lsmcube = EOI400_lsm_cube.copy(data=lsmdata)
lsmcube = change_attributes(lsmcube) 
iris.save(lsmcube,'LP_MinLSM_LSM_v1.0.nc', netcdf_format='NETCDF3_CLASSIC')

# check LSM has changed properly by comparing with PlioMIP2_enh and PI
#compare_with_EOI400(lsmcube)

# create ice mask file and vegetation file
icecube, vegcube = get_ice_and_veg(lsmcube)
iris.save(icecube,'LP_MinLSM_icemask_v1.0.nc',
          fill_value=0.0, netcdf_format='NETCDF3_CLASSIC')
iris.save(vegcube,'LP_MinLSM_mbiome_v1.0.nc',
          fill_value=0.0, netcdf_format='NETCDF3_CLASSIC')

print(new_plio_topo)
print(lsmcube)
print(icecube)
print(vegcube)
