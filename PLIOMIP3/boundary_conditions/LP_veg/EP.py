#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created May 2022 by Julia

This program will produce the boundary conditions for the 
PlioMIP3 Early Pliocene Experiment

We will be changing the land/sea mask and the bathymetry in the
CAS region
The icemask and vegetation mask will also need redoing
The CAS will be now at a depth of 500m.
"""

import numpy as np
import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import netCDF4
from datetime import date
#import cmocean to use topo colormap
import sys


def get_files(fileend, plio_fieldname):
    """
    gets the pliocene topography file
    """

    filestart = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'     
    plio_filename = 'Plio_enh/Plio_enh/Plio_enh_' + fileend + '.nc'

    plio_cube = iris.load_cube(filestart + plio_filename,plio_fieldname)

    return plio_cube

def change_cas(plio_cube):
    """
    inputs: pliocene_core topography
    output: pliocene topography with cas opened to 500m

    """

    cubelats = plio_cube.coord('latitude').points
    cubelons = plio_cube.coord('longitude').points
    cubedata = plio_cube.data

    cubedatanew=np.copy(cubedata)
  
    gridboxes_changed=0
    for i in range(0,len(cubelons)):
        if -85.0 < cubelons[i] < -79.0:
            for j in range(0,len(cubelats)):
                if 8.0 < cubelats[j] < 12.0:
                    if cubedatanew[j,i] > 0.0:
                        print('changing cube data',cubelons[i],cubelats[j])
                        # default depth
                        cubedatanew[j,i]=-500.0
                        
                        # -84.5 10.5
                        # average of gridboxes at either side are -203,-304,-250
                        # will set them all to -250
                        if -85.0 < cubelons[i] < -83.0:
                            cubedatanew[j,i]= -250.0
                                                        
                            print(cubelons[i],cubelats[j],cubedatanew[j,i],
                                  cubedatanew[j,i-1],cubedatanew[j,i+1])
                        else:
                            # check gridboxes to the south and north
                            print(cubelons[i],cubelats[j],cubedatanew[j,i],
                                  cubedatanew[j-1,i],cubedatanew[j+1,i])
                            


                        gridboxes_changed = gridboxes_changed+1
    print('gridboxes changed: ', gridboxes_changed)
               
    newcube = plio_cube.copy(data=cubedatanew)

    fig = plt.figure(figsize=[12,12])

  
    ax1=plt.subplot(2,2,1,projection=ccrs.PlateCarree())
    plio_reg_cube = plio_cube.extract(iris.Constraint(latitude = lambda cell: 0< cell < 20, longitude=lambda cell: -100 < cell < -70))
    cs = iplt.pcolormesh(plio_reg_cube, vmin=-5, vmax=5, cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    plt.gca().coastlines()
    gl = ax1.gridlines(draw_labels=True,dms=True,color='black')
    plt.title('plio 5-5')

    ax2=plt.subplot(2,2,2,projection=ccrs.PlateCarree())
    cs = iplt.pcolormesh(plio_reg_cube, vmin=-700, vmax=100, cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    gl = ax2.gridlines(draw_labels=True,dms=True,color='black')
    plt.gca().coastlines()
    plt.title('plio colormesh')

    ax3=plt.subplot(2,2,3,projection=ccrs.PlateCarree())
    plio_reg2_cube = newcube.extract(iris.Constraint(latitude = lambda cell: 0< cell < 20, longitude=lambda cell: -100 < cell < -70))
    cs = iplt.pcolormesh(plio_reg2_cube, vmin=-5, vmax=5, cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    plt.gca().coastlines()
    gl = ax3.gridlines(draw_labels=True,dms=True,color='black')
    plt.title('newcube 5-5')
 

    ax4=plt.subplot(2,2,4,projection=ccrs.PlateCarree())
    cs = iplt.pcolormesh(plio_reg2_cube, vmin=-700, vmax=100, cmap='Accent')
    plt.colorbar(cs,orientation='horizontal')
    plt.gca().coastlines()
    gl = ax4.gridlines(draw_labels=True,dms=True,color='black')
    plt.title('newcube')

    #plt.tight_layout()
    #plt.show()
    #sys.exit(0)
    return newcube



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
    
    print(np.sum(lsm_min_eoi400_cube.data),'is how many gridboxes have changed')
    # plot
    plt.subplot(2,1,1)
    qplt.contourf(lsm_min_eoi400_cube, cmap='bwr', levels=[-1.5,-0.5,0.5,1.5])
    plt.gca().coastlines()
    plt.title('EP - EOI400')
    plt.subplot(2,1,2)
    qplt.contourf(lsm_min_modern_cube, cmap='bwr', levels=[-1.5,-0.5,0.5,1.5])
    plt.gca().coastlines()
    plt.title('EP - modern')
    plt.show()
    sys.exit(0)
    plt.savefig('EP_anom.png')
   
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
    cube.attributes.update({"Code Used" : codeused + " EP.py"})
    try:
        del cube.attributes["Code_author"]
    except:
        print('no code author')
    print('long',cube.long_name)
    print('var',cube.var_name)
    print('name',cube.name)
    cube.var_name = cube.var_name + '_Open_CAS'
  
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


    #======================================================
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
                   # print('unknown veg code',i,j,vegcube.coord('longitude').points[i],vegcube.coord('latitude').points[j],vegcube.data[j,i],vegcube.data.mask[j,i],vegdata[j,i],vegdata.mask[j,i])
                    # need to try and correct
                   # print(i,j,np.shape(vegcube.data),np.shape(vegdata.mask))
                   # print('up above',vegcube.data[j-1,i],vegcube.data.mask[j-1,i],vegdata[j-1,i],vegdata.mask[j-1,i])
                   # print('down below',vegcube.data[j+1,i],vegcube.data.mask[j+1,i],vegdata[j+1,i],vegdata.mask[j+1,i])
                   # print('to the left',vegcube.data[j,i-1],vegcube.data.mask[j+1,i],vegdata[j,i-1],vegdata.mask[j,i-1])
                    #if i < 359:
                    #    print('to the right',vegcube.data[j,i+1],vegcube.data.mask[j,i+1],vegdata[j,i+1],vegdata.mask[j,i+1])
                    #print(' ')
                    # we have found that one of the problems is fiji
                    # set to tropical forest (biome0)
                    if i ==359 and j == 73:
                        vegdata[j,i]=1
                    # the other problem is near the Bering Strait
                    # set to 7
                    if j==161:
                        vegdata[j,i]=7

     # look to see where the LSM doesn't match and correct.
    # where lsm cube doesn't mach newvegcube
    # remember lsm=1: land 0: ocean

    veg_01mask = np.where(vegdata.mask==True, 0.0, 1.0)
    for j in range(0,ny):
        for i in range(0, nx):
            if j==151 and i==10:
                print(lsmcube.data[j,i], veg_01mask[j,i],vegdata[j,i])
            if ((lsmcube.data[j,i] != veg_01mask[j,i]) or
                (lsmcube.data[j,i] == 1.0 and vegdata[j,i] == 0.0)):
            
               # print('found one')
                if lsmcube.data[j,i] == 0.0: # ocean
                    print('setting land to ocean')
                    vegndata[j,i] = 0.0
                    vegdata.mask[j,i] = True
                else:
                    # we are going to set it to the biome at j-1
                    # tests show this is appropriate
                    vegdata[j,i] = vegdata[j-1,i]
                    vegdata.mask[j,i] = vegdata.mask[j-1,i]
               #     print('mask diff',j,i,
               #           vegcube.coord('longitude').points[i],
               #           vegcube.coord('latitude').points[j],
               #           lsmcube.data[j,i], 
               #           veg_01mask[j,i],vegdata[j,i])
               #     print('surrounddata',vegdata[j-1,i], vegdata[j+1,i],
               #           vegdata[j,i-1])
               #     if i < 359:
               #         print('surrounddata2',vegdata[j,i+1])
               #     else:
               #         print('surrounddata2',vegdata[j,0])

   
   
    newvegcube = vegcube.copy(data=vegdata)
    change_attributes(newvegcube)
   
    return newicecube, newvegcube

##############################################    
def get_soils_and_lakes(lsmcube): 
    """
    gets the soils cube and the lakes cube for this experiment.
    Note this is exactly the same as PlioMIP2 EOI400.  
    However some of the land gridboxes will now be ocean.  
    Therefore we just need to remove some gridboxes
    """

    soilfile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
                'Plio_enh/Plio_enh/Plio_enh_soil_v1.0.nc')
    print(iris.load(soilfile))
    soilcube = iris.load_cube(soilfile,'P4soil_P4lsm')
    soildata = soilcube.data
    soildata = np.ma.where(lsmcube.data == 0.0, 0.0, soildata)
    soildata.mask = np.where(soildata == 0.0, 1.0, 0.0)
   
    newsoilcube = soilcube.copy(data=soildata)
    change_attributes(newsoilcube)
    

    lakefile = ('/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/' +
               'Plio_enh/Plio_enh/Plio_enh_lake_v1.0.nc')
    lakecube = iris.load_cube(lakefile,'P4lake_P4lsm')
    lakedata=lakecube.data
    lakedata = np.ma.where(lsmcube.data == 0.0, 0.0, lakedata)
    lakedata.mask = np.where(lakedata == 0.0, 1.0, 0.0)
   
    newlakecube = lakecube.copy(data=lakedata)
    change_attributes(newlakecube)
   
    print(newsoilcube)
    return newsoilcube, newlakecube

###############################################
# main program


# get Pliocene EOI400 experiment

pliocore_topo_cube = get_files('topo_v1.0','p4_topo')


# change cas strait 
new_plio_topo = change_cas(pliocore_topo_cube)
new_plio_topo = change_attributes(new_plio_topo)
new_plio_topo.attributes.update({'Information2' : 'Altered for PlioMIP2 then opened CAS'})  

# save new topography

#iris.save(new_plio_topo,'EP_topo_v1.0.nc')

# convert to LSM and save

EOI400_enh = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
EOI400_lsm_cube = iris.load_cube(EOI400_enh,'p4_lsm')

newdata = new_plio_topo.data
lsmdata = np.where(newdata > 0.0, 1.0, 0.0)
lsmcube = EOI400_lsm_cube.copy(data=lsmdata)
lsmcube = change_attributes(lsmcube)
lsmcube.attributes.update({'Information2' : 'Altered for PlioMIP2 then opened CAS'})  
    
#iris.save(lsmcube,'EP_LSM_v1.0.nc')

# check LSM has changed properly by comparing with PlioMIP2_enh and PI
#compare_with_EOI400(lsmcube)

# create ice mask file and vegetation file
icecube, vegcube = get_ice_and_veg(lsmcube)
#iris.save(icecube,'EP_icemask_v1.0.nc',
#          fill_value=0.0)
iris.save(vegcube,'EP_mbiome_v1.0.nc',
          fill_value=0.0)


# create soils and lakes files
soilscube, lakescube = get_soils_and_lakes(lsmcube)
#iris.save(soilscube,'EP_soil_v1.0.nc',
#          fill_value=0.0)
#iris.save(lakescube,'EP_lake_v1.0.nc',
#          fill_value=0.0)

print(new_plio_topo)
print(lsmcube)
print(icecube)
print(vegcube)
print(soilscube)
print(lakescube)
