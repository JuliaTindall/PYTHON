#NAME
#
#  This program is based on the calculations of the AMOC but will actually
#  calculate ocean heat transport for the Global, Atlantic, Pacific and Indian
#  oceans
#
#  it will write all the OHT out to a file.
#
#

# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import sys
#from netCDF4 import Dataset, MFDataset
#from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



def get_mask(firstfile):
    """
    gets the mask based on the basin file and the land sea mask
    """

    # read the basin file

    filename = 'basin_hadcm3'
    f=open(filename,'r')
    #discard first 3 lines
    content=f.readline()
    content=f.readline()
    content=f.readline()

    Atlantic_mask = np.ma.zeros((20,144,288))
    Pacific_mask = np.ma.zeros((20,144,288))
    Indian_mask = np.ma.zeros((20,144,288))
    Global_mask = np.ma.zeros((20,144,288))
  

    for j in range(144,0,-1):
        content=f.readline()
        # we have 4 basins - each are split into two divisions with a start
        # and end point
        (rowno,
         bas1_d1_s,bas1_d1_e,bas1_d2_s,bas1_d2_e,
         bas2_d1_s,bas2_d1_e,bas2_d2_s,bas2_d2_e,
         bas3_d1_s,bas3_d1_e,bas3_d2_s,bas3_d2_e,
         bas4_d1_s,bas4_d1_e,bas4_d2_s,bas4_d2_e)=content.split()
        
        for i in range(int(bas1_d1_s),int(bas1_d1_e)+1):
            Indian_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0
        for i in range(int(bas1_d2_s),int(bas1_d2_e)+1):
            Indian_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0

        for i in range(int(bas2_d1_s),int(bas2_d1_e)+1):
            Pacific_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0
        for i in range(int(bas2_d2_s),int(bas2_d2_e)+1):
            Pacific_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0

        for i in range(int(bas3_d1_s),int(bas3_d1_e)+1):
            Atlantic_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0
        for i in range(int(bas3_d2_s),int(bas3_d2_e)+1):
            Atlantic_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0
    
        for i in range(int(bas4_d1_s),int(bas4_d1_e)+1):
            Global_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0
        for i in range(int(bas4_d2_s),int(bas4_d2_e)+1):
            Global_mask[:,int(rowno)-1,np.mod(i-1,288)]=1.0

    f.close

    # get LSM from temperature:


    # put the data on a land sea mask
    cube_full = iris.load_cube(firstfile,'TEMPERATURE (OCEAN)  DEG.C')
    cube = cube_full[0,:,:,:]
    cube_data = cube.data
   
    Atlantic_mask.mask = cube.data.mask
    Atlantic_cube = cube.copy(data=Atlantic_mask)
    Pacific_mask.mask = cube.data.mask
    Pacific_cube = cube.copy(data=Pacific_mask)
    Indian_mask.mask = cube.data.mask
    Indian_cube = cube.copy(data=Indian_mask)
    Global_mask.mask = cube.data.mask
    Global_cube = cube.copy(data=Global_mask)

    #plt.subplot(2,2,1)
    #qplt.contourf(Atlantic_cube[0,:,:])
    #plt.subplot(2,2,2)
    #qplt.contourf(Pacific_cube[0,:,:])
    #plt.subplot(2,2,3)
    #qplt.contourf(Indian_cube[0,:,:])
    #plt.subplot(2,2,4)
    #qplt.contourf(Global_cube[0,:,:])
    #plt.show()
    #sys.exit(0)
   
  

    return Atlantic_cube, Pacific_cube, Indian_cube, Global_cube



####################################################################
def get_params(filename):
    """
    gets the spacing
    """

    # thickness
    cube = iris.load_cube(filename,'W')
    depths=np.copy(cube.coord('depth').points)
    depths2=np.zeros(21) # depths used to calculate thickness
    ndepths = len(depths)
    thickness=np.zeros(20)
    
    for k in range(ndepths,0,-1):
        depths2[k]=depths[k-1]
    depths2[0]=0.0
    depths2[20]=depths2[19]+(depths2[19]-depths2[18])

    for k in range(0,20):
        thickness[k]=depths2[k+1]-depths2[k]

    # get latitude and dx from the v field (also save v)
    v_cube = iris.load_cube(filename,'TOTAL OCEAN V-VELOCITY      CM S**-1')
    #v_cube = iris.util.squeeze(v_cube)
    v_cube.data = v_cube.data / 100.  # convert to m/2
    lats = v_cube.coord('latitude').points
    coslats= np.cos(lats * 2. * np.pi / 360.)
    lons=v_cube.coord('longitude').points

    # length of a longitude box at the equator
    dx=(lons[1]-lons[0]) * 111320.
 
    Tcube_origgrid = iris.load_cube(filename,
                    'POTENTIAL TEMPERATURE (OCEAN)  DEG.C')
    Tcube_origgrid = Tcube_origgrid + 273.15 # convert to kelvin
    Tcube = Tcube_origgrid.regrid(v_cube,iris.analysis.Linear())
   
    return (dx,thickness,lats, coslats, v_cube,Tcube,Tcube_origgrid)

###################################################
def V_grid_mask(orig_cube,vgrid_cube):
    """
    put mask onto a v_grid 
    v is staggered in longitude compared to temperature.
    we will only use the mask if all surrounding points are masked

    longitude:
    original (longitude goes from 0, 1.25 etc.
    vgrid    (longitude goes from 0.625, 1.875 etc.

    so we will set vgrid_mask(i) as 1 if orig(i) and orig(i+1) are also 1

    latitude:
    original (latitude goes from -89.375, -88.125  etc
    vgrid    (latitude goes from -88.75, -87,5 etc

    so we will set vgrid_mask[j] to 1 if orig[j] and orig[j+1] are also 1
    
    """
    lon_orig=orig_cube.coord('longitude').points
    nlonorig=len(lon_orig)
    lon_vgrid=vgrid_cube.coord('longitude').points
    nlonv=len(lon_vgrid)
    lat_orig=orig_cube.coord('latitude').points
    nlatorig=len(lat_orig)
    lat_vgrid=vgrid_cube.coord('latitude').points
    nlatv=len(lat_vgrid)

    origdata=orig_cube.data
    vmaskdata=np.ma.masked_array(np.zeros(np.shape(vgrid_cube.data)),
                                 mask=np.zeros(np.shape(vgrid_cube.data)))

    for k in range(0,20):
        for j in range(0,nlatv):
            for i in range(0,nlonv):
                ip1=np.mod(i+1,nlonorig)
                if (origdata[k,j,i] ==1 and origdata[k,j+1,i]==1 and
                    origdata[k,j,ip1]==1 and origdata[k,j+1,ip1]==1):
                    vmaskdata[k,j,i]=1.0
                else: #mask
                    vmaskdata.mask[k,j,i]=1.0

    vmaskcube = vgrid_cube.copy(data=vmaskdata)
                
    
    return vmaskcube

def calc_OHT(V_cube,dx,dz,coslats,T_cube,T_cube_origgrid,basin):
    """
    this calculates OHT using PJV streamfuction code
    """
    #calculate zonal integral m2/s-1 (note everything we dont want to use
    #                                 should be masked
    #we multiply it by the width of the gridbox (ie dx * cos latitude)
    #
    #  Note T_cube is temperature interpolated onto the V grid
    #       T_cube_origgrid is temperauture on the T grid
    #

    density=1025
    cp=4000

    T_mean = np.squeeze(np.average(T_cube.data, axis=3))

    # calculate T_mean_depth manually np.average does not work well with missing
    # values
    nlats = len(T_cube.coord('latitude').points)

    if Tref == 'mean':
        T_mean_depth = np.zeros(nlats)
        for j in range(0,nlats):
            Tacc=0
            depacc=0
            for k in range(0,len(dz)):
                if not T_mean.mask[k,j]:
                    Tacc = Tacc + (T_mean[k,j] * dz[k])
                    depacc = depacc + dz[k]
                    
            if depacc > 0:
                T_mean_depth[j] = Tacc / depacc
            else:
                T_mean_depth[j] = np.nan

        T_mean_depth = T_mean_depth.reshape(1,1,nlats,1)
        Tanom=T_cube.data - T_mean_depth

        
    if Tref == 'zero':
        Tanom=T_cube.data - 273.15


    # here we are going to try upstream sampling
    # if the flow is northwards (+ve) use the temperature from the south
    # if the flow is southwards (-ve) use the temperature from the north
    # use the eram
    #  I have decided that upstream sampling does not help and so have
    #  removed it
    
    #VTdata = np.copy(V_cube.data)
    #Vdata = V_cube.data

    #print(np.shape(T_mean_depth))
                
    #Torig_data = T_cube_origgrid.data
    #for k in range(0,20):
    #    for j in range(0,143):
    #        for i in range(0,288):
    #            if Vdata[0,k,j,i] < 0.0: # southward flow use index j+1
    #                                      # which is the northern gridbox
    #                Tuse = (Torig_data[0,k,j+1,i] +
    #                        Torig_data[0,k,j+1,(i+1)%288]) / 2.0
    #                VTdata[0,k,j,i] = (Vdata[0,k,j,i] *
    #                                   (Tuse - T_mean_depth[0,0,j,0]))
    
    #            else: # northward flow use index j
    #                Tuse = (Torig_data[0,k,j,i] +
    #                          Torig_data[0,k,j,(i+1)%288]) / 2.0
    #                VTdata[0,k,j,i] = (Vdata[0,k,j,i] *
    #                                   (Tuse - T_mean_depth[0,0,j,0]))

                    
    
    
    #ztotalVT=np.sum(V_cube.data * T_cube.data,axis=3) * dx * coslats
    ztotalVT=np.sum(V_cube.data * Tanom,axis=3) * dx * coslats
    #ztotalVT=np.nansum(VTdata,axis=3) * dx * coslats
    ztotalV_cube=(V_cube.collapsed('longitude',iris.analysis.SUM)
                  * dx * coslats)
    ztotalVT_cube = ztotalV_cube.copy(data=ztotalVT)

    
    vdzT_cube = ztotalVT_cube.copy()
    for k in range(0,20):
        vdzT_cube.data[:,k,:]=vdzT_cube.data[:,k,:] * dz[k]
    vdzT_cube.long_name = 'VdzT cube'
    
    # now calculate vertical integral of OHT


    totOHT_cube = vdzT_cube.collapsed('depth_1',iris.analysis.SUM)
    totOHT_cube.data = totOHT_cube.data * density * cp / 1.0E15
    totOHT_cube.long_name = 'total ocean heat transport'



#    if basin == 'Atlantic' or basin == 'Pacific':
#        if expt == 'xqbwe' or expt == 'xqbwg':
#            latreq=15.
#        else:
#            latreq=-30.
#        latitude_constraint = iris.Constraint(latitude=lambda lat: lat > latreq)
#        totOHT_cube = totOHT_cube.extract(latitude_constraint)
    
    #iris.save([totOHT_cube,V_cube,T_cube,vdzT_cube,OHT_cube],'temporary2.nc',fill_value =-99999.)


    return totOHT_cube
   
   
#####################################################################

# gets the basins over which we calculate
filestart = '/uolstore/Research/a/hera1/earjcti/um/'
expt = 'xqbwj'
yearstart=3900
yearend=4000
Tref = 'zero'   # this could be 'zero' which is zerodegC or 'mean' which is
                # the mean temperature for that latitude and basin

firstfile = (filestart+expt+'/pg/'+expt+'o#pg00000'+str(yearstart).zfill(4)+'c1+.nc')
Atlantic_cube, Pacific_cube, Indian_cube, Global_cube =get_mask(firstfile)
OHT_global_cubes = CubeList([])
OHT_atlantic_cubes = CubeList([])
OHT_pacific_cubes = CubeList([])


# get a v_cube for a v_grid
(dx,dz,lats,coslats,v_cube,Tcube,Tcube_origgrid) = get_params(firstfile)

# put masks onto a V-grid
Atlmask_cube = V_grid_mask(Atlantic_cube,v_cube[0,:,:,:])
Globmask_cube = V_grid_mask(Global_cube,v_cube[0,:,:,:])
Pacmask_cube = V_grid_mask(Pacific_cube,v_cube[0,:,:,:])


for year in range(yearstart,yearend):
    filename = (filestart+expt+'/pg/'+expt+'o#pg00000'+
                str(year).zfill(4)+'c1+.nc')

    # get other parameters we need
    (dx,dz,lats,coslats,v_cube,Tcube,Tcube_origgrid) = get_params(filename)
  
    # show V in the Atlantic mask
    V_atl = Atlmask_cube * v_cube # this is on grid 20 * 144 * 288
    V_glob = Globmask_cube * v_cube # 
    V_Pac = Pacmask_cube * v_cube # 
    T_atl = Atlmask_cube * Tcube # this is on grid 20 * 144 * 288
    T_glob = Globmask_cube * Tcube # 
    T_Pac = Pacmask_cube * Tcube # 
    
    # calculate OHT
    OHT_glob_cube = calc_OHT(V_glob,dx,dz,coslats,T_glob,Tcube_origgrid,'Global')
    OHT_Atl_cube = calc_OHT(V_atl,dx,dz,coslats,T_atl,Tcube_origgrid,'Atlantic')
    OHT_Pac_cube = calc_OHT(V_Pac,dx,dz,coslats,T_Pac,Tcube_origgrid,'Pacific')
  
    OHT_global_cubes.append(OHT_glob_cube)
    OHT_atlantic_cubes.append(OHT_Atl_cube)
    OHT_pacific_cubes.append(OHT_Pac_cube)


iris.util.equalise_attributes(OHT_global_cubes)
iris.util.equalise_attributes(OHT_atlantic_cubes)
iris.util.equalise_attributes(OHT_pacific_cubes)


print(OHT_global_cubes[0])
print(OHT_global_cubes[1])
OHT_global_timeseries = OHT_global_cubes.concatenate_cube()
OHT_global_timeseries.long_name = 'Ocean heat transport (global)'
OHT_global_timeseries.units ='PW'
OHT_global_timeseries.data = np.ma.where(OHT_global_timeseries.data > 1.0E20,
                                         -99999., OHT_global_timeseries.data)

OHT_atlantic_timeseries = OHT_atlantic_cubes.concatenate_cube()
OHT_atlantic_timeseries.long_name = 'Ocean heat transport (atlantic)'
OHT_atlantic_timeseries.units ='PW'
OHT_atlantic_timeseries.data = np.ma.where(OHT_atlantic_timeseries.data > 1.0E20,
                                         -99999., OHT_atlantic_timeseries.data)


OHT_pacific_timeseries = OHT_pacific_cubes.concatenate_cube()
OHT_pacific_timeseries.long_name = 'Ocean heat transport (pacific)'
OHT_pacific_timeseries.units ='PW'
OHT_pacific_timeseries.data = np.ma.where(OHT_pacific_timeseries.data > 1.0E20,
                                         -99999., OHT_pacific_timeseries.data)

fileout = (filestart + expt + '/OHT_' + expt + '_' + str(yearstart) + '_' + str(yearend) + 'Tref_is_'+ Tref + '.nc')
iris.save([OHT_global_timeseries,OHT_atlantic_timeseries,OHT_pacific_timeseries],fileout,fill_value = -99999.)
