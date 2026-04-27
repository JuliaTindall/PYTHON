#NAME
#    AMOC_density
#PURPOSE 
#
#  This program is based on AMOC_based_on_PJV
#  however it will use density as the z_coordinate
#  it will also do the PMOC the GMOC and the IMOC.  It will write them
#  all out to a file

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



def get_mask(filename,basin):
    """
    gets the mask from the Atlantic mask file on a v grid
    """

    cube = iris.load_cube(filename,basin + ' mask V_grid')

    return cube



####################################################################
def get_params(filename,mask_cube):
    """
    gets the spacing
    """

    # thickness
    thick_cube = iris.load_cube(filename,'thickness (dz)')
    
    # get latitude and dx from the v field (also save v)
    #v_cube = iris.load_cube(filename,'V on sigma coordinates')
    v_cube = iris.load_cube(filename,'alternative V cube')
    v_cube = iris.util.squeeze(v_cube)
    v_cube.data = v_cube.data 
    lats = v_cube.coord('latitude').points
    coslats= np.cos(lats * 2. * np.pi / 360.)
    lons=v_cube.coord('longitude').points

    # length of a longitude box at the equator
    dx=(lons[1]-lons[0]) * 111320.

    # mask V cube
    sigma_coord = v_cube.coord('sigma').points
    for k in range(0,len(sigma_coord)):
        v_cube.data[k,:,:]=v_cube.data[k,:,:] * mask_cube.data
        thick_cube.data[k,:,:]=thick_cube.data[k,:,:] * mask_cube.data

    return (dx,thick_cube,lats, coslats, v_cube)


def calc_stream(V_cube,dx,thickness_cube,coslats,basin,exptname):
    """
    this calculates the streamfunction in the same way that PJV did
    """
    #calculate zonal integral m2/s-1 (note everything we dont want to use
    #                                 should be masked
    #we multiply it by the width of the gridbox (ie dx * cos latitude)

    nz=len(V_cube.coord('sigma').points)
    ny=len(V_cube.coord('latitude').points)
    vdz_cube=V_cube.copy(data=V_cube.data * thickness_cube.data)
    vdz_cube.long_name = 'vdz_cube'


    vdz_cube.data.mask = np.where(vdz_cube.data < -999,1.0,0.0)
    tot_thick_cube = thickness_cube.collapsed('sigma',iris.analysis.SUM)
    tot_thick_cube.long_name = 'total thickness'
    
    iris.util.promote_aux_coord_to_dim_coord(vdz_cube, 'sigma')

    ztotalV=np.sum(vdz_cube.data,axis=2) * dx * coslats
    ztotalV_cube=(vdz_cube.collapsed('longitude',iris.analysis.SUM)
                  * dx * coslats)

    
    # now calculate vertical integral (this is the streamfunction).
    # note that if streamfunction is zero it means that at this level
    # everything that has gone north has also gone south and we have a closed
    # circuit

    phi_data=np.ma.zeros(np.shape(ztotalV_cube.data))
    phi_data[0,:]=0.0
    for k in range(1,nz):
        for j in range(0,ny):
            if ztotalV_cube.data.mask[k-1,j]:
                pass
            else:
                phi_data[k,j]=(ztotalV_cube.data[k-1,j]+phi_data[k-1,j])
                if lats[j] == 17.5:
                    print('forward',
                          k,
                          ztotalV_cube.data[k-1,j]/1.0E6,phi_data[k-1,j] / 1.0E6)
                   
    #phi_data.mask = np.where(ztotalV_cube.data.mask == 1.0, 1.0, 0.0)
    phi_cube=ztotalV_cube.copy(data=phi_data / 1.0E6)

  
    # try calculating in reverse
    phi_data_rev=np.zeros(np.shape(ztotalV_cube.data))
    for k in range(nz-1,-1,-1):
        if k == nz-1:
            phi_data_rev[k,:]= (-1.0) *  ztotalV_cube.data[k,:]
        else:
            phi_data_rev[k,:]=(phi_data_rev[k+1,:] - ztotalV_cube.data[k,:])
    phi_rev_cube = phi_cube.copy(data=phi_data_rev/1.0E6)
    phi_rev_cube.long_name = 'reversed'
    phi_rev_cube.data = np.ma.where(phi_rev_cube.data > 1E20,-99999.,
                                    phi_rev_cube.data)

    if basin == 'Atlantic' or basin == 'Pacific':
        xmin=-30
    else:
        xmin=-90

    phi_cube.long_name = 'Density MOC : ' + basin
    phi_cube.data = np.ma.where(phi_cube.data > 1E20,-99999.,phi_cube.data)
    iris.save([phi_cube,phi_rev_cube],'temporary.nc',fill_value = -99999.)
    sys.exit(0)
   
#####################################################################

# file where V and dz are stored on density coordinates
exptname = 'xqbwg'
startyear=3970
endyear=4000
basins = ['Atlantic','Pacific','Global']

filestart = '/home/earjcti/um/'+exptname + '/Vdz/' + exptname + '_vdz_'
for year in range(startyear,endyear):
    filename_density = filestart + str(year) +'.nc'
    for basin in basins:

        # gets the basins over which we calculate
        basinmask_cube=get_mask('masks.nc',basin) # get mask on V grid
        #print('got mask')
     
        # get other parameters we need
        (dx,thickness_cube,lats,
         coslats,V_cube) = get_params(filename_density,basinmask_cube)
        #print('got other params')
       
        # calculate stream function for the basin
        calc_stream(V_cube,dx,thickness_cube,coslats,basin,exptname)
        sys.exit(0)


plt.subplot(2,2,1)
vals=np.arange(-0.1,0.12,0.02)
qplt.contourf(V_atl[0,:,:],levels=vals,cmap='RdBu_r')
plt.subplot(2,2,2)
qplt.contourf(V_atl[5,:,:],levels=vals,cmap='RdBu_r')
plt.subplot(2,2,3)
qplt.contourf(V_atl[10,:,:],levels=vals,cmap='RdBu_r')
plt.subplot(2,2,4)
qplt.contourf(V_atl[15,:,:],levels=vals,cmap='RdBu_r')
plt.show()
sys.exit(0)

qplt.contourf(Atlantic_cube)
plt.title('Atlantic')
plt.show()
plt.close()
