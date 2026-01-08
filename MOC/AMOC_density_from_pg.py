#NAME
#    AMOC_density_from_pg
#PURPOSE 
#
#  This program will use a pg file.
#  It will regrid V and dz onto density coordinates
#  and will calculate the AMOC/PMOC/IMOC/GMOC
#
#  If the input is a file called:
#    xqbwco#pg000003999c1+.nc
#  The output will be a file called
#    xqbwco#pr000003999c1+.nc  (I have used r to stand for rho).
#
#

# Import necessary libraries
import iris
import iris.quickplot as qplt
from iris.cube import CubeList
import numpy as np
import matplotlib.pyplot as plt
import gsw
from scipy.interpolate import interp1d
import sys

                        

#============================================================================
def calculate_density(filename):
    """
    reads in temperature and salinity and converts to density (at 2000m)
    """

    # read in data
    T_cube = iris.load_cube(filename,'insitu_T')
    #print(T_cube)
    S_cube = iris.load_cube(filename,'salinity')
    S_cube.data = (S_cube.data * 1000.) + 35.0 # convert to psu
    latitude = T_cube.coord('latitude').points
    longitude = T_cube.coord('longitude').points
    depth = T_cube.coord('depth_1').points
    latmesh,depmesh,lonmesh = np.meshgrid(latitude,depth,longitude)

    # convert to absolute salinity and conservative temperature
    pressmesh=gsw.p_from_z(depmesh * (-1.0), latmesh)
    SA = gsw.SA_from_SP(S_cube.data, pressmesh,lonmesh,latmesh)
    CT = gsw.CT_from_t(SA, T_cube.data, pressmesh)

    # calculate density (rho potential) referenced to 2000m
    rho_potential = gsw.rho(SA,CT,2000) - 1000.
    #print(rho_potential.mask)
    rho_potential = np.where(rho_potential.mask,-99999,rho_potential)
    
    rho_potential_cube = T_cube.copy(data=rho_potential)
    iris.util.mask_cube(rho_potential_cube,T_cube.data.mask,in_place=True)
    rho_potential_cube.long_name = 'density (calculated from T and S'
    rho_potential_cube.units=None
    rho_potential_cube.attributes=None

    return rho_potential_cube
#==================================================
def convert_to_vgrid(filename,T_grid_cube):
    """
    converts from T_grid_cube to V_grid_cube (only interpolate horizontally)
    """
    
    temporary = iris.load_cube(filename,'TOTAL OCEAN V-VELOCITY      CM S**-1')
    grid_cube = temporary[0,0,:,:]

    V_grid_cube = T_grid_cube.regrid(grid_cube,iris.analysis.Nearest())
    iris.util.mask_cube(V_grid_cube,temporary.data.mask,in_place=True)
   
    V_grid_cube.long_name = 'density on v-grid'

    #iris.save([T_grid_cube,V_grid_cube],'temporary.nc',fill_value = -99999.)

    return V_grid_cube

#==================================================
def get_sigma_levels():
    """
    get sigma levels. The user can change this
    """

    sigma_levels_a = np.arange(29.0,34.0,0.5)
    sigma_levels_b = np.arange(34.0,36.4,0.2)
    sigma_levels_c = np.arange(36.4,37.0,0.05)
    sigma_levels_d = np.arange(37.01,37.26,0.01)
    sigma_levels_e = np.arange(37.26,37.281,0.002)
    sigma_levels_f = np.arange(37.29,37.4,0.01)
    sigma_levels_g = np.arange(37.4,37.7,0.05)
    sigma_levels_h = np.arange(37.7,38.5,0.2)


    #sigma_levels_a = np.arange(29.0,36.0,0.25)
    #sigma_levels_b = np.arange(36.0,36.4,0.1)
    #sigma_levels_c = np.arange(36.4,37.2,0.02)
    #sigma_levels_d = np.arange(37.21,37.4,0.01)
    #sigma_levels_e = np.arange(37.4,37.7,0.02)
    #sigma_levels_f = np.arange(37.71,38.5,0.1)
    
    sigma_levels = np.concatenate((sigma_levels_a,sigma_levels_b,sigma_levels_c,
                                   sigma_levels_d,sigma_levels_e,sigma_levels_f,
                                   sigma_levels_g,sigma_levels_h))

    return sigma_levels

#================================================================
def interpolate_to_rhopot(filename,field,density_V_grid_cube,
                                             sigma_levels):
    """
    this will interpolate the field to the sigma levels
    """
    if field == 'V':
        field_cube = iris.load_cube(filename,
                                    'TOTAL OCEAN V-VELOCITY      CM S**-1')/100.
        field_data = field_cube.data

    if field == 'depth':
        field_cube = iris.load_cube(filename,
                                    'TOTAL OCEAN V-VELOCITY      CM S**-1')/100.
        depths=field_cube.coord('depth_1').points

        field_data = np.ma.copy(field_cube.data)
        
        for k in range(0,len(depths)):
            field_data[:,k,:,:]=depths[k]

        field_data.mask = field_cube.data.mask
            
          
    density_V_grid_data = density_V_grid_cube.data
    lons = field_cube.coord('longitude').points
    lats = field_cube.coord('latitude').points

    interpolated_variable=np.ma.zeros((len(sigma_levels),len(lats),len(lons)))
    interpolated_variable[:,:,:]=-99999.


    for j in range(0,len(lats)):
        for i in range(0,len(lons)):
            rho_profile = density_V_grid_data[0,:,j,i]
            var_profile = field_data[0,:,j,i]
            if not var_profile.mask[0]:
                rho_prof_red=[]
                var_prof_red=[]
                for k in range(0,20):
                    if not var_profile.mask[k]:
                        rho_prof_red.append(rho_profile[k])
                        var_prof_red.append(var_profile[k])
                #print('here',i,j)
                new_var = np.interp(sigma_levels,rho_prof_red,
                                    var_prof_red,
                                    left=-99999.,
                                    right=-99999.)
                #print(new_var)
                #print(var_profile)
                #print(rho_profile)
                #print(var_profile.mask[0])
                #sys.exit(0)
                interpolated_variable[:,j,i]=new_var

    #sys.exit(0)
    interpolated_variable.mask = np.where(interpolated_variable < -9999.,
                                          1.0,0.0)
    # set up a cube of the correct shape
    interpolated_cube_2d = field_cube[0,0:1,:,:].copy()
    interpolated_cube_2d.attributes=None
    interpolated_cube_2d.remove_coord('t')
    cube_list = CubeList([])
    for k in range(0,len(sigma_levels)):
        interpolated_cube_2d.coord('depth_1').points = sigma_levels[k]
        cube_list.append(interpolated_cube_2d.copy())
    
    iris.util.equalise_attributes(cube_list)
    interpolated_cube=cube_list.concatenate_cube()
    interpolated_cube.coord('depth_1').rename('sigma')
    interpolated_cube.coord('sigma').units='kg.m-3'

    # set up the cube we want with the correct data
    interpol_final_cube = interpolated_cube.copy(data=interpolated_variable)
    interpol_final_cube.long_name = field + ' on sigma coordinates'
    interpol_final_cube.units='m'
  
    
    return interpol_final_cube

#============================================================================
def get_thickness(dep_rho_potential_cube,filename):
    """
    gets the thickness of each layer for every point there is a velocity
    """

    # stuff we need
    lons=dep_rho_potential_cube.coord('longitude').points
    lats=dep_rho_potential_cube.coord('latitude').points
    sigma=dep_rho_potential_cube.coord('sigma').points
    dep_data = dep_rho_potential_cube.data
  
    # first calculate the base of each level from the standard um file - w grid
    w_cube = iris.load_cube(filename,'VERT.VEL. ON OCEAN HALF LEVELS  CM/S')
    depths = w_cube.coord('depth').points
    top_depths = np.zeros(21)
    for k in range(0,19):
        top_depths[k+1] = depths[k]
    top_depths[20]=depths[18]+depths[18]-depths[17]
  
    #now get the depth of the ocean
    v_cube = iris.load_cube(filename,'TOTAL OCEAN V-VELOCITY      CM S**-1')
    ocean_depth = np.zeros((len(lats),len(lons)))
    for j in range(0,len(lats)):
        for i in range(0,len(lons)):
            ocean_depth[j,i]=top_depths[20] # for all are found
            for k in range(19,-1,-1):
                if v_cube.data.mask[0,k,j,i]:
                    ocean_depth[j,i]=top_depths[k]

   
    # loop over all points
    #   if v is set
    #      thickness[k] = (depth[k]-depth[k-1]) / 2.0 + (depth[k+1]-depth[k])/2
    #   if k+1 is not set assume it is at the bottom of the ocean so use
    #                                                         base_depths
    #   if k-1 is not set assume it is top of ocean so use 0

    thickness=np.ma.zeros((len(sigma),len(lats),len(lons)))
    for j in range(0,len(lats)):
        for i in range(0,len(lons)):
            for k in range(1,len(sigma)-1):
                if dep_data[k,j,i] > -900:
                    
                    # get thickness of lower part of gridbox
                    if dep_data[k+1,j,i] < 0.0 or dep_data.mask[k+1,j,i]:
                        # to base of ocean
                        kp1_dz = ocean_depth[j,i] - dep_data[k,j,i]
                    else:
                        # halfway to next level
                        kp1_dz = (dep_data[k+1,j,i] - dep_data[k,j,i])/2.0
                        
                    # get thickness of upper part of gridbox
                    if dep_data[k-1,j,i] < 0.0 or dep_data.mask[k-1,j,i]:
                        # to top of ocean
                        km1_dz = dep_data[k,j,i]
                    else:
                        # half way to previous level
                        km1_dz = (dep_data[k,j,i] - dep_data[k-1,j,i])/2.0
                    
                    thickness[k,j,i] = kp1_dz + km1_dz
                    #print(k,j,i,thickness[k,j,i],kp1_dz,km1_dz)
                    #print(dep_data.mask[k-1,j,i],dep_data[k,j,i])
                    #sys.exit(0)

    #thickness.mask = dep_rho_potential_cube.data.mask
    thickness_cube = dep_rho_potential_cube.copy(data=thickness)
    thickness_cube.long_name = 'thickness (dz)'
    thickness_cube.units = 'm'
        

    return thickness_cube
            

def get_V(filename,thickness_cube):
    """
    this is an way of calculating V on density levels
    which is designed to maintain vdz it will do the following for each i and j
    
    1. find the base of each depth on the v grid (called depth_base)
    2. for each the thickess cube and find the depth of each layer
       on the sigma grid
    3. find out where this thickness grid corresponds to on the depth grid
    4. calculate a weighted average of v
    """

    # read in V cube from standard file
    V_orig_cube = iris.load_cube(filename,
                                    'TOTAL OCEAN V-VELOCITY      CM S**-1')/100.
    V_data = V_orig_cube.data
    lats = V_orig_cube.coord('latitude').points
    lons = V_orig_cube.coord('longitude').points

    V_dens = np.ma.zeros(np.shape(thickness_cube))
 

    
    # first calculate the base of each level from the standard um file - w grid
    w_cube = iris.load_cube(filename,'VERT.VEL. ON OCEAN HALF LEVELS  CM/S')
    depths = w_cube.coord('depth').points
    top_depths = np.zeros(21)
    for k in range(0,19):
        top_depths[k+1] = depths[k]
    top_depths[20]=depths[18]+depths[18]-depths[17]

    # loop over layers to see where sigma corresponds to depths
    sigma = thickness_cube.coord('sigma').points
    thickness = thickness_cube.data
    for j in range(0,len(lats)):
        for i in range(0,len(lons)):
            depmin=0.0 # how much has been accounted for
                                    # by layers above
            for k in range(0,len(sigma)):
                if thickness[k,j,i] > 0.0:
                    #print(i,j)
                    # find min and max depth of this layer
                    depmax=depmin+thickness[k,j,i]
                    depreq=np.zeros(20)
                    # find how many metres is in each layer on orig grid
                    for k2 in range(0,20):
                        if top_depths[k2] <= depmin <= top_depths[k2+1]:
                            # top band of sigma layer is in this depth band
                            # find number of meters it contributes
                            depreq[k2]=np.min([top_depths[k2+1],depmax])-depmin
                            #print('case A',depmin,depmax,top_depths[k2],
                                  #top_depths[k2+1])
                        elif top_depths[k2] <= depmax <= top_depths[k2+1]:
                            # bottom band of sigma layer is in this depth band
                            # find number of meters it contributes
                            depreq[k2]=depmax - top_depths[k2]
                            #print('case B')
                        elif (depmax > top_depths[k2+1] >
                              top_depths[k2] > depmin):
                            # our sigma level encompases this whole band
                            # so depth required from this band is full depth
                            depreq[k2]=top_depths[k2+1] - top_depths[k2]
                            #print('case C')
                    # calculate V_dens as a weighted average of the
                    # depths which contribute to this layer
                    vavg=0
                    for k2 in range(0,20):
                        if depreq[k2] > 0.1:
                            #print('j1',k,V_data[0,k2,j,i], depreq[k2])
                            vavg = vavg + (V_data[0,k2,j,i] * depreq[k2])
                    vavg=vavg / thickness[k,j,i]
                    if np.isfinite(vavg):
                        V_dens[k,j,i]=vavg
                        #print('j2',k,vavg*100.,thickness[k,j,i])
                        #print(' ')
                    else:
                        print('not finite',i,j,k,vavg,thickness[k,j,i])
                        for k in range(0,20):
                            print(k,depreq[k],V_data[0,k,j,i])
                        sys.exit(0)
                    
                   #print('vavg for ',k,vavg)
                    depmin = depmin + thickness[k,j,i]

                    
    V_dens = np.ma.where(thickness_cube.data < -9999., -99999., V_dens)
    V_dens.mask = thickness_cube.data.mask
    V_dens_cube = thickness_cube.copy(data=V_dens)
    V_dens_cube.units='m.s-1'


    return V_dens_cube
                         

def get_mask(maskfile):
    """
    get the masks for all the basin.  Calculated using make_basin_mask.py
    get them all on a V-grid
    """

    mask_cubes = CubeList([])
    mask_names = []
    cube_Atlantic = iris.load_cube(maskfile,'Atlantic mask V_grid')
    mask_cubes.append(cube_Atlantic)
    mask_names.append('Atlantic')
    
    cube_Pacific = iris.load_cube(maskfile,'Pacific mask V_grid')
    mask_cubes.append(cube_Pacific)
    mask_names.append('Pacific')
   
    cube_Indian = iris.load_cube(maskfile,'Indian mask V_grid')
    mask_cubes.append(cube_Indian)
    mask_names.append('Indian')
   
    cube_Global = iris.load_cube(maskfile,'Global mask V_grid')
    mask_cubes.append(cube_Global)
    mask_names.append('Global')


    return mask_cubes,mask_names


def get_params(V_rho_potential_cube,Atlmask_cube):
    """
    get grid spacings
    """
    
    lats = V_rho_potential_cube.coord('latitude').points
    coslats= np.cos(lats * 2. * np.pi / 360.)
    lons=V_rho_potential_cube.coord('longitude').points

    # length of a longitude box at the equator
    dx=(lons[1]-lons[0]) * 111320.

    return (dx,lats,coslats)

################################################################
def calc_stream(V_cube,dx,thickness_cube,coslats,mask_name, mask_cube):
    """
    calculate the stream function for this basin
    """
    
    # mask the data as appropriate
    
    sigma_coord = V_cube.coord('sigma').points
    vdz_data = np.ma.zeros(np.shape(V_cube.data))
    for k in range(0,len(sigma_coord)):
        vdz_data[k,:,:]=(V_cube.data[k,:,:] *
                         thickness_cube.data[k,:,:] * mask_cube.data)

    #calculate zonal integral of v*dz
    # (note everything we dont want to use should be masked)
    #and multiply it by the width of the gridbox (ie dx * cos latitude)

    nz=len(V_cube.coord('sigma').points)
    ny=len(V_cube.coord('latitude').points)
    vdz_cube=V_cube.copy(data=vdz_data)
    vdz_cube.long_name = 'vdz_cube'    
    vdz_cube.data.mask = np.where(vdz_cube.data < -999,1.0,0.0)


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
    phi_data.mask = np.where(ztotalV_cube.data.mask == 1.0, 1.0, 0.0)
    phi_cube=ztotalV_cube.copy(data=phi_data / 1.0E6)
    phi_cube.long_name = mask_name + ' Meridional Overturning Circulation'
    phi_cube.units = 'Sv'
    
    # checks - no longer needed
    # read in the original for plotting
    #orig_AMOC = iris.load_cube('xqbwco#pk000003999c1+.nc',
    #                'Meridional Overturning Stream Function (Global)')
    #orig_AMOC = iris.util.squeeze(orig_AMOC)

    # try calculating in reverse
    phi_data_rev=np.zeros(np.shape(ztotalV_cube.data))
    for k in range(nz-1,-1,-1):
        if k == nz-1:
            phi_data_rev[k,:]= (-1.0) *  ztotalV_cube.data[k,:]
        else:
            phi_data_rev[k,:]=(phi_data_rev[k+1,:] - ztotalV_cube.data[k,:])
    phi_rev_cube=ztotalV_cube.copy(data=phi_data_rev / 1.0E6)
    phi_rev_cube.long_name = mask_name + ' Meridional Overturning Circulation (reversed'
    phi_rev_cube.units = 'Sv'
  

    #for j,lat in enumerate(lats):
    #    print(lat,j,phi_data_rev[0,j]/1.0E6)
        
    #sys.exit(0)
    #print('total2',np.sum(ztotalV_cube.data[:,91] / 1.0E6))

    return phi_cube, phi_rev_cube
                       

#=============================================================================

filename = 'xqbwco#pg000003999c1+.nc'
fileout = filename[0:8]+'r' + filename[9::]

# 1. calculate density on Temp, Salinity grid
density_T_grid_cube = calculate_density(filename)

# 2. put density cube on V grid using interpolation (horizontally)
density_V_grid_cube = convert_to_vgrid(filename,density_T_grid_cube)

# 3. get sigmal levels on which we will interpolate the fields

sigma_levels = get_sigma_levels()

# 3. interpolate depth to rho_potential

dep_rho_potential_cube = interpolate_to_rhopot(filename,'depth',
                                             density_V_grid_cube,
                                             sigma_levels)
print('got depths')

#4. We need to find a dz value (call this thickness)
thickness_cube = get_thickness(dep_rho_potential_cube,filename)
print('got thickness')
output_cubes = CubeList([dep_rho_potential_cube,
                         thickness_cube])      

#5. Convert velocity on depth levels to velocity on rho levels.  We do this
#   by finding which depths each rho level encompasses.  We then get the
#   V values from those depth levels

V_rho_potential_cube = get_V(filename,thickness_cube)
V_rho_potential_cube.long_name = 'V on density levels'
output_cubes.append(V_rho_potential_cube)
print('got V')
iris.save(output_cubes,fileout,fill_value=-99999.)
sys.exit(0)


#6.  We now have all the fields we need so find the AMOC / PMOC / IMOC / GMOC
#6a get the parameters we need
mask_cubes, mask_names = get_mask('masks.nc')
(dx,lats,coslats) = get_params(V_rho_potential_cube,mask_cubes[0])

#6b. find the overturning circulation for each of the masks
for i in range(0,len(mask_names)):
    MOC_cube, MOC_rev_cube = calc_stream(V_rho_potential_cube,dx,
                           thickness_cube,coslats,mask_names[i],
                           mask_cubes[i])
    output_cubes.append(mask_cubes[i])
    output_cubes.append(MOC_cube)
    output_cubes.append(MOC_rev_cube)



iris.save(output_cubes,fileout,fill_value=-99999.)


  
