#NAME
#    regrid density coordinates
#PURPOSE 
#
#  This program will regrid a pg file onto density coordinates
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
def get_sigma_levels(exptname):
    """
    get sigma levels. The user can change this
    """
    sigma_levels=[]
    try:
        # Open the file in read mode
        with open(exptname + '_sigmas.txt', 'r') as file:
            for line_number, line in enumerate(file, start=1):
                 # Strip newline characters and whitespace
                 line = line.strip()
            
                 # Split the line by colon
                 split_line = line.split(':')
                 value = split_line[1].strip()
                 file_sigma = float(value)
                 sigma_levels.append(file_sigma)
    except FileNotFoundError:
        print(f"The file '{exptname}'_sigmas.txt was not found.")
        sys.exit(0)
    except IOError:
        print(f"An error occurred while reading the file '{filename}'.")
        sys.exit(0)

    sigma_array = np.array(sigma_levels)
   
    return sigma_array

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

    interpolated_variable=np.zeros((len(sigma_levels),len(lats),len(lons)))
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
    interpolated_cube_1d = field_cube[0,0:1,:,:].copy()
    
    interpolated_cube_1d.attributes=None
    interpolated_cube_1d.remove_coord('t')
    cube_list = CubeList([])
    for k in range(0,len(sigma_levels)):
        interpolated_cube_1d.coord('depth_1').points = sigma_levels[k]
        cube_list.append(interpolated_cube_1d.copy())

    
    iris.util.equalise_attributes(cube_list)
    interpolated_cube=cube_list.concatenate_cube()
    interpolated_cube.coord('depth_1').rename('sigma')

    interpol_final_cube = interpolated_cube.copy(data=interpolated_variable)
    interpol_final_cube.long_name = field + ' on sigma coordinates'

    
    return interpol_final_cube

#============================================================================
def get_thickness(dep_rho_potential_cube,v_rho_potential_cube,
                               filename):
    """
    gets the thickness of each layer for every point there is a velocity
    """

    # stuff we need
    lons=v_rho_potential_cube.coord('longitude').points
    lats=v_rho_potential_cube.coord('latitude').points
    sigma=v_rho_potential_cube.coord('sigma').points
    v_data = v_rho_potential_cube.data
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

    thickness=np.zeros((len(sigma),len(lats),len(lons)))
    for j in range(0,len(lats)):
        for i in range(0,len(lons)):
            for k in range(1,len(sigma)-1):
                if v_data[k,j,i] > -900:
                    
                    # get thickness of lower part of gridbox
                    if dep_data[k+1,j,i] < 0.0:
                        # to base of ocean
                        kp1_dz = ocean_depth[j,i] - dep_data[k,j,i]
                    else:
                        # halfway to next level
                        kp1_dz = (dep_data[k+1,j,i] - dep_data[k,j,i])/2.0
                        
                    # get thickness of upper part of gridbox
                    if dep_data[k-1,j,i] < 0.0:
                        # to top of ocean
                        km1_dz = dep_data[k,j,i]
                    else:
                        # half way to previous level
                        km1_dz = (dep_data[k,j,i] - dep_data[k-1,j,i])/2.0
                    
                    thickness[k,j,i] = kp1_dz + km1_dz
                    #print(k,j,i,thickness[k,j,i],kp1_dz,km1_dz)
                    #print(dep_data[k-1,j,i],dep_data[k,j,i])

                    #sys.exit(0)

                    #print('found',sigma[k],lats[j],lons[i],
                    #          dep_data[k+1,j,i],
                    #          dep_data[k,j,i],dep_data[k-1,j,i],dz,
                    #          v_data[k,j,i],v_data[k-1,j,i],v_data[k+1,j,i])

    thickness_cube = dep_rho_potential_cube.copy(data=thickness)
    thickness_cube.long_name = 'thickness (dz)'
        

    return thickness_cube
            

def get_alt_V(filename,thickness_cube,V_orig_density_cube):
    """
    this is an alternative way of calculating V on density levels
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

    V_alt = np.zeros(np.shape(V_orig_density_cube))
 

    
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
    #for j in range(0,len(lats)):
    #    for i in range(0,len(lons)):
    for j in range(14,15):
        for i in range(138,139):
            depmin=0.0 # how much has been accounted for
                                    # by layers above
            for k in range(0,len(sigma)):
                if thickness[k,j,i] > 0.0:
                    # find min and max depth of this layer
                    depmax=depmin+thickness[k,j,i]
                    print(k,thickness[k,j,i],depmin,depmax)
                    depreq=np.zeros(20)
                    # find how many metres is in each layer on orig grid
                    for k2 in range(0,20):
                        if top_depths[k2] <= depmin <= top_depths[k2+1]:
                            # top band of sigma layer is in this depth band
                            # find number of meters it contributes
                            depreq[k2]=np.min([top_depths[k2+1],depmax])-depmin
                            print('case A',k2,depmin,depmax,top_depths[k2],
                                  top_depths[k2+1])
                        elif top_depths[k2] <= depmax <= top_depths[k2+1]:
                            # bottom band of sigma layer is in this depth band
                            # find number of meters it contributes
                            depreq[k2]=depmax - top_depths[k2]
                            print('case B',k2,top_depths[k2],depmax,top_depths[k2+1])
                        elif (depmax > top_depths[k2+1] >
                              top_depths[k2] > depmin):
                            # our sigma level encompases this whole band
                            # so depth required from this band is full depth
                            depreq[k2]=top_depths[k2+1] - top_depths[k2]
                            print('case C')
                    # calculate V_alt as a weighted average of the
                    # depths which contribute to this layer
                    vavg=0
                    for k2 in range(0,20):
                        if depreq[k2] > 0.1:
                            #print('j1',k,V_data[0,k2,j,i], depreq[k2])
                            vavg = vavg + (V_data[0,k2,j,i] * depreq[k2])
                    vavg=vavg / thickness[k,j,i]
                    if np.isfinite(vavg):
                        V_alt[k,j,i]=vavg
                        #print('j2',k,vavg*100.,thickness[k,j,i])
                        #print(' ')
                    else:
                        print('not finite',i,j,lons[i],lats[j],k,
                              vavg,thickness[k,j,i])
                        for k in range(0,20):
                            print(k,top_depths[k+1],top_depths[k],V_data[0,k,j,i])
                        sys.exit(0)
                    
                   #print('vavg for ',k,vavg)
                    depmin = depmin + thickness[k,j,i]

    V_alt = np.where(V_orig_density_cube.data < -9999., -99999., V_alt)
    V_alt_cube = V_orig_density_cube.copy(data=V_alt)


    return V_alt_cube
                         


                    
 
#=============================================================================

exptname = 'xqbwd'
startyear=3971
endyear=3972
filestart = '/home/earjcti/um/' + exptname + '/'

for year in range(startyear,endyear):
    print(year)
    filename = filestart + 'pg/' + exptname + 'o#pg00000' + str(year) + 'c1+.nc'

    # 1. calculate density on Temp, Salinity grid
    density_T_grid_cube = calculate_density(filename)

    # 2. put density cube on V grid using interpolation (horizontally)
    density_V_grid_cube = convert_to_vgrid(filename,density_T_grid_cube)

    # 3. get sigmal levels on which we will interpolate the fields
    
    sigma_levels = get_sigma_levels(exptname)

    # 3. interpolate fields to rho_potential (depth and v)

    dep_rho_potential_cube = interpolate_to_rhopot(filename,'depth',
                                                   density_V_grid_cube,
                                                   sigma_levels)

    v_rho_potential_cube = interpolate_to_rhopot(filename,'V',
                                                 density_V_grid_cube,
                                                 sigma_levels)

    #4. If we have a velocity we want a thickness
    thickness_cube = get_thickness(dep_rho_potential_cube,v_rho_potential_cube,
                                   filename)

    #5. Alternative way of interpolating V using thickness

    V_alt_cube = get_alt_V(filename,thickness_cube,v_rho_potential_cube)
    V_alt_cube.long_name = 'alternative V cube'

    fileout = filestart + 'Vdz/' + exptname + '_vdz_' + str(year) + '.nc'    
    iris.save([thickness_cube,dep_rho_potential_cube,v_rho_potential_cube,V_alt_cube],fileout,fill_value=-99999.)


  
