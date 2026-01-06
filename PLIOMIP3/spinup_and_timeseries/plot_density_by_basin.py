#NAME
#    plot density by basin
#PURPOSE 
#
#  This program will do two things.
#  1. It will do a latitude-depth plot of density temperature and salinity
#     for each basin
#  2. It will try and optimise the number of density classes so that
#     we have the same number of gridpoints in each density class

# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import gsw
import sys
import pandas as pd
from pathlib import Path
#from netCDF4 import Dataset, MFDataset
#from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

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
    S_cube.attributes.pop("valid_min")
    S_cube.attributes.pop("valid_max")
    
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

    return rho_potential_cube,T_cube,S_cube



def get_mask(filename,basin):
    """
    gets the mask from the Atlantic mask file on a v grid
    """

    cube = iris.load_cube(filename,basin + ' mask T_grid')

    return cube



####################################################################
def mask_function(cube,mask_cube):
    """
    masks the cubes
    """

    # mask V cube
    cube=iris.util.squeeze(cube)
    
    cubelev=np.copy(cube.data)
    depth_coord = cube.coord('depth_1').points
    for k in range(0,len(depth_coord)):
        cubelev[k,:,:] = np.ma.where(mask_cube.data > 0.5, cube.data[k,:,:],
                              -99999.)

    returncube = cube.copy(data=cubelev)
    returncube.data = np.ma.where(returncube.data > 1.0E20,-99999.,
                                  returncube.data)
    returncube.data.mask = np.where(returncube.data < -9999,1.0,0.0)
    #cube.data.mask = np.where(cube.data ==-99999., 1.0,0.0)
    

    return returncube


def basin_avg(cube):
    """
    calculates the basin average across the cube
    """

    cubedata = cube.data
    #print(cubedata.mask)
    zonalmean_data = np.ma.mean(cube.data,axis=2)
    #print(np.shape(cube.data))
    #print(np.shape(zonalmean_data))
    basin_avg_cube = cube.collapsed('longitude',iris.analysis.MEAN)
    basin_avg_cube = basin_avg_cube.copy(data=zonalmean_data)
    #print(zonalmean_data)
    #print(basin_avg_cube)
    #print(basin_avg_cube.data[0,90])
    #print(zonalmean_data[0,90])
    #print(cube.data[0,90,:])
    #print(cube.data.mask[0,90,:])
    #sys.exit(0)

    return basin_avg_cube



def quantile_binning(full_data, num_bins,filename):
    """
    Copilot wrote this subsection
    Bins data into quantile-based bins with roughly equal number of points 
    per bin.
    
    Parameters:
        data (array-like): The input data to bin.
        num_bins (int): Number of bins to create.
    
    Returns:
        bin_edges (np.ndarray): The edges of the bins.
        bin_labels (pd.Series): Labels indicating which bin each data point 
        belongs to.
    """
    print(len(full_data.flatten()))
    data = []
    for point in full_data.flatten():
        if point > 0:
            data.append(point)
    print(len(data))
    data = pd.Series(data)
    
    # Create quantile-based bins
    bin_labels = pd.qcut(data, q=num_bins, labels=False, duplicates='drop')
    
    # Extract bin edges
    bin_edges = pd.qcut(data, q=num_bins, duplicates='drop').unique().categories

    #  Calculate bin centres
    bin_centres = np.array([(interval.left + interval.right) / 2 for interval in bin_edges])

    # write bin centres to a text file
    f=open(filename, 'w') 
    for i, centre in enumerate(bin_centres):
        f.write(f"Bin {i} centre: {centre:.6f}\n")
    f.close()
    
    return bin_edges, bin_labels, bin_centres


#####################################################################
def process_data(filename,basin,year):
    """
    processess all the data for each year
    """


    # file where V and dz are stored on density coordinates

    # 1. calculate density on Temp, Salinity grid
    (density_T_grid_cube,
     temperature_cube,
     salinity_cube) = calculate_density(filename)

    # gets the basins over which we calculate
    mask_cube=get_mask('masks.nc',basin) # get mask on T grid

    density_T_grid_cube = mask_function(density_T_grid_cube,mask_cube)
    temperature_cube = mask_function(temperature_cube,mask_cube)
    salinity_cube = mask_function(salinity_cube,mask_cube)

    # calculate stream function for the basin
    density_basin_avg_cube = basin_avg(density_T_grid_cube)
    temperature_basin_avg_cube = basin_avg(temperature_cube)
    salinity_basin_avg_cube = basin_avg(salinity_cube)


    # put data into unequally spaced bins if global and we haven't already done this

    filename = exptname + '_sigmas.txt'
    filepath = Path(filename)

    if basin == 'Global' and not filepath.exists():
        edges, labels, bin_centres = quantile_binning(density_T_grid_cube.data,
                                                      90,filename)

        print("Bin edges:")
        for i, interval in enumerate(edges):
            print(f"Bin {i}: {interval}")

    # save everything to a file

    fileout = ('/uolstore/Research/a/hera1/earjcti/um/' + exptname +
               '/basin_diagnostics/' + exptname + '_' +
               basin + str(year) + '.nc')

    density_basin_avg_cube.long_name='density basin'
    temperature_basin_avg_cube.long_name='temperature basin'
    salinity_basin_avg_cube.long_name='salinity basin'
    iris.save([density_T_grid_cube,
               density_basin_avg_cube,temperature_basin_avg_cube,
               salinity_basin_avg_cube],
               fileout,fill_value = -99999.)


###########################################################    
def plot_anomaly(exptname,cntlname,startyear,endyear,cntlstart,cntlend,basin):
    """
    plots the difference in density salinity and temperature between two 
    experiments
    """

    exptfile=('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/mean_' + exptname + '_' +
                  basin + str(startyear) + '_' + str(endyear-1)+'.nc')
    cntlfile=('/home/earjcti/um/' + cntlname +
                  '/basin_diagnostics/mean_' + cntlname + '_' +
                  basin + str(cntlstart) + '_' + str(cntlend-1)+'.nc')

    expt_temp_cube = iris.load_cube(exptfile,'temperature basin')
    cntl_temp_cube = iris.load_cube(cntlfile,'temperature basin')

    expt_dens_cube = iris.load_cube(exptfile,'density basin')
    cntl_dens_cube = iris.load_cube(cntlfile,'density basin')

    expt_sal_cube = iris.load_cube(exptfile,'salinity basin')
    cntl_sal_cube = iris.load_cube(cntlfile,'salinity basin')

    temp_anom = expt_temp_cube - cntl_temp_cube
    dens_anom = expt_dens_cube - cntl_dens_cube
    sal_anom = expt_sal_cube - cntl_sal_cube

    vals = np.arange(-1.0,7.5,0.5)
    if cntlname !='xqbwc':
        vals = np.arange(-2.0,2.5,0.5)
    qplt.contourf(temp_anom,levels=vals,extend='both')
    plt.title(basin + ' Temperature anomaly: '+ exptname + '-'+ cntlname)
    plt.savefig('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/meanT_' + exptname + '-' + cntlname + '_'
                + basin + str(startyear) + '_' + str(endyear-1)+'.png')
    plt.close()

    vals = np.arange(-0.5,0.55,0.05)
    qplt.contourf(dens_anom,levels=vals,extend='both',cmap='RdBu_r')
    plt.title(basin + ' Density anomaly: '+ exptname + '-'+ cntlname)
    plt.savefig('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/meandens_' + exptname + '-' +
                cntlname + '_'
                + basin + str(startyear) + '_' + str(endyear-1)+'.png')
    plt.close()

    qplt.contourf(sal_anom,levels=vals,extend='both',cmap='RdBu_r')
    plt.title(basin + ' Salinity anomaly: '+ exptname + '-'+ cntlname)
    plt.savefig('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/meansal_' + exptname + '-' +
                cntlname + '_' + basin + str(startyear) + '_' +
                str(endyear-1)+'.png')
    plt.close()
   
    
    sys.exit(0)


def plot_Pacific_and_Atlantic(exptname,startyear,endyear,lev):
    """
    plots the difference in SST and SSS between the Pacific
    and the Atlantic for the experiment
    """
    Atlfile=('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/mean_' + exptname + '_Atlantic' +
                  str(startyear) + '_' + str(endyear-1)+'.nc')

    Atl_temp_cube = iris.load_cube(Atlfile,'temperature basin')

    Atl_dens_cube = iris.load_cube(Atlfile,'density basin')

    Atl_sal_cube = iris.load_cube(Atlfile,'salinity basin')


    Pacfile=('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/mean_' + exptname + '_Pacific' +
                  str(startyear) + '_' + str(endyear-1)+'.nc')

    Pac_temp_cube = iris.load_cube(Pacfile,'temperature basin')

    Pac_dens_cube = iris.load_cube(Pacfile,'density basin')

    Pac_sal_cube = iris.load_cube(Pacfile,'salinity basin')

    Atlcntl=('/home/earjcti/um/' + cntlname +
                  '/basin_diagnostics/mean_' + cntlname + '_Atlantic' +
                  str(startyear) + '_' + str(endyear-1)+'.nc')
    Paccntl=('/home/earjcti/um/' + cntlname +
                  '/basin_diagnostics/mean_' + cntlname + '_Pacific' +
                  str(startyear) + '_' + str(endyear-1)+'.nc')
    Atlcntl_cube = iris.load_cube(Atlcntl,'temperature basin')
    Paccntl_cube = iris.load_cube(Paccntl,'temperature basin')
    Atl_diff_cube = Atl_temp_cube - Atlcntl_cube
    Pac_diff_cube = Pac_temp_cube - Paccntl_cube

    # plot temperature
    plt.subplot(221)
    plt.plot(Atl_temp_cube.coord('latitude').points,Atl_temp_cube[lev,:].data,
             label='Atlantic')
    plt.plot(Pac_temp_cube.coord('latitude').points,Pac_temp_cube[lev,:].data,
             label='Pacific')
    plt.legend()
    plt.title('Basin temperature for : ' + period.get(exptname,exptname) + ' lev = ' + str(lev))
    plt.xlabel('latitude')                                     
    plt.ylabel('degC')
    plt.axvline(x=10)
    plt.grid(True)
  
    # plot salinity
    plt.subplot(222)
    plt.plot(Atl_sal_cube.coord('latitude').points,Atl_sal_cube[lev,:].data,
             label='Atlantic')
    plt.plot(Pac_sal_cube.coord('latitude').points,Pac_sal_cube[lev,:].data,
             label='Pacific')
    plt.legend()
    plt.title('Basin salinity for : ' + period.get(exptname,exptname))
    plt.xlabel('latitude')                                     
    plt.ylabel('psu')
    plt.axvline(x=10)
    plt.grid(True)


    # plot density
    plt.subplot(223)
    plt.plot(Atl_dens_cube.coord('latitude').points,Atl_dens_cube[lev,:].data,
             label='Atlantic')
    plt.plot(Pac_dens_cube.coord('latitude').points,Pac_dens_cube[lev,:].data,
             label='Pacific')
    plt.legend()
    plt.title('Basin dens for : ' + period.get(exptname,exptname))
    plt.xlabel('latitude')                                     
    plt.ylabel('density')
    plt.axvline(x=10)
    plt.grid(True)
    plt.tight_layout()   

    # plot anomaly
    plt.subplot(224)
    plt.plot(Atl_diff_cube.coord('latitude').points,Atl_diff_cube[lev,:].data,
             label='Atlantic')
    plt.plot(Pac_diff_cube.coord('latitude').points,Pac_diff_cube[lev,:].data,
             label='Pacific')
    plt.legend()
    plt.title('Basin anomaly for : ' + period.get(exptname,exptname) + '-' +
              period.get(cntlname,cntlname))
    plt.xlabel('latitude')                                     
    plt.ylabel('degc diff')
    plt.axvline(x=10)
    plt.grid(True)
    plt.tight_layout()   

    plt.show()
    sys.exit(0)
    
    plt.title(basin + ' Temperature anomaly: '+ exptname + '-'+ cntlname)
    plt.savefig('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/meanT_' + exptname + '-' + cntlname + '_'
                + basin + str(startyear) + '_' + str(endyear-1)+'.png')
    plt.close()

    vals = np.arange(-1.0,1.1,0.1)
    qplt.contourf(dens_anom,levels=vals,extend='both',cmap='RdBu_r')
    plt.title(basin + ' Density anomaly: '+ exptname + '-'+ cntlname)
    plt.savefig('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/meandens_' + exptname + '-' +
                cntlname + '_'
                + basin + str(startyear) + '_' + str(endyear-1)+'.png')
    plt.close()

    vals = np.arange(-1.0,1.1,0.1)
    qplt.contourf(sal_anom,levels=vals,extend='both',cmap='RdBu_r')
    plt.title(basin + ' Salinity anomaly: '+ exptname + '-'+ cntlname)
    plt.savefig('/home/earjcti/um/' + exptname +
                  '/basin_diagnostics/meansal_' + exptname + '-' +
                cntlname + '_' + basin + str(startyear) + '_' +
                str(endyear-1)+'.png')
    plt.close()
   
    
    sys.exit(0)


#######################################################################
exptname = 'xpsid'
startyear=12
endyear=2000
basin='Pacific'

period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400','xpsig':'EP',
          'xpsic':'PI','xqbwd':'LP','xqbwj':'LP490','xqbwe':'EP400',
          'xqbwg':'EP',
          'xqbwc':'PI'}


# get individual years diagnostics for the basin

for year in range(startyear,endyear):
    print(year)
    filestart = '/uolstore/Research/a/hera1/earjcti/um/' + exptname + '/pg/' + exptname 
    filename = filestart + 'o#pg'+str(year).zfill(9)+'c1+.nc'
    process_data(filename,basin,year)

#################################################
# get the mean dianostics for the basin
#dens_cubelist = CubeList([])
#sal_cubelist = CubeList([])
#temp_cubelist = CubeList([])

#for year in range(startyear,endyear):
#     filename = ('/home/earjcti/um/' + exptname +
#                  '/basin_diagnostics/' + exptname + '_' +
#                  basin + str(year) + '.nc')

#     dens_cubelist.append(iris.load_cube(filename,'density basin'))
#     temp_cubelist.append(iris.load_cube(filename,'temperature basin'))
#     sal_cubelist.append(iris.load_cube(filename,'salinity basin'))
#     sal_cube = iris.load_cube(filename,'salinity basin')
   
#iris.util.equalise_attributes(dens_cubelist)
#iris.util.equalise_attributes(sal_cubelist)
#iris.util.equalise_attributes(temp_cubelist)
#dens_cubes = dens_cubelist.merge_cube()
#sal_cubes = sal_cubelist.merge_cube()
#temp_cubes = temp_cubelist.merge_cube()

#dens_avg_cube = dens_cubes.collapsed('t',iris.analysis.MEAN)
#sal_avg_cube = sal_cubes.collapsed('t',iris.analysis.MEAN)
#temp_avg_cube = temp_cubes.collapsed('t',iris.analysis.MEAN)

#fileout = ('/home/earjcti/um/' + exptname +
#                  '/basin_diagnostics/mean_' + exptname + '_' +
#                  basin + str(startyear) + '_' + str(endyear-1)+'.nc')
#iris.save([dens_avg_cube,temp_avg_cube,sal_avg_cube],
#          fileout,fill_value = -99999.)
  

################################################
# plot anomalies
#cntlname = 'xpsig'
#cntlstart=startyear
#cntlend=endyear
#cntlstart=12
#cntlend=2000
#plot_anomaly(exptname,cntlname,startyear,endyear,cntlstart,cntlend,basin)
#plot_Pacific_and_Atlantic(exptname,startyear,endyear,0)  # the last number is the level

sys.exit(0)

