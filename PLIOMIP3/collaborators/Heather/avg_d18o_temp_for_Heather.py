#!/usr/bin/env python2.7
#NAME
#   avg_d18o_and_temp_for_heather
#PURPOSE
#   Heather wants some d18o_sw and temperature data to compare with data
#   I am also going to manually adjust so that the d18o_sw averages 0permil
#   for present; 0.3permil for Pliocene
#  
# Julia September 2024



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import cartopy.crs as ccrs
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess



def get_avg(year,weights):
    """
    gets the average d18o and temperature for this year
    year = year
    weights = weights for averaging.  Will be calculated here if they are unset
    """
 
   
    yearuse = str(year).zfill(9)
    filename=('/nfs/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    cubelist = iris.load(filename)
    #for cube in cubelist:
    #    print(cube.var_name, cube.long_name)
    #sys.exit(0)
    temperaturecube = cubelist.extract('insitu_T')[0]
    salinitycube = cubelist.extract('salinity')[0]
    ratio18ocube = cubelist.extract('otracer1')[0]

    ratio18ocube.coord('latitude').bounds = None
    ratio18ocube.coord('longitude').bounds = None
    try:
        del ratio18ocube.attributes["valid_max"]
    except:
        print("no valid max attribute")
    try:
        del ratio18ocube.attributes["valid_min"]
    except:
        print("no valid min attribute")
    
    
    if np.ndim(weights) == 1: # this has not yet been setup 
        # get weights for calculating mean - scale by depth
        ratio18ocube.coord('latitude').guess_bounds()
        ratio18ocube.coord('longitude').guess_bounds()
        weights = iris.analysis.cartography.area_weights(ratio18ocube)
        new_weights = np.zeros(np.shape(weights))
        depths = ratio18ocube.coord('depth_1').points
        layer_depths = np.zeros(20)
        layer_depths[0] = depths[0] * 2.0
        layer_depths[19] = depths[19] - depths[18]
        for k in range(1,19):
            layer_depths[k]=((depths[k+1] - depths[k-1]) * 0.5) 
        for k in range(0,20):
            new_weights[:,k,:,:] = weights[:,k,:,:] * layer_depths[k]
        ratio18ocube.coord('latitude').bounds=None
        ratio18ocube.coord('longitude').bounds=None
    else:
        new_weights = weights

    # calculate mean and add to the 
    mean18oratio_cube = ratio18ocube.collapsed(['longitude','latitude','depth_1'],
                                       iris.analysis.MEAN,
                                       weights=new_weights)

    # if pliocene adjust so that the meanratio is -0.3permil 
    # if preindustrial adjust so that the mean ratio is 0permil

    if gm_d18o.get(exptname) == '0permil':
        ratio18ocube_adj = ratio18ocube + 2005.2E-6 - mean18oratio_cube.data
    elif gm_d18o.get(exptname) == '-0.3permil':
        ratio18ocube_adj = ratio18ocube + 2004.6E-6 - mean18oratio_cube.data
    else:
        print('need a ratio for correcting to')
        sys.exit(0)

    print(ratio18ocube.coord('latitude'))

    return (new_weights, ratio18ocube, ratio18ocube_adj, temperaturecube,
            salinitycube)

#######################################################
def get_d18o_and_temperature():
    """
    gets the ocean salratioerature drift for exeriment = exptname between
    years startyear
    """

    # arrays for storing mean salratioeratures
    cubelist_18oratio = CubeList([])
    cubelist_18oadj = CubeList([])
    cubelist_temperature = CubeList([])
    cubelist_salcube = CubeList([])

    # obtain means for each year and store in the arrays
    weights = [0.0]
    for year in range(startyear,startyear+nyears):
        print(year)
        (weights, ratio18o_cube, adj18o_cube, temp_cube,
         sal_cube) = get_avg(year,weights)
        cubelist_18oratio.append(ratio18o_cube)
        cubelist_18oadj.append(adj18o_cube)
        cubelist_temperature.append(temp_cube)
        cubelist_salcube.append(sal_cube)

    
    # put into a single cube and average
    #print('j1',cubelist_18oratio)
    iris.util.equalise_attributes(cubelist_18oratio)
    iris.util.equalise_attributes(cubelist_temperature)
    iris.util.equalise_attributes(cubelist_18oadj)
    iris.util.equalise_attributes(cubelist_salcube)
    #print('j2',cubelist_18oratio)
    cubes18o = cubelist_18oratio.concatenate_cube()
    cubestemp = cubelist_temperature.concatenate_cube()
    cubessal = cubelist_salcube.concatenate_cube()
    cubes_18oadj = cubelist_18oadj.concatenate_cube()

    cubeavg18o=cubes18o.collapsed('t',iris.analysis.MEAN)
    cubeavg18o_adj=cubes_18oadj.collapsed('t',iris.analysis.MEAN)
    cubeavgtemp = cubestemp.collapsed('t',iris.analysis.MEAN)
    cubeavgsal = cubessal.collapsed('t',iris.analysis.MEAN)

    # convert 18o ratio to permil
    
    cubeavgd18o = (cubeavg18o - 2005.2E-6) / 2005.2E-9
    cubeavgd18o_adj = (cubeavg18o_adj - 2005.2E-6) / 2005.2E-9
    cubeavgd18o.long_name = 'd18o permile (raw model output)'
    cubeavgd18o_adj.long_name = ('d18o adjusted:  global mean =' 
                                + gm_d18o.get(exptname))


    # mask missing values
    print('j1',cubeavgd18o_adj.data)
   
    cubeavgd18o_adj.data = np.ma.where(cubeavgd18o_adj.data < -900.0, -99999.,
                                    cubeavgd18o_adj.data)
   
    #iris.util.mask_cube(cubeavgd18o,np.where(cubeavg18o.data > 1.0E10))
    iris.util.mask_cube(cubeavgd18o_adj,np.where(cubeavg18o_adj.data < -900.0))
    #iris.util.mask_cube(cubeavgtemp,np.where(cubeavg18o.data > 1.0E10))
    
    print('j2',cubeavgd18o_adj.data)
   
    filename = (timeperiod.get(exptname) + '_' + exptname + '_d18o_and_temp'
                + str(startyear) + '_' + str(startyear+nyears) + '.nc')


    
    iris.save([cubeavgd18o, cubeavgd18o_adj, cubeavgtemp,cubeavgsal],filename, 
              fill_value = -99999.)


def diff_two_files():
    """
    just compares the files produced here fore the pliocene and the pi
    """

    file1 = 'PI_xpsia_d18o_and_temp2970_3000.nc'
    file2 = 'MP_3.205Ma_xpsib_d18o_and_temp2970_3000.nc'
    cubes1 = iris.load(file1)
    cubes2 = iris.load(file2)

    for cube in cubes1:
        print(cube.var_name)
    temp1 =  cubes1.extract('insitu_T')[0]
    d18o1 = cubes1.extract('d18o_adjusted__global_mean__0permil')[0]
    temp2 =  cubes2.extract('insitu_T')[0]
    d18o2 = cubes2.extract('d18o_adjusted__global_mean___0_3permil')[0]

    print(temp2)
    print(temp1)
    tdiff=temp2-temp1
    tdiff.data = np.ma.where(tdiff.data < -900,-99999., tdiff.data)
    print(tdiff.data)
    d18odiff=d18o2-d18o1
    print(d18odiff)
    d18odiff.data = np.ma.where(d18odiff.data < -900,-99999., d18odiff.data)
   
    iris.save([tdiff,d18odiff],'plio-pi diff.nc',fill_value=-99999.)


def plot_ss_d18o():
    """
    just want to compare sea surface d18o to observations for pi
    as edmond didn't like the climate
    """

    #file1 = 'PI_xpsia_d18o_and_temp2970_3000.nc'
    file1 = 'MP_3.205Ma_xpsib_d18o_and_temp2970_3000.nc'

    cubes1 = iris.load(file1)
   
    for cube in cubes1:
        print(cube.var_name)
    #d18o_cube = cubes1.extract('d18o_adjusted__global_mean__0permil')[0]
    d18o_cube = cubes1.extract('d18o_adjusted__global_mean___0_3permil')[0]

    print(d18o_cube)
    d18o1_surf_cube = d18o_cube[0,:,:]
    print(d18o1_surf_cube)
    qplt.contourf(d18o1_surf_cube,cmap='RdBu_r',
                  levels=np.arange(-2.5,2.6,0.1),extend='both')
    plt.show()

################################
# main program

# annual mean
figureno=0

timeperiod = {'xpsia':'PI','xpsib':'MP_3.205Ma','xqbwa':'PI','xqbwb':'MP_3.205Ma'}
gm_d18o = {'xpsia':'0permil','xpsib':'-0.3permil','xqbwa':'0permil','xqbwb':'-0.3permil'}

HadCM3='y'
exptname='xqbwa'
startyear=3970  # can't start before year 12 because we aren't outputting d18o
nyears=30
plt.figure(figureno)
get_d18o_and_temperature()


# difference the two files
#diff_two_files()

#test plot sea surface d18o for comparison with data
#plot_ss_d18o()


sys.exit(0)

####

