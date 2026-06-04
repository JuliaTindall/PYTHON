#!/usr/bin/env python2.7
#NAME
#   avg_v_w_for_heather
#PURPOSE
#   Heather wants v and w values
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



def get_avg(year):
    """
    gets the average v and w for this year
    """
 
   
    yearuse = str(year).zfill(9)
    filename=('/uolstore/Research/a/hera1/earjcti/um/'+exptname+'/pg/'+
              exptname+'o#pg' + yearuse + 'c1+.nc')
    
    cubelist = iris.load(filename)
    #for cube in cubelist:
    #    print(cube.var_name, cube.long_name)
    #sys.exit(0)
    Vcube = cubelist.extract('field704')[0]
    Wcube = cubelist.extract('W')[0]


    return (Vcube,Wcube)

#######################################################
def get_v_and_w():
    """
    as described
    """

    # arrays for storing mean salratioeratures
    cubelist_v = CubeList([])
    cubelist_w = CubeList([])

    # obtain means for each year and store in the arrays
    for year in range(startyear,startyear+nyears):
        print(year)
        (v_cube, w_cube) = get_avg(year)
        cubelist_v.append(v_cube)
        cubelist_w.append(w_cube)

    
    # put into a single cube and average
    iris.util.equalise_attributes(cubelist_v)
    iris.util.equalise_attributes(cubelist_w)
    
    cubesV = cubelist_v.concatenate_cube()
    cubesW = cubelist_w.concatenate_cube()

    avgV_cube = cubesV.collapsed('t',iris.analysis.MEAN)
    avgW_cube = cubesW.collapsed('t',iris.analysis.MEAN)

    avgV_cube.data = np.ma.where(avgV_cube.data < -900.0,
                                 -99999.,avgV_cube.data)
    avgW_cube.data = np.ma.where(avgW_cube.data < -900.0,
                                 -99999.,avgW_cube.data)
   
    #iris.util.mask_cube(cubeavgd18o,np.where(cubeavg18o.data > 1.0E10))
    #iris.util.mask_cube(cubeavgd18o_adj,np.where(cubeavg18o_adj.data < -900.0))
    #iris.util.mask_cube(cubeavgtemp,np.where(cubeavg18o.data > 1.0E10))
    
    #print('j2',cubeavgd18o_adj.data)
   
    filename = (timeperiod.get(exptname) + '_' + exptname + '_V_and_W_'
                + str(startyear) + '_' + str(startyear+nyears) + '.nc')


    
    iris.save([avgV_cube, avgW_cube],filename,fill_value = -99999.)


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

HadCM3='y'
exptname='xqbwb'
startyear=3970  # can't start before year 12 because we aren't outputting d18o
nyears=30
plt.figure(figureno)
get_v_and_w()


# difference the two files
#diff_two_files()

#test plot sea surface d18o for comparison with data
#plot_ss_d18o()


sys.exit(0)

####

