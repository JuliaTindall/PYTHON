#!/usr/bin/env python2.7
#NAME
#   avgAMOC from pk files (also does GMOC PMOC IMOC)
#PURPOSE
#writest to a netcdf file for reading in elsewehere
#  b
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
    
    cubelist = iris.load(filename)
    #for cube in cubelist:
    #    print(cube.var_name, cube.long_name)
    #sys.exit(0)
    Vcube = cubelist.extract('field704')[0]
    Wcube = cubelist.extract('W')[0]


    return (Vcube,Wcube)

#######################################################
def get_data():
    """
    as described
    """

    # arrays for storing mean salratioeratures
    cubelists = CubeList([])
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



################################
# main program

# annual mean
figureno=0

timeperiod = {'xpsia':'PI','xpsib':'MP_3.205Ma','xqbwa':'PI','xqbwb':'MP_3.205Ma','xqbwj':'LP490'}

HadCM3='y'
exptname='xqbwj'
centuary='39'  # currently doing one centuary
basins = ['Atlantic','Pacific','Indian','Global']

mean_cubelist = CubeList([])

for basin in basins:

    files=('/uolstore/Research/a/hera1/earjcti/um/'+exptname+'/pk2/'+
           exptname+'o#pk00000' + centuary + '??c1+.nc')
    field = 'Meridional Overturning Stream Function (' + basin + ')'
    allcubes = iris.load(files,field)

    conccube = allcubes.concatenate_cube()
    meancube = conccube.collapsed('time',iris.analysis.MEAN)

    mean_cubelist.append(meancube)

fileout = ('/uolstore/Research/a/hera1/earjcti/um/'+exptname+'/MOC/'+
           exptname+'pkmean_' + centuary + '00_' + centuary + '99.nc')
iris.save(mean_cubelist,fileout)


# difference the two files
#diff_two_files()

#test plot sea surface d18o for comparison with data
#plot_ss_d18o()


sys.exit(0)

####

