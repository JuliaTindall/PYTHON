#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created July 2023 by Julia

This program will sort out the vegetation conditions for the MP.  
"""

import numpy as np
import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib
import matplotlib.pyplot as plt
from iris.cube import CubeList
#import matplotlib.ticker as mticker
#import netCDF4
#from datetime import date
#import cmocean to use topo colormap
#import sys


def compare_mod_HadCM3_HadGEM2():
    """
    looks to see how the modern veg differs in HadCM3 and HadGEM2
    """

    pihg2_cube = iris.util.squeeze(iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/qrparm.veg.frac_1860.nc'))
    pseudo_dim = pihg2_cube.coord('pseudo')
   
    hcm3_cubes = CubeList([])
    for lev in range(1,10):
        cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/preind2/qrfrac_hcm3_lev' + str(lev) + '.nc')
        cube.coord('Surface').rename('pseudo')
        cube.coord('pseudo').points = lev
        hcm3_cubes.append(cube)
    iris.util.equalise_attributes(hcm3_cubes)
    pihc3_cube_lr = hcm3_cubes.concatenate_cube()
   
    # regrid hadcm3 cube to hg2 grid
    pihc3_cube = iris.util.squeeze(pihc3_cube_lr.regrid(pihg2_cube,iris.analysis.Linear()))
    #pihc3_cube.add_dim_coord(pseudo_dim,0)

    print(pihc3_cube)
    print(pihg2_cube)

    diffcube = pihg2_cube - pihc3_cube
    for lev in range(0,9):
        fig = plt.figure(figsize=[11.0,4.0])
        plt.subplot(1,3,1)
        qplt.contourf(pihc3_cube[lev,:,:],levels=np.arange(0,1.1,0.1))
        plt.title('hadcm3 lev' + str(lev))
        plt.subplot(1,3,2)
        qplt.contourf(pihg2_cube[lev,:,:],levels=np.arange(0,1.1,0.1))
        plt.title('hadgem2 lev' + str(lev))
        plt.subplot(1,3,3)
        qplt.contourf(diffcube[lev,:,:],levels=np.arange(-0.5,0.6,0.1),
                      extend='both')
        plt.title('diff')
        plt.savefig('vegtype_' + str(lev) + '.eps')
        plt.close()
                             
#########################################################################
def compare_megabiomes_pi_plio():

    vegnames = {1:'Tropical Forest',2:'Warm-Temperate forest',
                3:'Savannah',4:'Grassland',5:'Desert',6:'temperate forest',
                7:'Boreal Forest',8:'Tundra',9:'Dry Tundra',28: 'Landice'}
    vegtype = [1,2,3,4,5,6,7,8,9,28]
    hg2_cube = iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/qrparm.mask.nc','LAND MASK (No halo) (LAND=TRUE)')
    pi_triffid_cube = iris.util.squeeze(iris.load_cube('/nfs/hera1/earjcti/UM_ANCIL/HadGEM/00.0k/qrparm.veg.frac_1860.nc'))
   

    # for this we will regrid onto HadGEM2 grid using nearest neighbour

    start = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/'
    pi_mb_cube_origgrid = iris.load_cube(start + 'Modern_std/Modern_std/' + 
                                         'Modern_std_mbiome_v1.0.nc')
    pi_mb_cube = pi_mb_cube_origgrid.regrid(hg2_cube,iris.analysis.Nearest())
  
    pi_data = pi_mb_cube.data
    pi_data_tiled = np.tile(pi_data,[9,1,1]) # used for getting triffid
      
  
   
    plio_mb_cube_origgrid = iris.load_cube(start + 'Plio_enh/Plio_enh/'+
                                  'Plio_enh_mbiome_v1.0.nc')
    plio_mb_cube = plio_mb_cube_origgrid.regrid(hg2_cube,iris.analysis.Nearest())
 
    plio_data = plio_mb_cube.data

    for i in vegtype:
        pi_thisbiome = np.where(pi_data == i, 1.0, 0.0)
        plio_thisbiome = np.where(plio_data == i, 1.0, 0.0)
        
        pi_thisbiome_cube = pi_mb_cube.copy(data = pi_thisbiome)
        plio_thisbiome_cube = plio_mb_cube.copy(data = plio_thisbiome)

        diffcube = pi_mb_cube.copy(data=plio_thisbiome - pi_thisbiome)

        fig = plt.figure(figsize=[11.0,4.0])
        plt.subplot(1,3,1)
        qplt.contourf(pi_thisbiome_cube)
        plt.gca().coastlines()
        plt.title('pi' + vegnames.get(i))
        plt.subplot(1,3,2)
        qplt.contourf(plio_thisbiome_cube)
        plt.gca().coastlines()
        plt.title('plio' + vegnames.get(i))
        plt.subplot(1,3,3)
        qplt.contourf(diffcube)
        plt.gca().coastlines()
        plt.title('diff')
        plt.savefig('mbiome_' + vegnames.get(i) + '.eps')
        plt.close()
        

        # extract the TRIFFID TYPES FOR THIS BIOME IN THE PREINDUSTRIAL
        triffid_thismb = np.ma.masked_where(pi_data_tiled != i, pi_triffid_cube.data)
        triffid_thismb = np.where(triffid_thismb > 1.0E10, 99999,triffid_thismb)
        #triffid_thismb = np.ma.masked_where(pi_data != i, 1.0)
        triffid_thismb_cube = pi_triffid_cube.copy(data=triffid_thismb)
        iris.save([triffid_thismb_cube,pi_thisbiome_cube,plio_thisbiome_cube,
                  diffcube],vegnames.get(i) + '.nc',fill_value = 99999)
    
    sys.exit(0)
    
##########################################################################
# setup what to change

#######################################
# compare pi veg.  HadCM3 vs HadGEM2

#compare_mod_HadCM3_HadGEM2()

#####################################
#  compare preindustrial and pliocene megabiomes

compare_megabiomes_pi_plio()

# next add in lakes



# this is what needs changing (lat,lon)
#gridbox_chg = [[12.5,273.75],[10.0,277.5],[7.5,281.25]]
#gridbox_lsm= {'[12.5,273.75]' : 0., '[10.0,277.5]': 0., '[7.5,281.25]': 0.}
#gridbox_missing= {'[12.5,273.75]' : 2.0E20, '[10.0,277.5]': 2.0E20, '[7.5,281.25]': 2.0E20}
#gridbox_missing_alt= {'[12.5,273.75]' : -1.0737E+9, '[10.0,277.5]': -1.0737E+9, '[7.5,281.25]': -1.0737E+9}
#gridbox_zero= {'[12.5,273.75]' : 0.0, '[10.0,277.5]': 0.0, '[7.5,281.25]': 0.0}


# check all the files have the same land sea mask
#make_lsm_EP()
#make_other_ancil_EP()
#check_lsm('EP_mbiome_v1.0.nc','veg', 'm')
#check_lsm('EP_soil_v1.0.nc','soil', 'm')
#check_lsm('EP_lake_v1.0.nc','lake', 'm')
