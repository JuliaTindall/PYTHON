#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 18 14:13:50 2019

#@author: earjcti1
#
#  This program will extract fileds and put in a timeseries file
#  this is what we will upload to globus for PlioMIP3
#  it is daily variables so we need core only
#

import os
import numpy as np
import scipy as sp
#import cf
import iris
from iris.cube import CubeList
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import iris.coord_categorisation
import sys
import warnings

     
def monthly_data_single_level(field,shortname):
    """
    this will get the database averages for the monthly data on a single level
    """

    allcubes = iris.load('/uolstore/home/users/earjcti/hera1/um/' + expt + '/pcpd/' + expt + 'a#pd000003[8-9]*',field)
    iris.util.equalise_attributes(allcubes)
    cubes = allcubes.concatenate_cube()
    cubes.var_name = shortname

    

    #if shortname == 'clt':
        # convert from cloud fraction to percent
    #    cubes.data = cubes.data * 100.
    #    cubes.units = "%"

    #    print(cubes)

    return cubes
   


def monthly_data_all_levels(field,shortname):
    """
    this will get the database averages for the monthly data on all levels
    (we will keep them all)
    """

    allcubes = iris.load('/uolstore/home/users/earjcti/hera1/um/' + expt + '/pcpd/' + expt + 'a#pc000003[8-9]*',field)
    #print('allcubes',allcubes,field)
    iris.util.equalise_attributes(allcubes)
    cubes = allcubes.concatenate_cube()
    cubes.var_name = shortname

    if shortname == 'ta':
        # convert to degC
        cubes.convert_units("degC")
        print(cubes)


    return cubes

def monocn_data_single_level(field,shortname):
    """
    this will get the database averages for the monthly data on a single level
    for the ocean
    """

    allcubes = iris.load('/uolstore/home/users/earjcti/hera1/um/' + expt + '/pf/' + expt + 'o#pf000003[8-9]*',field)
    print('allcubes',allcubes,field)
    iris.util.equalise_attributes(allcubes)
    cubes = allcubes.concatenate_cube()
    cubes.var_name = shortname

    cubes.data = np.ma.where(cubes.data.mask == 1.0, -999.999, cubes.data)

    if shortname == 'tos':
        cubes.long_name = "OCN TOP-LEVEL TEMPERATURE"
        cubes.attributes["title"] = cubes.long_name

    if shortname == 'siconc':
        cubes.units = "%"
        cubes.data = cubes.data * 100.

    #    print(cubes)

    return cubes
   


##########################################################
# main program

# this is regridding where all results are in a single file
# create a dictionary with the long field names in and the field names we want
# we are also using dictionaries so that we only have to change timeperiod name
# when rerunning
            

expt = 'xqbwg'


fields = [ "SNOW AMOUNT AFTER TIMESTEP     KG/M2",
          "SURFACE RUNOFF RATE          KG/M2/S",
          "SUB-SURFACE RUNOFF RATE      KG/M2/S",
          "SOIL MOISTURE CONTENT"
           "INCOMING SW RAD FLUX (TOA): ALL TSS",
           "OUTGOING SW RAD FLUX (TOA)",
           "OUTGOING LW RAD FLUX (TOA)",
           "CLEAR-SKY (II) UPWARD LW FLUX (TOA)",
           "DOWNWARD LW RAD FLUX: SURFACE",
           "TOTAL DOWNWARD SURFACE SW FLUX",
           "NET DOWN SURFACE LW RAD FLUX",
           "NET DOWN SURFACE SW FLUX: SW TS ONLY",
           "CLEAR-SKY (II) UP SURFACE SW FLUX",
           "CLEAR-SKY (II) DOWN SURFACE SW FLUX",
           "CLEAR-SKY (II) DOWN SURFACE LW FLUX",
           "CLEAR-SKY (II) UPWARD SW FLUX (TOA)",
           "SURFACE LATENT HEAT FLUX        W/M2",
           "SURFACE & B.LAYER HEAT FLUXES   W/M2"]

#expt = 'xqfmg'


###############################
### DICTIONARIES FOR NAMING
# expt name
alt_expt = {'xqfmg': 'F_EP280',
            'xqfmh': 'F_EP',
            'xqfme' : 'F_LP280',
            'xqfmf' : 'F_LP',
            'xqfmb' : 'F_PI280',
            'xqfmc' : 'F_PI400',
            'xqfmd' : 'F_PI490',
            'xqbwc' : 'PI', 'xqbwd' : 'LP',
            'xqbwg':'EP',
            'xqbwn':'high_NH',
            'xqbwo':'high_SH'}

# CMIP name  ; this is from fernandas spreadsheet

cmip_name = {"SNOW AMOUNT AFTER TIMESTEP     KG/M2" : "snw",
             "SURFACE RUNOFF RATE          KG/M2/S" : "mrros",
             "SUB-SURFACE RUNOFF RATE      KG/M2/S": "mrrob",
             "SOIL MOISTURE CONTENT" : "mrso",
             "INCOMING SW RAD FLUX (TOA): ALL TSS" : "rsdt",
             "OUTGOING SW RAD FLUX (TOA)": "rsut",
             "OUTGOING LW RAD FLUX (TOA)":"rlut",
             "CLEAR-SKY (II) UPWARD LW FLUX (TOA)" : "rlutcs",
             "DOWNWARD LW RAD FLUX: SURFACE":"rlds",
             "TOTAL DOWNWARD SURFACE SW FLUX":"rsds",
             "NET DOWN SURFACE LW RAD FLUX":"rlntds",
             "NET DOWN SURFACE SW FLUX: SW TS ONLY":"rsntds",
             "CLEAR-SKY (II) UP SURFACE SW FLUX":"rsuscs",
             "CLEAR-SKY (II) DOWN SURFACE SW FLUX":"rsdscs",
             "CLEAR-SKY (II) DOWN SURFACE LW FLUX":"rldscs",
             "CLEAR-SKY (II) UPWARD SW FLUX (TOA)":"rsutcs",
             "SURFACE LATENT HEAT FLUX        W/M2":"hfls",
             "SURFACE & B.LAYER HEAT FLUXES   W/M2":"hfss"}

levels_req = {'ua':'y','va':'y','ta':'y','wa':'y','zg':'y','hus':'y'}
ocean_req = {'tos':'y', 'siconc':'y','sithick':'y'}
            

for field in fields:
    shortname = cmip_name.get(field)
    ocn = ocean_req.get(shortname,'n')
    multlev = levels_req.get(shortname,'n')
    
    if multlev == 'n' and ocn == 'n':  # single level only atmosphere
        cubes = monthly_data_single_level(field,shortname)

    if multlev == 'y' and ocn == 'n':  # multiple level atm
        cubes = monthly_data_all_levels(field,shortname)

    if  multlev == 'n' and ocn == 'y': #single level ocean
        cubes = monocn_data_single_level(field,shortname)

    # find climatological mean
    iris.coord_categorisation.add_month_number(cubes,  't',  name = 'month')
    iris.coord_categorisation.add_year(cubes,  't',  name = 'year')
    meanmonthcube = cubes.aggregated_by('month', iris.analysis.MEAN)
    iris.util.promote_aux_coord_to_dim_coord(meanmonthcube,'month')
  
  
        
    fileout = ('/uolstore/Research/a/hera1/earjcti/um/' + expt + '/globus/HadCM3_')
    fileout = (fileout + alt_expt.get(expt) + '_' + expt + '_' + 
               cmip_name.get(field) + '_monthly_clim_3800_4000.nc')
    if ocn == 'y':
        iris.save(meanmonthcube,fileout,fill_value=-999.999)
    else:
        iris.save(meanmonthcube,fileout)
    sys.exit(0)



