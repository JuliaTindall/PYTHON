#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Thu Mar 18 14:13:50 2019

#@author: earjcti1
#
# This program will extract the fields needed for the PlioMIP2 database
#  It will extract the monthly averages from
# steve hunters processed data 
#
#
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
import sys
import warnings

warnings.filterwarnings("ignore")

def simplify_cube(cube):
    """
    gets cube and makes sure dimensions are longitude, latitude surface and
    t. 
    """    
    for coord in cube.coords():
        if coord.var_name == 'level275':
            coord.var_name = 'surface'
    
    cube.coord('surface').points = 0.0
    cube.coord('surface').units = 'm'
    cube.coord('surface').attributes = None
    
    cube.data = np.where(cube.data > 1.0E10, 0., cube.data)
    return cube

def get_evap(filename_):
    """
    will add up all the fluxes that make evaporation and returns
    total evaporation within a cube
    The fluxes are:
      evaporation from canopy
      evaporation from sea
      transpiration (this is exactly the same as evaporation from soil
                     I have checked some examples)
      sublim from surface
    """

    varnames_sec = ["EVAPORATION FROM SEA (GBM)   KG/M2/S",
                "TRANSPIRATION RATE           KG/M2/S"]
    
    varnames_ts = ["EVAP FROM CANOPY - AMOUNT   KG/M2/TS",
                "SUBLIM. FROM SURFACE (GBM)  KG/M2/TS"]
    
    
    for i, var in enumerate(varnames_sec):
        cube = iris.load_cube(filename_,var)
        cube = simplify_cube(cube)
        if i == 0:
            cubetot = cube
        else:
            cubetot = cubetot + cube
        
        
    for i, var in enumerate(varnames_ts):
        cube = iris.load_cube(filename_,var)
        cube = simplify_cube(cube)
        cube.data = cube.data / (30. * 60.)
        cube.units = 'kg m-2 s-1'
        cubetot = cubetot + cube
     
    return cubetot

#####################################
def extract_fields(filestart,expt,filetype,extra,startyear,endyear,timeperiod,
                   fileoutstart,varnamein,varnameout):

    # load required cubes
    #cubes=iris.load(filename)
    #print(cubes)
    #sys.exit(0)
    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    yearextra = {'0':'00', '1':'01','2':'02', '3':'03','4':'04', '5':'05',
                 '6':'06', '7':'07','8':'08', '9':'09','a':'10', 'b':'11',
                 'c':'12', 'd':'13','e':'14', 'f':'15','g':'16', 'h':'17',
                 'i':'18', 'j':'19','k':'20', 'l':'21','m':'22', 'n':'23',
                 'o':'24', 'p':'25','q':'26', 'r':'27','s':'28', 't':'29',
                 'u':'30', 'v':'31','w':'32', 'x':'33','y':'34', 'z':'35'}



    
    # loop over years
    
   
    for year in range(startyear,endyear):
        allcubes=CubeList([])
        stringyear=np.str(year).zfill(2)
        for mon in range(0,len(monthnames)):
            print(mon)
            if format == '#':
                filename=(filestart + expt + filetype + 
                      extra + stringyear + monthnames[mon] + '+.nc')
            else:
                print(filestart,expt,filetype, extra, stringyear,monthnames[mon])
                filename=(filestart + expt + filetype + 
                      extra + stringyear + monthnames[mon] + '.nc')
         
            if varnamein == 'TOTAL EVAPORATION':
                cube = get_evap(filename)
            else:
                print(filename,varnamein)
                cube=iris.load_cube(filename,varnamein)
                cube.coord('t').points=mon+1
            allcubes.append(cube)
            if mon==11:
                print(allcubes)
            
        #make sure the metadata on all cubes are the same
        iris.util.equalise_attributes(allcubes)
        catcube=allcubes.concatenate()[0]
        
        # if precipitation convert to mm/day
        if varnameout=='TotalPrecipitation':
            catcube.data=catcube.data * 60.*60.*24.
            catcube.long_name='TOTAL PRECIPITATION RATE    MM/DAY'
            catcube.units='mm/day'
       
        if varnameout=='evap':
            catcube.data=catcube.data * 60.*60.*24.
            catcube.long_name='TOTAL EVAPORATION    MM/DAY'
            catcube.units='mm/day'
      
        if varnameout=='so':
            catcube.data=(catcube.data * 1000.) + 35.
            catcube.long_name='SALMINTY (OCEAN)     (PSU)'
            catcube.units='ppt'
      
            
        # if temperature convert to Celcius
        print(varnameout)
        if varnameout=='SurfaceTemperature':
            print('convert to celcius')
            catcube.convert_units('celsius')
       
    
#        stringyear=np.str(year).zfill(2)
        stringyear=np.str(year).zfill(3)
        print(timeperiod,extra,stringyear)
        #fileout=fileoutstart+varnameout+'/'+timeperiod+'.'+varnameout+'.'+yearextra.get(extra) + stringyear+'.nc'        
        fileout=fileoutstart+varnameout+'/'+timeperiod+'.'+varnameout+'.'+ stringyear+'.nc'        
        iris.save(catcube,fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
     
    
   

##########################################################
# main program

# this is regridding where all results are in a single file
# create a dictionary with the long field names in and the field names we want
# we are also using dictionaries so that we only have to change timeperiod name
# when rerunning
            
fileextra = {
		"SURFACE TEMPERATURE AFTER TIMESTEP" : "a@pd",
        "TOTAL PRECIPITATION RATE     KG/M2/S" : "a#pd",
        "TEMPERATURE AT 1.5M": "a#pd",
        "OCN TOP-LEVEL TEMPERATURE          K" : "o@pf",
        "AICE : ICE CONCENTRATION" : "o@pf", 
                "TOTAL CLOUD AMOUNT - RANDOM OVERLAP" : "a#pd",
                "TEMPERATURE (OCEAN)  DEG.C" : "o@pf",
		"U COMPNT OF WIND ON PRESSURE LEVELS" : "a@pc",
		"V COMPNT OF WIND ON PRESSURE LEVELS" : "a@pc", 
        	"NET DOWN SURFACE SW FLUX: SW TS ONLY" : "a@pd",
		"NET DOWN SURFACE LW RAD FLUX" : "a@pd",
        	"SURFACE & B.LAYER HEAT FLUXES   W/M2" : "a@pd",
		"SURFACE LATENT HEAT FLUX        W/M2" : "a@pd",
        "OMEGA ON PRESSURE LEVELS" : "a@pc",
		"SPECIF HUM;P LEVS;U GRID.  USE MACRO" : "a@pc",
		"GEOPOTENTIAL HEIGHT: PRESSURE LEVELS" : "a@pc",
		"TEMPERATURE ON PRESSURE LEVELS" : "a@pc",
		"SURFACE LATENT HEAT FLUX        W/M2" : "a@pd",
        "DOWNWARD LW RAD FLUX: SURFACE" : "a@pd",
        "TOTAL DOWNWARD SURFACE SW FLUX" : "a@pd",
        "INCOMING SW RAD FLUX (TOA): ALL TSS" : "a@pd",
        "OUTGOING SW RAD FLUX (TOA)" : "a@pd",
        "OUTGOING LW RAD FLUX (TOA)" : "a@pd",
        "SURFACE & B.LAYER HEAT FLUXES   W/M2" : "a@pd",
		"PRESSURE AT MEAN SEA LEVEL" : "a@pd",
		"PSTAR AFTER TIMESTEP" : "a@pd",
		"TOTAL EVAPORATION" : "a@pd",
		"X-COMP OF SURF & BL WIND STRESS N/M2" : "a@pd",
                "Y-COMP OF SURF & BL WIND STRESS N/M2" : "a@pd",
		"TOTAL OCEAN U-VELOCITY      CM S**-1" : "o@pf",
		"TOTAL OCEAN V-VELOCITY      CM S**-1" : "o@pf",
		"VERT.VEL. ON OCEAN HALF LEVELS  CM/S" : "o@pf",
                "POTENTIAL TEMPERATURE (OCEAN)  DEG.C" : "o@pf",
                "SALINITY (OCEAN)       (PSU-35)/1000" : "o@pf"
	}

shortname = {
		"SURFACE TEMPERATURE AFTER TIMESTEP" : "SurfaceTemperature",
        "TOTAL PRECIPITATION RATE     KG/M2/S" : "TotalPrecipitation",
        "TEMPERATURE AT 1.5M": "NearSurfaceTemperature",
        "TEMPERATURE (OCEAN)  DEG.C":"OceanTemp",
        "TOTAL CLOUD AMOUNT - RANDOM OVERLAP" : "totcloud",
        "OCN TOP-LEVEL TEMPERATURE          K" : "SST",
        "AICE : ICE CONCENTRATION" : "SeaIceConc", 
		"U COMPNT OF WIND ON PRESSURE LEVELS" : "ua",
		"V COMPNT OF WIND ON PRESSURE LEVELS" : "va",
        	"NET DOWN SURFACE SW FLUX: SW TS ONLY" : "fsns",
		"NET DOWN SURFACE LW RAD FLUX" : "flns",
        "OMEGA ON PRESSURE LEVELS" : "wap",
		"SPECIF HUM;P LEVS;U GRID.  USE MACRO" : "spechumid",
		"GEOPOTENTIAL HEIGHT: PRESSURE LEVELS" : "zg",
		"TEMPERATURE ON PRESSURE LEVELS" : "ta",
		"SURFACE LATENT HEAT FLUX        W/M2" : "hfls",
        "DOWNWARD LW RAD FLUX: SURFACE" : "rlds",
        "TOTAL DOWNWARD SURFACE SW FLUX" : "rsds",
        "INCOMING SW RAD FLUX (TOA): ALL TSS" : "rsdt",
        "OUTGOING SW RAD FLUX (TOA)" : "rsut",
        "OUTGOING LW RAD FLUX (TOA)" : "rlut",
        "SURFACE & B.LAYER HEAT FLUXES   W/M2" : "surfheatflux",
        "CLEAR-SKY (II) UPWARD LW FLUX (TOA)" : "rlutcs",
                "CLEAR-SKY (II) UPWARD SW FLUX (TOA)" : "rsutcs",
		"PRESSURE AT MEAN SEA LEVEL" : "MSLP",
		"PSTAR AFTER TIMESTEP" : "ps",
		"TOTAL EVAPORATION" : "evap",
		"X-COMP OF SURF & BL WIND STRESS N/M2" : "tauu",
                "Y-COMP OF SURF & BL WIND STRESS N/M2" : "tauv",
		"TOTAL OCEAN U-VELOCITY      CM S**-1" : "uo",
		"TOTAL OCEAN V-VELOCITY      CM S**-1" : "vo",
		"VERT.VEL. ON OCEAN HALF LEVELS  CM/S" : "wo",
                "POTENTIAL TEMPERATURE (OCEAN)  DEG.C" : "thetao",
                "SALINITY (OCEAN)       (PSU-35)/1000" : "so"
	}


exptname = {
        "e280" : "tenvo",
        "e400" : "tenvq",
        "e560":"tenvs",
        "eoi400" : "tenvj",
        "eoi350" : "tenvk",
        "eoi450" : "tenvl",
        "eoi280" : "tenvm",
        "e280_corr" : "xozzz"
      
}

extraname = {
        "e280" : "t",
        "e400" : "t",
        "eoi400" : "o",
        "eoi350" : "o",
        "eoi450" : "o",
        "eoi280" : "o",
        "e280_corr" : "s",
        "e560" : "v",
        "xozza" : "p",
        "xozzc" : "o",
        "xozzd" : "o",
        "xozze" : "o",
        "xozzf" : "o",
        "tenvl" : "o",
        "tenvk" : "o",
        "tenvm" : "o",
        "xozzz" : "p",
        "xozzb" : "p",  # anything before 14.3.2023 used "o" after that 
                        # it should be 'p'
                        # in outputfilename
        "xpkma" : "00000",
        "xpkmb" : "00000",
        "xpkmc" : "00000"}

#fieldname=["Y-COMP OF SURF & BL WIND STRESS N/M2","X-COMP OF SURF & BL WIND STRESS N/M2" 
#           ]
fieldname = [#"SURFACE TEMPERATURE AFTER TIMESTEP",
#         "TOTAL PRECIPITATION RATE     KG/M2/S",
#        "TEMPERATURE AT 1.5M",
       #  "OCN TOP-LEVEL TEMPERATURE          K",
      #  "TEMPERATURE (OCEAN)  DEG.C",
       #  "AICE : ICE CONCENTRATION", 
#          "TOTAL CLOUD AMOUNT - RANDOM OVERLAP",
       
	#	"U COMPNT OF WIND ON PRESSURE LEVELS",
#		"V COMPNT OF WIND ON PRESSURE LEVELS",
#        	"NET DOWN SURFACE SW FLUX: SW TS ONLY",
#		"NET DOWN SURFACE LW RAD FLUX",
#        	"SURFACE & B.LAYER HEAT FLUXES   W/M2",
#		"SURFACE LATENT HEAT FLUX        W/M2",
#        "OMEGA ON PRESSURE LEVELS",
#		"SPECIF HUM;P LEVS;U GRID.  USE MACRO",
#		"GEOPOTENTIAL HEIGHT: PRESSURE LEVELS",
#		"TEMPERATURE ON PRESSURE LEVELS",
#		"SURFACE LATENT HEAT FLUX        W/M2",
#        "DOWNWARD LW RAD FLUX: SURFACE",
#        "TOTAL DOWNWARD SURFACE SW FLUX",
#        "INCOMING SW RAD FLUX (TOA): ALL TSS",
#        "OUTGOING SW RAD FLUX (TOA)",
#        "OUTGOING LW RAD FLUX (TOA)",
#         "CLEAR-SKY (II) UPWARD SW FLUX (TOA)",
#         "CLEAR-SKY (II) UPWARD LW FLUX (TOA)",
#        "SURFACE & B.LAYER HEAT FLUXES   W/M2",
#		"PRESSURE AT MEAN SEA LEVEL",
#		"PSTAR AFTER TIMESTEP",
		"TOTAL EVAPORATION",
#		"X-COMP OF SURF & BL WIND STRESS N/M2",
#                "Y-COMP OF SURF & BL WIND STRESS N/M2",
#		"TOTAL OCEAN U-VELOCITY      CM S**-1",
#		"TOTAL OCEAN V-VELOCITY      CM S**-1",
#		"VERT.VEL. ON OCEAN HALF LEVELS  CM/S",
          "POTENTIAL TEMPERATURE (OCEAN)  DEG.C",
            "SALINITY (OCEAN)       (PSU-35)/1000" 
	]
	       
#fieldname=["V COMPNT OF WIND ON PRESSURE LEVELS"]

linux_win='l'
startyear=0
endyear=100
timeperiod='tenvm'   
expt=exptname.get(timeperiod,timeperiod)
extra=extraname.get(timeperiod,'00000')
format = '@'  # is the filename xxxxx#pd or xxxxx@pd

if linux_win=='w':
    filestart='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\HadCM3\\'+exptname.get(timeperiod)+'/'
    fileoutstart='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\HadCM3_UPLOAD\\'+timeperiod+'/'
else:
    filestart='/nfs/hera1/earjcti/um/'+expt+'/pd/'
    fileoutstart='/nfs/hera1/earjcti/um/'+expt + '/'
    #fileoutstart='/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/'+timeperiod+'/'

for i in range(0,len(fieldname)):
    varnamein=fieldname[i]
    varnameout=shortname.get(varnamein)
    filetype=fileextra.get(varnamein,'a'+ format+ 'pd')

    extract_fields(filestart,expt,filetype,extra,startyear,endyear,timeperiod,fileoutstart,varnamein,varnameout)

#sys.exit(0)
