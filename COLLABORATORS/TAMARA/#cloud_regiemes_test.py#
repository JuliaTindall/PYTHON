#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29.03.2022 by Julia

We are trying to plot cloud regiemes based on Williams and Webb 2008.
This is an experimental program
1. plot mean cloud albedo, 
2. cloud top pressure
3. total cloud cover
"""
import numpy as np
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import sys

  
def get_cloud_regieme_diagnostics():
    """
    for this test I think we need the following diagnostics
    1. total cloud cover - TRY TOTAL CLOUD AMOUNT - RANDOM OVERLAP
    2. cloud top pressure
    3. mean cloud albedo
    """
    tempcube = iris.util.squeeze(iris.load_cube(
            TESTFILE, 'TOTAL CLOUD AMOUNT - RANDOM OVERLAP'))
    TCCcube=tempcube.collapsed('t_1',iris.analysis.MEAN)

    tempcube = iris.util.squeeze(iris.load_cube(
            TESTFILE, 'INCOMING SW RAD FLUX (TOA): ALL TSS'))
    SW_IN_TOP_cube = tempcube.collapsed('t_1',iris.analysis.MEAN)
    #mask = np.where(SW_IN_TOP_cube.data == 0, 1.0, 0.0)
    #SW_IN_TOP_cube.data.mask = mask
   
    tempcube = iris.util.squeeze(iris.load_cube(
            TESTFILE, 'TOTAL DOWNWARD SURFACE SW FLUX'))
    SW_DOWN_SURF_cube = tempcube.collapsed('t',iris.analysis.MEAN)
    
    cloud_albedo_cube_data = 1.0 - (SW_DOWN_SURF_cube.data / SW_IN_TOP_cube.data)
    cloud_albedo_cube = SW_DOWN_SURF_cube.copy(data=cloud_albedo_cube_data)
    cloud_albedo_cube.long_name = 'albedo'   
    # mask where there is not much incoming sw radiation
    cloud_albedo_cube.data.mask = np.where(SW_IN_TOP_cube.data<20.0, 1.0, 0.0)
    

  
    tempcube = iris.util.squeeze(iris.load_cube(
            TESTFILE,'TOTAL CLOUD TOP HEIGHT (KFT)'))
    cloud_top_pressure_cube = tempcube.collapsed('t_1',iris.analysis.MEAN)
    #TCCcube.data.mask = mask
    #cloud_albedo_cube.data.mask = mask
    #cloud_top_pressure_cube.data.mask = mask

    # try and find alternative albedo as above but remove effects of surface albedo
    # I think this is (1.0 - downward solar surface / upwards surface cs)
    
    tempcube = iris.util.squeeze(iris.load_cube(
            TESTFILE, 'CLEAR-SKY (II) UPWARD SW FLUX (TOA)'))
    SW_UPCS_TOA_cube = tempcube.collapsed('t',iris.analysis.MEAN)

    alt_alb_cube_data = 1.0 - (SW_DOWN_SURF_cube.data / SW_UPCS_TOA_cube.data)
    alt_alb_cube = SW_DOWN_SURF_cube.copy(data = alt_alb_cube_data)
    alt_alb_cube.long_name = 'alternative albedo'
    alt_alb_cube.data = np.where(alt_alb_cube.data<0.0, 0.0, alt_alb_cube.data)
    alt_alb_cube.data = np.where(alt_alb_cube.data>1.0, 1.0, alt_alb_cube.data)
    
    plt.subplot(221)
    qplt.contourf(cloud_top_pressure_cube, levels=np.arange(10000, 90000, 10000), extend='both')
    plt.gca().coastlines()
    plt.subplot(222)
    qplt.contourf(TCCcube, levels=np.arange(0.20, 1.05, 0.05), extend='both')
    plt.gca().coastlines()
    plt.subplot(223)
    qplt.contourf(cloud_albedo_cube, levels=np.arange(0.20, 0.60, 0.05), extend='both')
    plt.gca().coastlines()
    plt.subplot(224)
    qplt.contourf(alt_alb_cube, levels=np.arange(0.2, 0.6, 0.05), extend='both')
    plt.gca().coastlines()
    #plt.show()
    

   
    return alt_alb_cube, cloud_top_pressure_cube, TCCcube,cloud_albedo_cube
                                 
  
##########################################################
def get_tropical_regiemes(albedo_cube, CTP_cube, TCC_cube,orig_alb_cld):
    """
    estimates the cloud regieme for each tropical gridpoint
    """
    REGIEME_NAMES_TR = {0:"Shallow cumulus", 1:"Congestus",
                 2:"Thin cirrus", 3:"Stratocu./Cu. Transition",
                 4:"Anvil cirrus", 5:"Deep Convection",
                 6:"Stratocumulus"}
    REGIEME_CHARACTERISTICS = {0:[0.261, 0.652, 0.314], 
                           1:[0.339, 0.483, 0.813],
                           2:[0.211, 0.356, 0.740],
                           3:[0.338, 0.784, 0.640],
                           4:[0.532, 0.285,0.979],
                           5:[0.446, 0.722, 0.824]}
   
    # upper and lower value of bin
    bin_lower = [0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0] 
    bin_upper = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1000000]
    binnamesp = []
    binnames = []
    for i, low in enumerate(bin_lower):
        binnamesp.append(np.str(np.int(low * 1000.)))
        binnames.append(np.str(np.int(low)))
    
   # binnames.append("{:.1f}".format(bin_upper[-1]))
    pressure_bins = np.zeros(11)
    TCC_bins = np.zeros(11)
    albedo_bins = np.zeros(11)
    alt_albedo_bins = np.zeros(11)

    for j, lat in enumerate(albedo_cube.coord('latitude').points):
        if -20 <= lat <= 20:
            for i, lon in enumerate(albedo_cube.coord('longitude').points):
                # get pressure
                pressure=CTP_cube.data[j,i] / 100000.
                cloud_frac = TCC_cube.data[j,i] 
                if cloud_frac == 0.0:
                    albedo=0.0
                else:
                    albedo = albedo_cube.data[j,i] / cloud_frac

                o_albedo = orig_alb_cld.data[j,i]
                print("{:.3f}".format(pressure))
                for b, low in enumerate(bin_lower):
                    if low < pressure < bin_upper[b]:
                        pressure_bins[b] = pressure_bins[b] + 1
                    if low < cloud_frac < bin_upper[b]:
                        TCC_bins[b] = TCC_bins[b] + 1
                    if low < o_albedo < bin_upper[b]:
                        alt_albedo_bins[b] = alt_albedo_bins[b] + 1
             

                #print(lon,lat,"{:.2f}".format(albedo),"{:.2f}".format(o_albedo),
#                      "{:.3f}".format(pressure),"{:.2f}".format(cloud_frac))
               
                                        

    # plot bar charts of where everything is
    plt.subplot(3,1,1)                   
    plt.bar(binnamesp,pressure_bins,align = 'edge')
    plt.xlabel('Pressure (hPa)')
    plt.ylabel('Number')
    plt.title('Tropical cloud top pressure')
  
    plt.subplot(3,1,2)                   
    plt.bar(binnames,TCC_bins,align = 'edge')
    plt.xlabel('Fraction')
    plt.ylabel('Number')
    plt.title('Cloud Fraction (tropics)')
    #plt.show()
########################################################
# main program

#regiemes are defined in Williams and web characteristics are 
# cloud albedo (frac), cloud top pressure (hpa) / 1000, total cloud cover (frac)
# tropics
# extratropics
REGIEME_NAMES_ET = {0: "Shallow cumulus", 1:"Congestus (TR)",
                 2:"Stratocu./Cu. Transition (TR)", 3: "cirrus",
                 4:"Stratocumulus", 5:"Frontal",
                 6:"Thin Cirrus"}
# extratropics - ice covered
REGIEME_NAMES_ETice = {0: "Shallow cumulus", 1:"Stratocumulus",
                 2:"Thick mid-level", 3: "Frontal",
                 4:"Thin mid-level", 5:"Thin Cirrus"}

MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'

EXPT = 'tenvj'
#TESTFILE = '/nfs/hera1/earjcti/um/tenvo/pd/tenvoa@pdt99sp.nc'
TESTFILE = 'xozzaa@pap00dc.nc'
(albedo_cube, cloud_top_pressure_cube, TCCcube,orig_cloud_albedo) = get_cloud_regieme_diagnostics() 

iris.save([albedo_cube, cloud_top_pressure_cube, TCCcube,orig_cloud_albedo],'test.nc',
          fill_value = -999.999)


#topical_cloud_regiemes = get_tropical_regiemes(albedo_cube, cloud_top_pressure_cube, TCCcube,orig_cloud_albedo)


