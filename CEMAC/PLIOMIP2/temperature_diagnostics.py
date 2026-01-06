#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 29th


#@author: earjcti
#
# This program will read all the means from the regridded files and plot them


import os
import numpy as np
import scipy as sp
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
#from mpl_toolkits.basemap import Basemap, shiftgrid
#import Basemap
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
import cf_units as unit
from iris.experimental.equalise_cubes import equalise_attributes
import sys

def resort_coords(cube):
    """
    this will make all the dimensions of the cube match. 
    """
    
    for coord in cube.coords():        
        name=coord.standard_name
        if name !='latitude' and name!='longitude':
            if name==None:
                if coord.long_name==None:
                    cube.remove_coord(coord.var_name)
                else:
                    cube.remove_coord(coord.long_name)
            else:
                cube.remove_coord(name)
                
    for coord in cube.coords():   # now this will be longitude or latitude
        coord.points=coord.points.astype('float32') 
        coord.var_name = coord.standard_name
        coord.long_name = coord.standard_name
       
    return cube  

def get_pliomip1_data(fieldreq):
    """
    we will get the data from pliomip1
    returns an array of the mean data and the min, max mean, of the seasonal cycle
    """
    
    PLIOMIP1_FILE = (FILESTART + 'PLIOMIP1/means_for_' 
                     + fieldreq + '.txt')
    
    f1 = open(PLIOMIP1_FILE)
    
    lines = f1.readlines()
    lines[:] = [line.rstrip('\n') for line in lines]
            
    # find line which contains 'modelnameglobal[mean_ocean_eoi400'
    # which is the start of the land sea contrast
    string = 'modelnameglobal[mean_ocean_eoi400'
    land_amplification = []
    for i, line in enumerate(lines):
        if string in line:
            index = i
            
    for i in range(index + 1, len(lines)):
        line = lines[i]
        if 'modelname' in line:
            break   
        (model, mean_ocean_eoi400, meanocean_e280, meanocean_anom, 
        mean_land_eoi400, mean_land_e280, mean_land_anom) = line.split(',')
        amp = np.float(mean_land_anom) / np.float(meanocean_anom)
        if model != 'MEAN':
            land_amplification.append(amp)
            
    land_amp_arr = np.asarray(land_amplification, dtype=float)
    
    # find line which contains 'modelname20N-20S[mean_ocean_eoi400'
    # which is the start of the land sea contrast over the tropics
    string = 'modelname20N-20S[mean_ocean_eoi400'
    index=0
    land_amplification_20 = []
    for i, line in enumerate(lines):
        if string in line:
            index = i
     
            
    for i in range(index + 1, len(lines)):
        line = lines[i]
        if 'modelname' in line:
            break   
        (model, mean_ocean_eoi400, meanocean_e280, meanocean_anom, 
        mean_land_eoi400, mean_land_e280, mean_land_anom) = line.split(',')
        amp = np.float(mean_land_anom) / np.float(meanocean_anom)
        if model != 'MEAN':
            land_amplification_20.append(amp)
            
    land_amp20_arr = np.asarray(land_amplification_20, dtype=float)
    
    # find line which contains '45N-90N_anom'
    # which is the start of the fields averaged over certain regions
    string = '45N-90N_anom'
    index=0
    NH_SH_ratio45 = []
    PA_NH_60 = []
    PA_SH_60 = []
    for i, line in enumerate(lines):
        if string in line:
            index = i
     
            
    for i in range(index + 1, len(lines)):
        line = lines[i]
        if 'modelname' in line:
            break   
        (model, anom_45_90N, anom_45_90S, 
         anom_60_90N, anom_60_90S) = line.split(',')
        NH_SH_ratio = np.float(anom_45_90N) / np.float(anom_45_90S)
      
        if model != 'MEAN':
            NH_SH_ratio45.append(NH_SH_ratio)
            PA_NH_60.append(anom_60_90N)
            PA_SH_60.append(anom_60_90S)
   
    NH_SH_ratio45_arr = np.asarray(NH_SH_ratio45, dtype=float)
    PA_NH_60_arr = np.asarray(PA_NH_60, dtype = float)
    PA_SH_60_arr = np.asarray(PA_SH_60, dtype = float)
    
        
    return (land_amp_arr, land_amp20_arr, NH_SH_ratio45_arr, 
            PA_NH_60_arr, PA_SH_60_arr)

#####################################
def plotmean(modelnames,field,exptname,cntlname,linux_win,units):
   
    if linux_win=='w':
        fileoutstart='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\allplots\\'+field+'\\'
        exptlsm=FILESTART+'regridded/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
        cntllsm=FILESTART+'regridded/PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc'
    else:
        fileoutstart='/nfs/hera1/earjcti/regridded/allplots/'+field+'/'
        dataoutstart='/nfs/hera1/earjcti/regridded/alldata/'
        exptlsm=FILESTART+'PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
        cntllsm=FILESTART+'PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc'
     
 
    ########################################################
    # setup: get the lsm for the land sea contrast plot
    
    tempcube=iris.load_cube(exptlsm)
    cubegrid=iris.load_cube('one_lev_one_deg.nc')
    exptlsmcube=tempcube.regrid(cubegrid,iris.analysis.Linear())
    exptsea=(exptlsmcube.data - 1.0)*(-1.0)
    
   
    tempcube=iris.load_cube(cntllsm)
    cntllsmcube=tempcube.regrid(cubegrid,iris.analysis.Linear())
    cntlsea=(cntllsmcube.data - 1.0)*(-1.0)
    
 
    #########################################################
    # need to get data from annual mean plot
    
    nh_anomaly=np.zeros(len(modelnames))
    sh_anomaly=np.zeros(len(modelnames))
    nh_anomaly_extratropics=np.zeros(len(modelnames))
    sh_anomaly_extratropics=np.zeros(len(modelnames))
    nh_anomaly_polar=np.zeros(len(modelnames))
    sh_anomaly_polar=np.zeros(len(modelnames))
    all_anomaly=np.zeros(len(modelnames))
    land_anomaly=np.zeros(len(modelnames))
    sea_anomaly=np.zeros(len(modelnames))
    land_anomaly_tropics=np.zeros(len(modelnames))
    sea_anomaly_tropics=np.zeros(len(modelnames))
    
    for modelno in range(0,len(modelnames)):
        modeluse=modelnames[modelno]
        exptfile=FILESTART+'regridded/'+modeluse+'/'+exptname+'.'+field+'.allmean.nc'
        cntlfile=FILESTART+'regridded/'+modeluse+'/'+cntlname+'.'+field+'.allmean.nc'
        cube=iris.load_cube(exptfile)
        exptcube = resort_coords(cube)
        cube=iris.load_cube(cntlfile)
        cntlcube = resort_coords(cube)
    
        #######################################
        # compare NH vs SH temperature
        exptcube.coord('latitude').guess_bounds()
        exptcube.coord('longitude').guess_bounds()
        cntlcube.coord('latitude').guess_bounds()
        cntlcube.coord('longitude').guess_bounds()
        grid_areas_expt = iris.analysis.cartography.area_weights(exptcube)
        grid_areas_cntl = iris.analysis.cartography.area_weights(cntlcube)
        
        # exptcube
        nlat=len(exptcube.coord('latitude').points)
        lats=exptcube.coord('latitude').points
        
        grid_areas_nh=np.zeros(grid_areas_expt.shape)
        grid_areas_sh=np.zeros(grid_areas_expt.shape)
        grid_areas_nh_extratropics=np.zeros(grid_areas_expt.shape)
        grid_areas_sh_extratropics=np.zeros(grid_areas_expt.shape)
        grid_areas_nh_polar=np.zeros(grid_areas_expt.shape)
        grid_areas_sh_polar=np.zeros(grid_areas_expt.shape)
        polarval=60.0
        for j in range(0,nlat):
            if lats[j] <0:
                grid_areas_sh[j,:]=grid_areas_expt[j,:]
            else:
                grid_areas_nh[j,:]=grid_areas_expt[j,:]
            if lats[j] < -45.0:
                grid_areas_sh_extratropics[j,:]=grid_areas_expt[j,:]
            if lats[j] > 45.0:
                grid_areas_nh_extratropics[j,:]=grid_areas_expt[j,:]
            if lats[j] <= -1.0*polarval:
                grid_areas_sh_polar[j,:]=grid_areas_expt[j,:]
            if lats[j] >= polarval:
                grid_areas_nh_polar[j,:]=grid_areas_expt[j,:]
           
        expt_nh = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_nh)
        expt_sh = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sh)
        expt_anom = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN,weights=grid_areas_expt)
        expt_nh_extratropics = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_nh_extratropics)
        expt_sh_extratropics = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sh_extratropics)
        expt_nh_polar = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_nh_polar)
        expt_sh_polar = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sh_polar)
        
        
        
        #cntlcube
        nlat=len(cntlcube.coord('latitude').points)
        lats=cntlcube.coord('latitude').points
        grid_areas_nh=np.zeros(grid_areas_cntl.shape)
        grid_areas_sh=np.zeros(grid_areas_cntl.shape)
        for j in range(0,nlat):
            if lats[j] <0:
                grid_areas_sh[j,:]=grid_areas_cntl[j,:]
            else:
                grid_areas_nh[j,:]=grid_areas_cntl[j,:]
           
        cntl_nh = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_nh)
        cntl_sh = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sh)
        cntl_anom = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN,weights=grid_areas_cntl)
        cntl_nh_extratropics = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_nh_extratropics)
        cntl_sh_extratropics = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sh_extratropics)
        cntl_nh_polar = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_nh_polar)
        cntl_sh_polar = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sh_polar)
        
        
        nh_anomaly[modelno]=expt_nh.data-cntl_nh.data
        sh_anomaly[modelno]=expt_sh.data-cntl_sh.data
        nh_anomaly_extratropics[modelno]=expt_nh_extratropics.data-cntl_nh_extratropics.data
        sh_anomaly_extratropics[modelno]=expt_sh_extratropics.data-cntl_sh_extratropics.data
        nh_anomaly_polar[modelno]=expt_nh_polar.data-cntl_nh_polar.data
        sh_anomaly_polar[modelno]=expt_sh_polar.data-cntl_sh_polar.data
        all_anomaly[modelno]=expt_anom.data-cntl_anom.data
      
        
        #######################################
        # compare land with sea (globally and for tropics)
        
        
        # expt
        # first check grid
        for i in range(0,len(exptcube.coord('latitude').points)):
            if exptcube.coord('latitude').points[i] !=exptlsmcube.coord('latitude').points[i]:
                print('differences in lsm and gridded data',i)
                sys.exit(0)
       
    
        for i in range(0,len(exptcube.coord('longitude').points)):    
            if exptcube.coord('longitude').points[i] != exptlsmcube.coord('longitude').points[i]:
                print('differences in lsm and gridded data',i,exptcube.coord('longitude').points[i],
                  exptlsmcube.coord('longitude').points[i])
                sys.exit(0)
            
        grid_areas_land=grid_areas_expt * exptlsmcube.data  
        grid_areas_sea=grid_areas_expt * exptsea
        
        expt_land = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_land)
        expt_sea = exptcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sea)
        
        # get grid areas and experiment means for tropics
        nlat=len(exptcube.coord('latitude').points)
        lats=exptcube.coord('latitude').points
        grid_areas_land_tropics=np.zeros(grid_areas_land.shape)
        grid_areas_sea_tropics=np.zeros(grid_areas_sea.shape)
        for j in range(0,nlat):
            if ((lats[j] < 20.) and (lats[j] > -20.)):
                grid_areas_land_tropics[j,:]=grid_areas_land[j,:]
                grid_areas_sea_tropics[j,:]=grid_areas_sea[j,:]
    
        expt_land_tropics = exptcube.collapsed(['longitude', 'latitude'],
                                               iris.analysis.MEAN, weights=grid_areas_land_tropics)
        expt_sea_tropics = exptcube.collapsed(['longitude', 'latitude'], 
                                              iris.analysis.MEAN, weights=grid_areas_sea_tropics)
    
        # cntl
        # first check grid
        for i in range(0,len(cntlcube.coord('latitude').points)):
            if cntlcube.coord('latitude').points[i] !=cntllsmcube.coord('latitude').points[i]:
                print('differences in lsm and gridded data',i)
                sys.exit(0)
       
    
        for i in range(0,len(cntlcube.coord('longitude').points)):    
            if cntlcube.coord('longitude').points[i] != cntllsmcube.coord('longitude').points[i]:
                print('differences in lsm and gridded data',i,cntlcube.coord('longitude').points[i],
                  cntllsmcube.coord('longitude').points[i])
                sys.exit(0)
            
        grid_areas_land=grid_areas_cntl * cntllsmcube.data  
        grid_areas_sea=grid_areas_cntl * cntlsea
    
        cntl_land = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_land)
        cntl_sea = cntlcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_sea)
    
        # get grid areas and experiment means for tropics
        nlat=len(cntlcube.coord('latitude').points)
        lats=cntlcube.coord('latitude').points
        grid_areas_land_tropics=np.zeros(grid_areas_land.shape)
        grid_areas_sea_tropics=np.zeros(grid_areas_sea.shape)
        for j in range(0,nlat):
            if ((lats[j] < 20.) and (lats[j] > -20.)):
                grid_areas_land_tropics[j,:]=grid_areas_land[j,:]
                grid_areas_sea_tropics[j,:]=grid_areas_sea[j,:]
    
        cntl_land_tropics = cntlcube.collapsed(['longitude', 'latitude'],
                                               iris.analysis.MEAN, weights=grid_areas_land_tropics)
        cntl_sea_tropics = cntlcube.collapsed(['longitude', 'latitude'], 
                                              iris.analysis.MEAN, weights=grid_areas_sea_tropics)
    
    
        land_anomaly[modelno]=expt_land.data-cntl_land.data
        sea_anomaly[modelno]=expt_sea.data-cntl_sea.data
    
        land_anomaly_tropics[modelno]=expt_land_tropics.data-cntl_land_tropics.data
        sea_anomaly_tropics[modelno]=expt_sea_tropics.data-cntl_sea_tropics.data
    
    # get data from PlioMIP1 if required
    if PLIOMIP1 == 'y':
        (pliomip1_landsea_amp, 
         pliomip1_tropics_landsea_amp,
         pliomip1_NH_SH_ratio45,
         pliomip1_PA_NH,
         pliomip1_PA_SH) = get_pliomip1_data(field)
   
    
    # plot NH SH contrast
    ax=plt.subplot(1,1,1)
    ax.plot(modelnames,nh_anomaly,'x',label='NH anomaly')
    ax.plot(modelnames,sh_anomaly,'x',label='SH anomaly')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+(0.3*box.height), box.width * 0.8, (0.8*box.height)])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel(units)
    plt.xticks(rotation='90')
    plt.title('mPWP - PI anomaly for each hemisphere')
    fileout=fileoutstart+'/hemisphere_difference.eps'
    plt.savefig(fileout)
    fileout=fileoutstart+'/hemisphere_difference.pdf'
    plt.savefig(fileout)
    plt.close()
    
    
    # plot NH SH contrast for extratropics
    ax=plt.subplot(2,1,1)
    ax.plot(modelnames,nh_anomaly_extratropics,'x',label='>45N anom')
    ax.plot(modelnames,sh_anomaly_extratropics,'x',label='<45S anom')
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+(0.3*box.height), box.width * 0.8, (0.8*box.height)])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel(units)
    #plt.xticks(rotation='45')
    ax.tick_params(axis='x',labelbottom='False')
    plt.title('Plio_core - PI_Ctl; extratropical NH/SH anomaly')
    plt.figtext(0.02, 0.97,'c)',
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    
    ax=plt.subplot(2,1,2)
    ax.plot(modelnames,nh_anomaly_extratropics/sh_anomaly_extratropics
            ,'x',label='>45N/<45S')
    ax.plot(modelnames,(np.zeros(len(modelnames))+1.0))
    
    if PLIOMIP1 == 'y':
       for mod_p1 in pliomip1_NH_SH_ratio45:
            ax.axhline(y=mod_p1, color='grey', alpha=0.4)
    
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+(0.4*box.height), box.width * 0.8, (0.9*box.height)])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('factor')
    plt.xticks(rotation='90')
    #plt.title('mPWP - PI anomaly')
    
    fileout=fileoutstart+'/hemisphere_difference_extratropics.eps'
    plt.savefig(fileout)
    fileout=fileoutstart+'/hemisphere_difference_extratropics.pdf'
    plt.savefig(fileout)
    plt.close()
    
    # write out hemisphere difference extratropics
    
    txtout = open(FILEOUT, "w+") 
    if field == 'NearSurfaceTemperature':
        txtout.write('data for 3c \n')
        
        writedata = ("model_name, nh_anom_et, sh_anom_et \n")
    
        txtout.write(writedata)
        for i, mod in enumerate(modelnames):
            writedata = (mod + ',' + 
                         np.str(np.around(nh_anomaly_extratropics[i],2)) + ',' + 
                         np.str(np.around(sh_anomaly_extratropics[i],2)) + '\n')
            txtout.write(writedata)
   
     # plot polar amplification 
    ax=plt.subplot(1,1,1)
    labelname='>'+np.str(np.int(polarval))+'N amplification'
    ax.plot(modelnames,nh_anomaly_polar/all_anomaly,'x',label=labelname)
    labelname='<'+np.str(np.int(polarval))+'S amplification'
    ax.plot(modelnames,sh_anomaly_polar/all_anomaly,'x',label=labelname)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+(0.3*box.height), box.width * 0.7, (0.7*box.height)])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel('factor')
    plt.xticks(rotation='90')
    ax.axhline(y=1.0, xmin=0.0, xmax=len(modelnames), color='r')
    #ax.tick_params(axis='x',labelbottom='False')
    plt.title('Plio_Core - PI_Ctl; polar amplification factor')
    plt.figtext(0.02, 0.97,'d)',
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    
    print('NH polar amplification', np.mean(nh_anomaly_polar) / np.mean(all_anomaly))
    print('SH polar amplification', np.mean(sh_anomaly_polar) / np.mean(all_anomaly))
    print('total polar amplification', np.mean(nh_anomaly_polar + sh_anomaly_polar) / (2.0 * np.mean(all_anomaly)))
    print('all polar amplification', (nh_anomaly_polar + sh_anomaly_polar) / (2.0 * all_anomaly))
    nh_amps = [x1 / x2 for (x1, x2) in zip(nh_anomaly_polar, all_anomaly)] 
    sh_amps = [x1 / x2 for (x1, x2) in zip(sh_anomaly_polar, all_anomaly)] 
    print('nh amps unsorted', nh_amps)
    print('nh amps',np.sort(nh_amps))
    print('nh median', np.median(nh_amps))
    print('sh median', np.median(sh_amps))
    print('nh percentiles 10/50/90',np.percentile(nh_amps, 10),
          np.percentile(nh_amps,50), np.percentile(nh_amps,90))
    print('sh percentiles 10/50/90',np.percentile(sh_amps, 10),
          np.percentile(sh_amps,50), np.percentile(sh_amps,90))
   
   
    
    fileout=fileoutstart+'/polar_amplification_'+np.str(np.int(polarval))+'.eps'
    plt.savefig(fileout)
    fileout=fileoutstart+'/polar_amplification_'+np.str(np.int(polarval))+'.pdf'
    plt.savefig(fileout)
    plt.close()
    
    if field == 'NearSurfaceTemperature':
        txtout.write('data for 3d \n')
        
        writedata = ("model_name, nh_anom_polar, sh_anom_polar, global_anom \n")
    
        txtout.write(writedata)
        for i, mod in enumerate(modelnames):
            writedata = (mod + ',' + 
                         np.str(np.around(nh_anomaly_polar[i],2)) + ',' + 
                         np.str(np.around(sh_anomaly_polar[i],2)) + ',' + 
                         np.str(np.around(all_anomaly[i],2)) + 
                         '\n')
            txtout.write(writedata)
    
    
    # plot land sea contrast
    if field != 'TotalPrecipitation':
        ax=plt.subplot(2, 1, 1)
    else:
        ax = plt.subplot(1, 1, 1)
    ax.plot(modelnames,land_anomaly,'x',label='Land anomaly')
    ax.plot(modelnames,sea_anomaly,'x',label='Sea anomaly')
    box = ax.get_position()
    if field != 'TotalPrecipitation':
        ax.set_position([box.x0, box.y0+(0.3*box.height), box.width * 0.8, (0.8*box.height)])
        ax.tick_params(axis='x',labelbottom='False')
        plt.figtext(0.02, 0.97,'b)',
                       horizontalalignment='left',
                       verticalalignment='top',
                       fontsize=20)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    
    plt.ylabel(units)
   
    #plt.xticks(rotation='45')
    plt.title('Plio_Core - PI_Ctrl global land-sea anomaly')
    
    
    # plot land sea contrast tropics]
    
    ax2=plt.subplot(2,1,2)
    ax2.plot(modelnames,land_anomaly_tropics,'x',label='Land anomaly')
    ax2.plot(modelnames,sea_anomaly_tropics,'x',label='Sea anomaly')
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0+(0.4*box.height), box.width * 0.8, (0.9*box.height)])
    ax2.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.ylabel(units)
    plt.xticks(rotation='90')
    plt.title('Plio_Core - PI_Ctrl 20N-20S land-sea anomaly')
    
    print('sea anomaly tropics',np.mean(sea_anomaly_tropics))

    fileout=fileoutstart+'/land_sea_contrast.pdf'
    plt.savefig(fileout)
    fileout=fileoutstart+'/land_sea_contrast.eps'
    plt.savefig(fileout)
    plt.close()
    
    # write out
    
    if field == 'NearSurfaceTemperature':
        txtout.write('data for 3b \n')
    else:
        txtout.write('data for 6b \n' )
        
    writedata = ("model_name, land_anomaly, sea_anomaly, " + 
                 "tropical_land_anomaly, tropical_sea_anomaly\n")
    
    
    txtout.write(writedata)
    for i, mod in enumerate(modelnames):
            writedata = (mod + ',' + 
                         np.str(np.around(land_anomaly[i],2)) + ',' + 
                         np.str(np.around(sea_anomaly[i],2)) + ',' + 
                         np.str(np.around(land_anomaly_tropics[i],2)) + ',' + 
                         np.str(np.around(sea_anomaly_tropics[i],2)) + '\n')
            print(writedata)
            txtout.write(writedata)
        
    txtout.close
    print('mean land anomaly', np.mean(land_anomaly))
    print('mean sea anomaly', np.mean(sea_anomaly))
    
    ######################################
    # plot land sea contrast as a factor
    factor_land=land_anomaly / sea_anomaly
    factor_land_tropics=land_anomaly_tropics / sea_anomaly_tropics
   
    if field != 'NearSurfaceTemperature':
        ax=plt.subplot(2,1,1)
    else:
        ax = plt.subplot(1, 1, 1)
    ax.plot(modelnames,factor_land,'x')
    ax.plot(modelnames,(np.zeros(len(modelnames))+1.0))
    if PLIOMIP1 == 'y':
        for mod_amp in pliomip1_landsea_amp:
            ax.axhline(y=mod_amp, color='grey', alpha=0.4)
           
    box = ax.get_position()
    if field != 'NearSurfaceTemperature':
        ax.set_position([box.x0 + (0.1 * box.width),
                         box.y0+(0.3*box.height), 
                         box.width * 0.8, (0.8*box.height)])
        ax.tick_params(axis='x',labelbottom='False')
    else:
        ax.set_position([box.x0, box.y0+(0.3*box.height), 
                         box.width, (0.7*box.height)])
    plt.figtext(0.02, 0.97,'b)',
                horizontalalignment='left',
                verticalalignment='top',
                fontsize=20)
    plt.xticks(rotation='90')
    plt.ylabel('land amplification')
   
    #plt.xticks(rotation='45')
    plt.title('mPWP - PI; land_anomaly / sea anomaly')
    
    if field != 'NearSurfaceTemperature':
        ax=plt.subplot(2,1,2)
        ax.plot(modelnames,factor_land_tropics,'x')
        ax.plot(modelnames,(np.zeros(len(modelnames))+1.0))
        if PLIOMIP1 == 'y':
            for mod_amp in pliomip1_tropics_landsea_amp:
                ax.axhline(y=mod_amp, color='grey', alpha=0.4)
        
        box = ax.get_position()
        ax.set_position([box.x0 + (0.1 * box.width), 
                         box.y0+(0.3*box.height), 
                        box.width * 0.8, (0.9*box.height)])
        plt.ylabel('land amplification')
        plt.xticks(rotation='90')
        plt.title('mPWP - PI (20N-20S); land_anomaly / sea anomaly')
        plt.figtext(0.02, 0.97,'b)',
                   horizontalalignment='left',
                   verticalalignment='top',
                   fontsize=20)
    fileout=fileoutstart+'/land_sea_amplification.eps'
    plt.savefig(fileout)
    fileout=fileoutstart+'/land_sea_amplification.pdf'
    plt.savefig(fileout)
    plt.close()
   
    
   
   
    
   
    
    
    
        
################################################################
def plotmap(modelnames,field,exptname,cntlname,linux_win,units):  
    # this subprogram will 
    # 1. read in the data from each map.  
    # 2. It will normalise the data by subtracting the mean and dividing by the spatial standard deviation
    # 3. it will then plot the field change by number of standard deviations above and below the m2an

    normalized_cubes=iris.cube.CubeList([])
    standard_dev=[]
    if linux_win=='w':
        FILESTART='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
        fileoutstart='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\allplots\\'+field+'\\'
        exptlsm=FILESTART+'regridded/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    else:
        FILESTART='/nfs/hera1/earjcti/'
        fileoutstart='/nfs/hera1/earjcti/regridded/allplots/'+field+'/'
        dataoutstart='/nfs/hera1/earjcti/regridded/alldata/'
        exptlsm=FILESTART+'PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
      
    
    lsmcube=iris.load_cube(exptlsm)
        
    for modelno in range(0,len(modelnames)):
        modeluse=modelnames[modelno]
        print(modeluse)
        exptfile=FILESTART+'regridded/'+modeluse+'/'+exptname+'.'+field+'.allmean.nc'
        cntlfile=FILESTART+'regridded/'+modeluse+'/'+cntlname+'.'+field+'.allmean.nc'
        cube=iris.load_cube(exptfile)
        exptcube = resort_coords(cube)
        exptcube.data=exptcube.data.astype('float32') # change to float32 for concatentation later
       
       
        for coord in exptcube.coords():
            if coord.standard_name !='longitude' and coord.standard_name !='latitude':
                exptcube.remove_coord(coord)
                
        cube=iris.load_cube(cntlfile)
        cntlcube = resort_coords(cube)
        cntlcube.data=cntlcube.data.astype('float32')
        
        for coord in cntlcube.coords():
            if coord.standard_name !='longitude' and coord.standard_name !='latitude':
                cntlcube.remove_coord(coord)
    
        
        diffcube=exptcube-cntlcube
       
        #print(exptcube.coord('surface'))
      
        #######################################
        # get mean and spatial standard deviation
        diffcube.coord('latitude').guess_bounds()
        diffcube.coord('longitude').guess_bounds()
       
        grid_areas_diff = iris.analysis.cartography.area_weights(diffcube)
       
        diffmean = diffcube.collapsed(['longitude', 'latitude'], iris.analysis.MEAN, weights=grid_areas_diff)
        diffsd = diffcube.collapsed(['longitude', 'latitude'], iris.analysis.STD_DEV)
        
        standard_dev.append(diffsd.data)
       
        ########################################
        # normalise the cube by dividing by the mean and std dev
        # and plot
        
        normcube=(diffcube-diffmean)/diffsd
        print(normcube)
      
        
        V=np.arange(-4,4.5,0.5)
        mycmap = plt.cm.get_cmap('RdBu_r',len(V+2))
        newcolors=mycmap(np.linspace(0,1,len(V+2)))
        white=([1,1,1,1])
        print((len(V)/2)-1,(len(V)/2)+2)
        print((np.ceil(len(V)/2)-1),np.floor((len(V)/2)+2))
        
        newcolors[np.int(np.ceil((len(V)/2)-1)):np.int(np.floor((len(V)/2)+2))
                  ,:]=white
        
        mycmap=ListedColormap(newcolors)
           
        #Draw the contour w
       
        qplt.contourf(normcube, V,extend='both',cmap=mycmap)
        #qplt.contourf(normcube, V,extend='both',cmap='RdBu_r')
        
       
        #sys.exit()
        qplt.contour(lsmcube,1,colors='black')
        plt.title(modeluse+': '+field+'\n No of stddev from mean'
                  +'('+(np.str(np.round(diffsd.data,1)))+units+')')

        if linux_win=='l':
            fileout=fileoutstart+'Map_normalised/'+modeluse+'.eps'
            plt.savefig(fileout)
            
        fileout=fileoutstart+'Map_normalised//'+modeluse+'.pdf'
        plt.savefig(fileout)
        plt.close()
        
        #normcube.Coord(1, standard_name='model', long_name='model', var_name='model', units='1', 
       #                       bounds=None, attributes=None, coord_system=None)
        tempcube=iris.util.new_axis(normcube)
        tempcube.add_dim_coord(iris.coords.DimCoord(modelno, 
                standard_name='model_level_number', long_name='model', 
                var_name='model', 
                units=None,
                bounds=None,
                coord_system=None, circular=False),0)
        # tempcube needs to be dtype=float64
        
        
      
        normalized_cubes.append(tempcube)
        
    
    ################################################
    # END OF MODEL LOOP
    equalise_attributes(normalized_cubes)
   
    print(normalized_cubes)
    print(normalized_cubes[0])
    print(normalized_cubes[1])
    print(normalized_cubes[0].dtype)
    print(normalized_cubes[1].dtype)
    allnormcube=normalized_cubes.concatenate_cube()
    
    meancube=allnormcube.collapsed(['model_level_number'], iris.analysis.MEAN)
    maxcube=allnormcube.collapsed(['model_level_number'], iris.analysis.MAX)
    mincube=allnormcube.collapsed(['model_level_number'], iris.analysis.MIN)
    mediancube=allnormcube.collapsed(['model_level_number'], iris.analysis.MEDIAN)
    
    ###########################
    # plot the mean value
    
    V=np.arange(-2.5,2.75,0.25)
    mycmap = plt.cm.get_cmap('RdBu_r',len(V+2))
    newcolors=mycmap(np.linspace(0,1,len(V+2)))
    white=([1,1,1,1])
    newcolors[np.int(np.ceil((len(V)/2)-1)):np.int(np.floor((len(V)/2)+2))
                  ,:]=white
    mycmap=ListedColormap(newcolors)
        
    qplt.contourf(meancube, V,extend='both',cmap=mycmap)
    qplt.contour(lsmcube,1,colors='black') 
    minval=np.str(np.round(np.amin(standard_dev),1))
    maxval=np.str(np.round(np.amax(standard_dev),1))
    print(standard_dev)
    print(np.amin(standard_dev))
    print(np.amax(standard_dev))
    print(minval,maxval)
    #sys.exit(0)
    plt.title(field+'\n Mean No of stddev from mean'
                  +'('+minval+'-'+maxval+units+')')

    if linux_win=='l':
            fileout=fileoutstart+'Map_normalised/meandiff.eps'
            plt.savefig(fileout)
            
    fileout=fileoutstart+'Map_normalised//meandiff.pdf'
    plt.savefig(fileout)
    plt.close()
    plt.show()
        
    
    ###########################
    # plot the median value
    
        
    qplt.contourf(mediancube, V,extend='both',cmap=mycmap)
    qplt.contour(lsmcube,1,colors='black') 
    minval=np.str(np.round(np.amin(standard_dev),1))
    maxval=np.str(np.round(np.amax(standard_dev),1))
    plt.title(field+'\n Median No of stddev from mean'
                  +'('+minval+'-'+maxval+units+')')

    if linux_win=='l':
            fileout=fileoutstart+'Map_normalised/mediandiff.eps'
            plt.savefig(fileout)
            
    fileout=fileoutstart+'Map_normalised//mediandiff.pdf'
    plt.savefig(fileout)
    plt.close()
    plt.show()
        
    
    ###########################
    # plot the maximum value
    
    V=np.arange(-4.5,4.75,0.25)
    mycmap = plt.cm.get_cmap('RdBu_r',len(V+2))
    newcolors=mycmap(np.linspace(0,1,len(V+2)))
    white=([1,1,1,1])
    newcolors[np.int(np.ceil((len(V)/2)-1)):np.int(np.floor((len(V)/2)+2))
                  ,:]=white
    mycmap=ListedColormap(newcolors)
        
    qplt.contourf(maxcube, V,extend='both',cmap=mycmap)
    qplt.contour(lsmcube,1,colors='black') 
    minval=np.str(np.round(np.amin(standard_dev),1))
    maxval=np.str(np.round(np.amax(standard_dev),1))
    plt.title(field+'\n Maximum No of stddev from mean'
                  +'('+minval+'-'+maxval+units+')')

    if linux_win=='l':
            fileout=fileoutstart+'Map_normalised/maxdiff.eps'
            plt.savefig(fileout)
            
    fileout=fileoutstart+'Map_normalised//maxdiff.pdf'
    plt.savefig(fileout)
    plt.close()
   
    
    ###########################
    # plot the minimum value
    
    V=np.arange(-2.5,2.75,0.25)
    mycmap = plt.cm.get_cmap('RdBu_r',len(V+2))
    newcolors=mycmap(np.linspace(0,1,len(V+2)))
    white=([1,1,1,1])
    newcolors[np.int(np.ceil((len(V)/2)-1)):np.int(np.floor((len(V)/2)+2))
                  ,:]=white
    mycmap=ListedColormap(newcolors)
        
    qplt.contourf(mincube, V,extend='both',cmap=mycmap)
    qplt.contour(lsmcube,1,colors='black') 
    minval=np.str(np.round(np.amin(standard_dev),1))
    maxval=np.str(np.round(np.amax(standard_dev),1))
    plt.title(field+'\n Minimum No of stddev from mean'
                  +'('+minval+'-'+maxval+units+')')

    if linux_win=='l':
            fileout=fileoutstart+'Map_normalised/mindiff.eps'
            plt.savefig(fileout)
            
    fileout=fileoutstart+'Map_normalised//mindiff.pdf'
    plt.savefig(fileout)
    plt.close()
   

        
         


##########################################################
# main program

filename=' '
linux_win='l'
if linux_win=='w':
    FILESTART='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
else:
    FILESTART='/nfs/hera1/earjcti/'
        

modelnames=['CESM2', 'IPSLCM6A', 'COSMOS', 
            'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
            'MIROC4m', 'IPSLCM5A2', 'HadCM3',
            'GISS2.1G', 'CCSM4', 
            'CCSM4-Utr', 'CCSM4-UoT', 
            'NorESM-L', 'MRI2.3', 'NorESM1-F'
            ]


#modelnames=['HadCM3','NorESM-L']
fieldnames=['NearSurfaceTemperature']
#units=['degC']
#
#fieldnames=['TotalPrecipitation']
units=['mm/day']
exptname='EOI400'
cntlname='E280'
PLIOMIP1 = 'n'

for field in range(0,len(fieldnames)):
    if fieldnames[field] == 'TotalPrecipitation':
        FILEOUT = FILESTART + 'regridded/alldata/data_for_6b.txt'
        
    if fieldnames[field] == 'NearSurfaceTemperature':
        FILEOUT = FILESTART + 'regridded/alldata/data_for_3b_3c_3d.txt'
       
    
    # will plot NH vs SHPlio_enh_LSM_v1.0Plio_enh_LSM_v1.0
    plotmean(modelnames,fieldnames[field],exptname,cntlname,linux_win,units[field])
     
    # plotmap will show where the field is larger or smaller than the mean  
    plotmap(modelnames,fieldnames[field],exptname,cntlname,linux_win,units[field])

#sys.exit(0)
