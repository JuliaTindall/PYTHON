#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on April 29th


#@author: earjcti
#
# This program will read all the means from the regridded files and plot them


import os
import numpy as np
import scipy as sp
#import cf
import iris
import iris.util
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
#from mpl_toolkits.basemap import Basemap, shiftgrid
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
import cf_units as unit
from iris.experimental.equalise_cubes import equalise_attributes
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

import sys


def plotmean_newaxis(cube, modelno_):
     tempcube=iris.util.new_axis(cube)
     tempcube.add_dim_coord(iris.coords.DimCoord(modelno_, 
            standard_name='model_level_number', long_name='model', 
            var_name='model', 
            units=None,
            bounds=None,
            coord_system=None, circular=False),0) 
     return tempcube
 

def resort_coords(cube,levelno):
    """
    this will make all the dimensions of the cube match.  They will all be
    longitude, latitude, level-no (ie 1 for first model, 2 for second model...)
    
    input is the cube and the level number
    output is the cube with the new dimensions
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
         
    newcube = plotmean_newaxis(cube, [levelno])
    # this will make sure cell_methods match and that cubes can
    # be concatenated
    newcube.cell_methods = None
    
        
    return newcube
###########################################
def get_NH_mean(modelname, expt, field):
     """
     gets the mean of the NH for model - modelname
                             and expt - experiment name (ie PI)
     returns a numpy array of length 12 with the average for each month
     """
     # read data into iris cube
     filename = (FILESTART + modelname + 
                  '/' + expt + '.' + field + 
                  '.mean_month.nc')
     
     if (modelname == 'GISS2.1G' or modelname == 'IPSLCM6A'
         or (modelname == 'IPSLCM5A' and field == 'SST')
         or (modelname == 'IPSLCM5A2' and field == 'SST')
         or (modelname == 'NorESM-L' and field == 'SST')
         or (modelname == 'NorESM1-F' and field == 'SST')):
          cubes = iris.load(filename)
          cube = cubes[0]
     else:
          cube = iris.load_cube(filename, field)
        
     # get weights and average over NH
     cube.coord('latitude').guess_bounds()
     cube.coord('longitude').guess_bounds()
     grid_areas = iris.analysis.cartography.area_weights(cube)
     grid_areas_nh = np.zeros(grid_areas.shape)
     for j, lat in enumerate(cube.coord('latitude').points):
          if lat > 0:
               grid_areas_nh[:, j, :] = grid_areas[:, j, :]

     cube_nh = (cube.collapsed(['longitude', 'latitude'],
                iris.analysis.MEAN, weights = grid_areas_nh))
     
     return cube_nh.data
 
def get_pliomip1_data(fieldreq):
    """
    we will get the data from pliomip1
    returns an array of the mean data and the min, max mean, of the seasonal cycle
    """
    
    PLIOMIP1_FILE = (FILESTART[:-10] + 'PLIOMIP1/means_for_' 
                     + fieldreq + '.txt')
    
    f1 = open(PLIOMIP1_FILE)
    
    lines = f1.readlines()
    lines[:] = [line.rstrip('\n') for line in lines]
    
    # means are at the start of the file
    allanoms_list = []
    for i in range(1, len(lines)):
        line = lines[i]
        if line[0:9] == 'modelname':
            break    # we have now got all the means
        modname, eoi400, e280, anom = line.split(',')
        if modname == 'MEAN':
            meananom=anom
        else:
            allanoms_list.append(anom)
            
    # find line which contains 'jan feb mar ' which is the start of the seasonal cycle
    string = 'jan feb mar'
    min_seas_cyc = np.zeros(12) + 1000.
    max_seas_cyc = np.zeros(12) - 1000.
    for i, line in enumerate(lines):
        if string in line:
            index = i
            
    for i in range(index + 1, len(lines)):
        line = lines[i]
        if 'modelname' in line:
            break   
        modname, anomstr = line.split(',')
        anom_arr = np.array(anomstr.strip('[]').split(), dtype=float)
        if modname == 'MEAN':
            mean_seas_cyc = anom_arr
        else:
            for j in range(0,12):
                min_seas_cyc[j] = np.min([anom_arr[j], min_seas_cyc[j]]) 
                max_seas_cyc[j] = np.max([anom_arr[j], max_seas_cyc[j]]) 
    

    allanoms = np.asarray(allanoms_list, dtype=float)
    
        
    return allanoms, min_seas_cyc, max_seas_cyc, mean_seas_cyc
  
   

#####################################
def plotmean(modelnames,field,exptname,cntlname,linux_win,units,individual_plot):
   
    
    
    names = {"EOI400" : "Plio_Core",
             "E280" : "PI_Ctrl"
            }
    namefield = {"NearSurfaceTemperature" : "SAT",
                 "TotalPrecipitation" : "Precipitation",
                 "SST" : "SST"
                 }
    
    # set up lists to store all values for both experiment and control    
    model_global_mean=[]
    model_global_sd=[]
    
    # store months in a numpy array (model, experiment,month)
    monmeans=np.zeros((len(modelnames),2,12))   
    monsd=np.zeros((len(modelnames),2,12))
    
    # store latitudes in a numpy array (model, experiment,latitude)
    latmeans=np.zeros((len(modelnames),2,180))   
    latsd=np.zeros((len(modelnames),2,180))
    lats=np.zeros(180)
     
    for modelno in range(0,len(modelnames)):
        modeluse=modelnames[modelno]
        filenames=[]
        filenames.append(FILESTART+modeluse+'/'+exptname+'.'+field+'.data.txt')
        filenames.append(FILESTART+modeluse+'/'+cntlname+'.'+field+'.data.txt')
       
        # set up temporary lists to store data from each experiment
        means=[]
        sds=[]
        nmon=12
       
        for fileno in range(0,len(filenames)):
            
            f=open(filenames[fileno],"r")
            f1=f.readlines()
            f2 = [x.replace('\n', '') for x in f1]
            
            # get the means according to their position in the file
            all_mean_sd=f2[2]
            all_mon_mean_sd=f2[5:5+12]
            all_lat_mean_sd=f2[20:20+180]
           
            # extract global mean
            mean,sd=all_mean_sd.split(',')
            means.append(mean)
            sds.append(sd)
            
            # extract monthly means 
            for x in all_mon_mean_sd:
                mon,mean,sd=x.split(',')
                monmeans[modelno,fileno,int(mon)-1]=float(mean)
                monsd[modelno,fileno,int(mon)-1]=float(sd)
            
            
            # extract latitude means
            for x in all_lat_mean_sd:
                lat,mean,sd=x.split(',')
                latss=int(float(lat)+89.5) # convert latitude to a subscript
               
                if mean != ' --' and mean != '--':
                    latmeans[modelno,fileno,latss]=float(mean) # stores latitudinal means
                else:
                    latmeans[modelno,fileno,latss]=np.nan
                if sd != ' --' and sd != '--' :
                    latsd[modelno,fileno,latss]=float(sd)
                else:
                    latsd[modelno,fileno,latss]=np.nan
                lats[latss]=lat # stores latitudes
            
            
        model_global_mean.append(means)
        model_global_sd.append(sds)
           
    ############################################################
    # get the monthly means for the NH
    # 
    # store months in a numpy array (model, experiment,month)
    # ss 1 is experiment ss 2 is control

    monmeans_NH = np.zeros((len(modelnames),2,12))   
    monsd = np.zeros((len(modelnames),2,12))

    for i, model in enumerate(modelnames):
         expt_mon_mean = get_NH_mean(model, exptname, field)
         monmeans_NH[i, 0, :] = expt_mon_mean

         cntl_mon_mean = get_NH_mean(model, cntlname, field)
         monmeans_NH[i, 1, :] = cntl_mon_mean

    
    #===============================================================
    # if pliomip1 is set get pliomip1 data
    if PLIOMIP1 == 'y':
        (mean_pliomip1, min_seas_pliomip1,
         max_seas_pliomip1, mean_seas_pliomip1) = get_pliomip1_data(field)
    
    #############################################################
    # plot the global mean and error bars from each model.
    
   
    expt_global_mean=[float(item[0]) for item in model_global_mean]
    expt_global_2sigma=([float(item[0])*2.0 for item in model_global_sd])

    cntl_global_mean=[float(item[1]) for item in model_global_mean]
    cntl_global_2sigma=([float(item[1])*2.0 for item in model_global_sd])
  
    
   
   
    #fig,ax=plt.subplots(2,1,1)
    ax=plt.subplot(2,1,1)
   
    ax.errorbar(modelnames,expt_global_mean,
                yerr=expt_global_2sigma,fmt='x',label=names.get(exptname))
    ax.errorbar(modelnames,cntl_global_mean,
                yerr=cntl_global_2sigma,fmt='x',label=names.get(cntlname))
    # Shrink current axis by 20% and put a legend to the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+(0.15*box.height), box.width * 0.8, box.height])
    plt.figtext(0.02, 0.97,'a)',
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    
    titlename='Global mean ' + namefield.get(field)
    plt.title(titlename)
    ax.tick_params(axis='x',labelbottom='False')
    plt.ylabel(units)
    
    
    
    
    ax=plt.subplot(2,1,2)
    ax.plot(modelnames,[x1 - x2 for (x1, x2) in zip(expt_global_mean, cntl_global_mean)],'x')
    print('means',[x1 - x2 for (x1, x2) in zip(expt_global_mean, cntl_global_mean)])
    anomalies = [x1 - x2 for (x1, x2) in zip(expt_global_mean, cntl_global_mean)]
    print('multimodelmean = ',
          np.mean(expt_global_mean) - np.mean(cntl_global_mean), np.mean(anomalies))
    print('multimodelmedian = ', np.median(anomalies))
    print('percentiles 10/50/90',np.percentile(anomalies, 10),
          np.percentile(anomalies,50), np.percentile(anomalies,90))
    sys.exit(0)
    
    sorted_anomalies = np.sort(anomalies)
    print(sorted_anomalies)
    #print('means % change',[((x1 - x2) *100 / x2) for (x1, x2) in zip(expt_global_mean, cntl_global_mean)])
    #print('multimodelmean = ',(np.mean(expt_global_mean) - np.mean(cntl_global_mean)) *100. / np.mean(cntl_global_mean))
     # if pliomip1 is set overplot pliomip1 means as grey horizontal bars
    if (PLIOMIP1 == 'y' and field !='SST'):
        for mean_mod in mean_pliomip1:
            ax.axhline(y=mean_mod, color='grey', alpha=0.4)
    # Shrink axis as appropriate
    box = ax.get_position()
    ax.set_position([box.x0, box.y0+(0.3*box.height), box.width * 0.8, box.height*0.8])
    plt.title('Global mean ' + namefield.get(field) + ' anomaly')
    plt.ylabel(units)
    plt.xticks(rotation='90', fontsize=8)
    #plt.xticks(x, labels, rotation='vertical')
   

    fileout=FILESTART+'allplots/'+field+'/global_means.eps'
    plt.savefig(fileout)
    fileout=FILESTART+'allplots/'+field+'/global_means.pdf'
    plt.savefig(fileout)
    plt.close()
    
    # write out data
    
    txtfile1 = open(FILEOUT,"w+")
    if field == 'NearSurfaceTemperature':
        txtfile1.write("data_for_1a \n")
    elif field == 'TotalPrecipitation':
        txtfile1.write("data_for_5a \n")
        
    txtfile1.write('modelname, Pliocore_global_mean, pliocore_global_mean_2sigma, ' + 
                   'picntl_global_mean, picntl_global_mean_2sigma \n')
    for i, model in enumerate(modelnames):
        txtfile1.write((model + ',' + np.str(np.around(expt_global_mean[i],2)) + 
                       ',' + np.str(np.around(expt_global_2sigma[i],2)) + 
                       ',' + np.str(np.around(cntl_global_mean[i],2)) + 
                       ',' + np.str(np.around(cntl_global_2sigma[i],2)) + '\n'))
        
    txtfile1.write('\n')
    
      
    
    ################################################################
    #  plot the NH seasonal cycle from each model
   

    plt.subplot(2,1,1)
    for i in range(0,len(monmeans_NH[:,0,0])):
        # plot experiment data
       
        plt.plot(monmeans_NH[i,0,:],color='r')
        # plot control data
        plt.plot(monmeans_NH[i,1,:],color='b')
        plt.title('NH annual cycle of '+field)
        plt.ylabel(units)
        
    plt.subplot(2,1,2)
    plt.plot(np.mean(monmeans_NH[:,0,:],axis=0),label=exptname,color='r')
    plt.plot(np.mean(monmeans_NH[:,1,:],axis=0),label=cntlname,color='b')
    plt.ylabel(units)
    plt.xlabel('month')    
    plt.legend()
    fileout=FILESTART+'allplots/'+field+'/seasonal_cycle_all_models.eps'
    plt.savefig(fileout)
    fileout=FILESTART+'allplots/'+field+'/seasonal_cycle_all_models.pdf'
    plt.savefig(fileout)
    plt.close
    
    
    ax=plt.subplot(1,1,1)
    labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    for i in range(0,len(monmeans_NH[:,0,0])):
        if i < len(latmeans[:,0,0]) / 2.0:
           ax.plot(labels,monmeans_NH[i,0,:]-monmeans_NH[i,1,:],label=modelnames[i])
        else:
            ax.plot(labels,monmeans_NH[i,0,:]-monmeans_NH[i,1,:],label=modelnames[i],
                    linestyle='dashed')
        
        
    ax.plot(labels,np.mean(monmeans_NH[:,0,:],axis=0)-np.mean(monmeans_NH[:,1,:],axis=0),
        color='black',linestyle='dashed',linewidth=2,label='avg')
    
    print('monthmeans',np.mean(monmeans_NH[:,0,:],axis=0)-np.mean(monmeans_NH[:,1,:],axis=0))
   
    # plot pliomip1 data if appropriate
    if PLIOMIP1 == 'y':
        ax.plot(labels, mean_seas_pliomip1, color='black', linestyle='dotted',
                linewidth=2, label='PlioMIP1')
        ax.fill_between(labels, min_seas_pliomip1, max_seas_pliomip1, alpha=0.2, 
                        color="grey")
        
    
    plt.title(exptname+'-'+cntlname+': NH '+ field + ' anomaly')
    plt.ylabel(units)
    #plt.xlabel('month') 
    plt.figtext(0.02, 0.97,'a)',
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    plt.title(names.get(exptname) + '-'
              + names.get(cntlname) + ': ' 
              + namefield.get(field) + ' NH anomaly')
    
    # Shrink current axis by 20% and put a legend to the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
   
    fileout=FILESTART+'allplots/'+field+'/seasonal_cycle_all_models_anomaly.eps'
    plt.savefig(fileout)
    fileout=FILESTART+'allplots/'+field+'/seasonal_cycle_all_models_anomaly.pdf'
    plt.savefig(fileout)
    plt.close()
    
    # write out to a file
    if field == 'NearSurfaceTemperature':
        txtfile1.write("data_for_3a \n")
    elif field == 'TotalPrecipitation':
        txtfile1.write("data_for_6a \n")
       
        
    txtfile1.write('modelname, Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct, Nov, Dec \n')
    for i, model in enumerate(modelnames):
        txtfile1.write((model + ',' + 
                        np.str(np.around((monmeans_NH[i,0,0] - monmeans_NH[i,1,0]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,1] - monmeans_NH[i,1,1]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,2] - monmeans_NH[i,1,2]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,3] - monmeans_NH[i,1,3]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,4] - monmeans_NH[i,1,4]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,5] - monmeans_NH[i,1,5]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,6] - monmeans_NH[i,1,6]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,7] - monmeans_NH[i,1,7]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,8] - monmeans_NH[i,1,8]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,9] - monmeans_NH[i,1,9]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,10] - monmeans_NH[i,1,10]),2)) + ',' + 
                        np.str(np.around((monmeans_NH[i,0,11] - monmeans_NH[i,1,11]),2)) + '\n'))
        
    txtfile1.write('\n')
    txtfile1.close

    ###################################################################
    # plot the latitudinal range from each model
    
    # won't print this out as it doesn't look very useful.
    # absolute value of temperature by latitude
    #plt.subplot(2,1,1)
    #for i in range(0,len(latmeans[:,0,0])):
    #    # plot experiment data
    #   
    #    plt.plot(latmeans[i,0,:],lats,color='r')
    #    # plot control data
    #    plt.plot(latmeans[i,1,:],lats,color='b')
    #    plt.title('latitudinal average of '+field)
    #    plt.xlabel(units)
    #    
    #plt.subplot(2,1,2)
    #plt.plot(np.mean(latmeans[:,0,:],axis=0),lats,label=exptname,color='r')
    #plt.plot(np.mean(latmeans[:,1,:],axis=0),lats,label=cntlname,color='b')
    #plt.xlabel(units)
    #plt.ylabel('latitude')    
    #plt.legend()
    
    
    ax=plt.subplot(1,1,1)
    plt.figtext(0.02, 0.97,'c)',
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    for i in range(0,len(latmeans[:,0,0])):
        if i < len(latmeans[:,0,0]) / 2.0:
            ax.plot(latmeans[i,0,:]-latmeans[i,1,:],lats,label=modelnames[i])
        else:
            ax.plot(latmeans[i,0,:]-latmeans[i,1,:],lats,label=modelnames[i],
                    linestyle='dashed')
    ax.plot(np.mean(latmeans[:,0,:],axis=0)-np.mean(latmeans[:,1,:],axis=0),
             lats,label='avg',color='black',linestyle='dashed',linewidth=2)
    plt.title(names.get(exptname) + '-'
              + names.get(cntlname) + ': ' 
              + namefield.get(field) + ' anomaly')
    plt.xlabel(units)
    plt.ylabel('latitude') 
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    
    fileout=FILESTART+'allplots/'+field+'/latitude_anomaly.eps'
    plt.savefig(fileout)
    fileout=FILESTART+'allplots/'+field+'/latitude_anomaly.pdf'
    plt.savefig(fileout)
    plt.close()
  
     # write out to a file
    if field == 'NearSurfaceTemperature':
        txtfile1.write("data_for_1c \n")
    elif field == 'SST':
        txtfile1.write("zonal mean data for SST \n")
   
   
        for i, model in enumerate(modelnames):
            txtfile1.write(model + '\n')
            txtfile1.write('latitude, Tanom \n')
            for j, lat in enumerate(lats):
                txtfile1.write(np.str(np.around(lat,2)) + ',' + 
                               np.str(np.around(latmeans[i,0,j] - latmeans[i,1,j],2)) + '\n')
        
        txtfile1.write('\n')
        
    txtfile1.close
   
    #########################################################
    # plot a picture of the monthly anomaly from each model
    
    
    if field=='SurfaceTemperature':
            valmin=0.0
            valmax=10.5
            incr=0.5
            cmapname='Reds'
            
            datamin = -30.
            datamax = 35.
            dataincr = 5.
    if field=='NearSurfaceTemperature':
            valmin=0.0
            valmax=10.5
            incr=0.5
            cmapname='Reds'
            
            datamin = -30.
            datamax = 35.
            dataincr = 5.
    if field=='SST':
            valmin=0.0
            valmax=10.5
            incr=0.5
            cmapname='Reds'
            
            datamin = -5.
            datamax = 32.
            dataincr = 2.
    if field=='TotalPrecipitation':
            valmin=-1.4
            valmax=1.6
            incr=0.2
            cmapname='RdBu'
            
            datamin = 0.
            datamax = 5.
            dataincr = 0.1
    Vanom=np.arange(valmin,valmax,incr)
    Vdata= np.arange(datamin, datamax, dataincr)
    
   
    
    #########################################################
    # plot a picture of the change from each model and the 
    # multimodelmean
    
    
    anom_cubes=iris.cube.CubeList([])
    expt_cubes=iris.cube.CubeList([])
    cntl_cubes=iris.cube.CubeList([])
    for modelno in range(0,len(modelnames)):
        modeluse=modelnames[modelno]
        exptfile=FILESTART+modeluse+'/'+exptname+'.'+field+'.allmean.nc'
        cntlfile=FILESTART+modeluse+'/'+cntlname+'.'+field+'.allmean.nc'
        
        exptcube=iris.load_cube(exptfile)
        cntlcube=iris.load_cube(cntlfile)
       
       
        if modeluse == 'EC-Earth3.1' and field == 'SST':
           cntlcube.coord('latitude').bounds = None
           cntlcube.coord('longitude').bounds = None
           
        if modeluse == 'CCSM4-UoT' and field == 'TotalPrecipitation':
           cntlcube.coord('latitude').var_name = 'latitude'
           cntlcube.coord('longitude').var_name = 'longitude'
           exptcube.coord('latitude').var_name = 'latitude'
           exptcube.coord('longitude').var_name = 'longitude'
           cntlcube.coord('latitude').long_name = None
           cntlcube.coord('longitude').long_name = None
           exptcube.coord('latitude').long_name = None
           exptcube.coord('longitude').long_name = None
           cntlcube.coord('latitude').points = cntlcube.coord('latitude').points.astype('float32')
           cntlcube.coord('longitude').points = cntlcube.coord('longitude').points.astype('float32')
           exptcube.coord('latitude').points =  exptcube.coord('latitude').points .astype('float32')
           exptcube.coord('longitude').points = exptcube.coord('longitude').points.astype('float32')
        
        diffcube=exptcube-cntlcube
        
        # check float 32 for concatenation
        diffcube.data=diffcube.data.astype('float32') 
        exptcube.data=exptcube.data.astype('float32') 
        cntlcube.data=cntlcube.data.astype('float32') 
        if field=='NearSurfaceTemperature' or field == 'SST':
            if (modeluse=='MIROC4m' or modeluse=='COSMOS'
                or (modeluse == 'CCSM4-Utr' and field =='SST')):
                diffcube.units='Celsius'
                exptcube.units='Celsius'
                cntlcube.units='Celsius'
            else:
                diffcube.convert_units('Celsius')
                exptcube.convert_units('Celsius')
                cntlcube.convert_units('Celsius')
        
        # remove scalar coordinates so that we can concatenate
        # also add a new axis with the model number
        
        newcube = resort_coords(diffcube,modelno)
        newdata = newcube.data
        anom_cubes.append(newcube)    
    
        newcube = resort_coords(exptcube,modelno)
        newcube.rename(field)
        print(modeluse)
        print('MINIMUM',np.min(newcube.data))
        expt_cubes.append(newcube)
        
        newcube = resort_coords(cntlcube,modelno)
        newcube.rename(field)
        cntl_cubes.append(newcube)          
         
        
        # plot individual values if required
        if individual_plot=='y':
            plt.subplot(1,2,1)
            cs=iplt.contourf(exptcube,Vdata,extend='both',
                             cmap='terrain')
            titlename=modeluse+' '+exptname+': '+field
            cbar=plt.colorbar(cs,orientation="horizontal")
            cbar.set_label(units)
            cbar.ax.tick_params(labelsize=8, labelrotation=60) 
            plt.title(titlename,fontsize=8)
            plt.gca().coastlines()
        
            plt.subplot(1,2,2)
            cs=iplt.contourf(cntlcube,Vdata,extend='both',
                             cmap='terrain')
            titlename=modeluse+' '+cntlname+': '+field
            cbar=plt.colorbar(cs,orientation="horizontal")
            cbar.ax.tick_params(labelsize=8, labelrotation=60) 
            cbar.set_label(units)
            plt.title(titlename,fontsize=8)
            plt.gca().coastlines()
       
            fileout=FILESTART+'allplots/'+field+'/individual_models/'+modeluse+'_'+field+'.eps'
            plt.savefig(fileout)
            fileout=FILESTART+'allplots/'+field+'/individual_models/'+modeluse+'_'+field+'.pdf'
            plt.savefig(fileout)
            plt.close()
        
            # plot a picture of the anomaly from each of the models
            plt.subplot(1,1,1)
            
            
            cs=iplt.contourf(exptcube-cntlcube,Vanom,extend='both',cmap=cmapname)
            titlename=modeluse+' '+exptname+'-'+cntlname+': '+field
            cbar=plt.colorbar(cs,orientation="horizontal")
            cbar.set_label(units)
            plt.title(titlename,fontsize=8)
            plt.gca().coastlines()
            fileout=FILESTART+'allplots/'+field+'/individual_models/'+modeluse+'_'+field+'_anomaly.eps'
            plt.savefig(fileout)
            fileout=FILESTART+'allplots/'+field+'/individual_models/'+modeluse+'_'+field+'_anomaly.pdf'
            plt.savefig(fileout)
            plt.close()
        
   
       
    #############################    
    # get the multi-modelmean and standard deviation
    
   
    iris.experimental.equalise_cubes.equalise_attributes(expt_cubes)
    iris.experimental.equalise_cubes.equalise_attributes(cntl_cubes)
    allpliocube = expt_cubes.concatenate_cube()
   
    meanplio = allpliocube.collapsed(['model_level_number'], iris.analysis.MEAN)
    meanplio.rename(field + 'mean_mPWP')
    maxplio = allpliocube.collapsed(['model_level_number'], iris.analysis.MAX)
    maxplio.rename(field + 'max_mPWP')
    minplio = allpliocube.collapsed(['model_level_number'], iris.analysis.MIN)
    minplio.rename(field + 'min_mPWP')
    stdplio = allpliocube.collapsed(['model_level_number'], iris.analysis.STD_DEV)
    stdplio.rename(field + 'std_mPWP')
    
    allpicube = cntl_cubes.concatenate_cube()
    meanpi = allpicube.collapsed(['model_level_number'], iris.analysis.MEAN)
    meanpi.rename(field + 'mean_pi')
    maxpi = allpicube.collapsed(['model_level_number'], iris.analysis.MAX)
    maxpi.rename(field + 'max_pi')
    minpi = allpicube.collapsed(['model_level_number'], iris.analysis.MIN)
    minpi.rename(field + 'min_pi')
    stdpi = allpicube.collapsed(['model_level_number'], iris.analysis.STD_DEV)
    stdpi.rename(field + 'std_pi')
    
    allmeancube=anom_cubes.concatenate_cube() # all model mean anomalies
    meancube=allmeancube.collapsed(['model_level_number'], iris.analysis.MEAN)
    meancube.rename(field + 'mean_anomaly')
    maxcube=allmeancube.collapsed(['model_level_number'], iris.analysis.MAX)
    maxcube.rename(field + 'max_anomaly')
    mincube=allmeancube.collapsed(['model_level_number'], iris.analysis.MIN)
    mincube.rename(field + 'min_anomaly')
    stdcube=allmeancube.collapsed(['model_level_number'], iris.analysis.STD_DEV)
    stdcube.rename(field + 'anomaly_multimodel_stddev')
    
    ###############################################
    # write out the mean and standard deviation to a netcdf file
    
    cubelist = iris.cube.CubeList([meanplio, stdplio, maxplio, minplio,
                                   meanpi, stdpi, maxpi, minpi,
                                   meancube, stdcube, maxcube, mincube])
    fileout = (FILESTART + field + '_multimodelmean.nc')
    iris.save(cubelist,fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=1.0E20)

    
    
    ###########################
    # plot the mean value and standard deviation
    
   
    ax = plt.axes(projection = ccrs.PlateCarree())
    V=np.arange(valmin,valmax,incr)
    mycmap = plt.cm.get_cmap(cmapname,len(V+2))
    qplt.contourf(meancube, V,extend='both',cmap=mycmap)
    if linux_win == 'l':
        exptlsm = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    else:
        exptlsm = FILESTART + 'PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
    lsmcube=iris.load_cube(exptlsm)
    qplt.contour(lsmcube,1,colors='black') 
    plt.figtext(0.02, 0.97,'b)',
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    plt.title(namefield.get(field) +' anomaly: multimodel mean')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-', draw_labels='true')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    fileout=(FILESTART+'allplots/'+field+'/multimodelmean.eps')
    plt.savefig(fileout)
    fileout=(FILESTART+'allplots/'+field+'/multimodelmean.pdf')
    plt.savefig(fileout)
    plt.close()
    
    if field=='NearSurfaceTemperature':
        V=np.arange(0,5.0,incr)
        textout = 'd)'
    if field=='TotalPrecipitation':
        textout = 'c)'
        V=np.arange(0,1.3, 0.1)
    ax = plt.axes(projection = ccrs.PlateCarree())
    qplt.contourf(stdcube, V,extend='both',cmap='plasma')
    lsmcube=iris.load_cube(exptlsm)
    qplt.contour(lsmcube,1,colors='black') 
    plt.title(namefield.get(field) +' anomaly:Standard Deviation')
    plt.figtext(0.02, 0.97,textout,
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-', draw_labels='true')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    fileout=(FILESTART+'allplots/'+field+'/multimodelstdev.eps')
    plt.savefig(fileout)
    fileout=(FILESTART+'allplots/'+field+'/multimodelstdev.pdf')
    plt.savefig(fileout)
    plt.close()
     

    cubelist = iris.cube.CubeList([meancube, stdcube])
    iris.save(cubelist,FILEOUTNC,netcdf_format='NETCDF3_CLASSIC',fill_value=1.0E20)           
        
        
    ##########################################################
    # plot the monthly anomaly from each model
    if linux_win=='l' and individual_plot=='y':
        for modelno in range(0,len(modelnames)):
            modeluse=modelnames[modelno]
            exptfile=FILESTART+modeluse+'/'+exptname+'.'+field+'.mean_month.nc'
            cntlfile=FILESTART+modeluse+'/'+cntlname+'.'+field+'.mean_month.nc'
            exptcube=iris.load_cube(exptfile)
            cntlcube=iris.load_cube(cntlfile)
     
            for mon in range(0,12):
                anom=exptcube.data[mon,:,:]-cntlcube.data[mon,:,:]
                lat=exptcube.coord('latitude').points
                lon=exptcube.coord('longitude').points
                lons,lats=np.meshgrid(lon,lat)
                map=Basemap(llcrnrlon=0.0,urcrnrlon=360.0,
                         llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',
                         resolution='c')
                x, y = map(lons, lats)
                map.drawcoastlines()
                V=np.arange(valmin,valmax,incr)
                cs = map.contourf(x,y,anom,V,cmap=cmapname,extend="both")
                cbar = plt.colorbar(cs,orientation="horizontal")
                plt.title(modeluse+':'+'month is '+str(mon+1))
                cbar.set_label(units)
                fileout=(FILESTART+'allplots/'+field+'/global_months/'+modeluse+'_'+field+
                     '_anomaly'+str(mon+1)+'.eps')
                plt.savefig(fileout)
                fileout=(FILESTART+'allplots/'+field+'/global_months/'+modeluse+'_'+field+
                     '_anomaly'+str(mon+1)+'.pdf')
                plt.savefig(fileout)
                plt.close()
          
        
  


##########################################################
# main program
        
filename=' '
linux_win='l'

modelnames=['CESM2', 'HadGEM3','IPSLCM6A', 'COSMOS', 
            'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
            'MIROC4m', 'IPSLCM5A2', 'HadCM3',
            'GISS2.1G', 'CCSM4', 
            'CCSM4-Utr', 'CCSM4-UoT', 
            'NorESM-L', 'MRI2.3', 'NorESM1-F'
            ]

#modelnames = ['CCSM4-Utr']
            
PLIOMIP1 = 'n'
#modelnames=['COSMOS',
#            'CCSM4-UoT']   

fieldnames=['TotalPrecipitation']
units=['mm/day']
#fieldnames=['NearSurfaceTemperature']
#units=['degC']
#fieldnames=['SST']
#units=['degC']
exptname='EOI400'
cntlname='E280'
individual_plot='n' # do you want to plot the anomalies for all of the individual models

if linux_win=='w':
   FILESTART='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded/'
else:
   FILESTART='/nfs/hera1/earjcti/regridded/'

FILEOUT = FILESTART + 'dummy.txt'
FILEOUTNC = FILESTART + 'dummy.nc'

for field in range(0,len(fieldnames)):

    if fieldnames[field] == 'TotalPrecipitation':
        FILEOUT = FILESTART + 'alldata/data_for_5a_6a.txt'
        FILEOUTNC = FILESTART + 'alldata/data_for_5b_5c.nc'
    if fieldnames[field] == 'NearSurfaceTemperature':
        FILEOUT = FILESTART + 'alldata/data_for_1a_1c_3a.txt'
        FILEOUTNC = FILESTART + 'alldata/data_for_1b_1d.nc'
    if fieldnames[field] == 'SST':
        FILEOUT = FILESTART + 'alldata/zonal_mean_SST_data.txt'
  
   
    plotmean(modelnames,fieldnames[field],exptname,cntlname,linux_win,units[field],individual_plot)

#sys.exit(0)
\
