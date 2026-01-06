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
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
import cf_units as unit
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#os.environ["PROJ_LIB"] = r'C:\Users\julia\Miniconda2\pkgs\proj4-5.2.0-hc56fc5f_1003\Library\share'
#from mpl_toolkits.basemap import Basemap, shiftgrid
import sys



#####################################
def  climate_sensitivity_analysis(modelnames,fieldname,exptname,cntlname,linux_win,units):
   
    if linux_win=='w':
        filestart='C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
        datatext = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\data_for_7a-b.txt'
        netcdfout = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\data_for_7c.nc'
    else:
        filestart='/nfs/hera1/earjcti/regridded/'
        datatext = '/nfs/hera1/earjcti/regridded/alldata/data_for_7a-b.txt'
        netcdfout = '/nfs/hera1/earjcti/regridded/alldata/data_for_7c.nc'
     
        
    # set up a dictionary for the climate sensitivity
   
    # from my stuff
    #clim_sens ={'NorESM-L': 3.1,
    #             'NorESM1-F':2.29,
    #             'IPSLCM6A': 4.8,
    #             'IPSLCM5A':3.4,
    #             'HadCM3': 3.7,
    #             'MIROC4m':3.9,
    #             'COSMOS':4.1,
    #             'UofT':3.8,
    #             'EC-Earth3.1':3.2,
    #             'MRI-CGCM2.3':2.8,
    #             'CESM1.0.5': 3.1,
    #             'GISS': 3.31
    #             }
        
    # from Alan's table provided by authors
    clim_sens ={'NorESM-L': 3.1,
                 'NorESM1-F':2.3,
                 'IPSLCM6A': 4.8,
                 'IPSLCM5A2':3.6,
                 'IPSLCM5A':4.1,
                 'HadCM3': 3.5,
                 'MIROC4m':3.9,
                 'COSMOS':4.7,
                 'CCSM4-UoT':3.2,
                 'EC-Earth3.1':3.2,
                 'EC-Earth3.3':4.3,
                 'MRI2.3':2.8, # from my investigation
                 'CCSM4-Utr': 3.2,
                 'GISS2.1G': 3.3,
                 'CESM2': 5.3,
                 'CESM1.2' :4.1,
                 'CCSM4' :3.2,
                 'CCSM4-avg' : 3.2
                 }
        
     
    # first get the data.  We need climate sensitivity, 
    # global temperature anomaly, latitude temperature anomaly
    # and gridbox by gridbox temperature anomaly
    climdiff=np.zeros(len(modelnames))
    climdiffmon=np.zeros((len(modelnames),12))
    climdifflat=np.zeros((len(modelnames),180))
    sensitivity_array=np.zeros(len(modelnames))
    alllats=np.arange(0,180,1)-89.5
    
    for mod in range(0,len(modelnames)):
        sensitivity_array[mod]=clim_sens.get(modelnames[mod])
        
        # get data from experiment file
        fileexpt=filestart+modelnames[mod]+'/'+exptname+'.'+fieldname+'.data.txt'
        file1= open(fileexpt,"r")
        lines=list(file1)
       
        
        meanexpt,sdexpt=lines[2].split(",")
        monmeanexpt=np.zeros(12)
        latmeanexpt=np.zeros(180)
        for l in range(5,17):
            index,mean,sd=lines[l].split(",")
            monmeanexpt[np.int(index)-1]=np.float(mean)
        for l in range(20,200):
            lat,mean,sd=lines[l].split(",")
            index=np.where(alllats==np.float(lat))
            latmeanexpt[index]=np.float(mean)
       
        
        filecntl=filestart+modelnames[mod]+'/'+cntlname+'.'+fieldname+'.data.txt'
        file2= open(filecntl,"r")
        lines=list(file2)
        meancntl,sdexpt=lines[2].split(",")
        monmeancntl=np.zeros(12)
        latmeancntl=np.zeros(180)
        for l in range(5,17):
            index,mean,sd=lines[l].split(",")
            monmeancntl[np.int(index)-1]=np.float(mean)
        for l in range(20,200):
            lat,mean,sd=lines[l].split(",")
            index=np.where(alllats==np.float(lat))
            latmeancntl[index]=np.float(mean)
        
        climdiff[mod]=np.float(meanexpt)-np.float(meancntl)
        
        climdiffmon[mod,:]=monmeanexpt-monmeancntl
        climdifflat[mod,:]=latmeanexpt-latmeancntl
        #print(modelnames[mod],climdifflat[mod])
    
     
    ########################################################
    # plot the climate sensitivity vs the global mean
   
    fig, ax1 = plt.subplots(figsize=(7.0, 4.0))
    ax1.plot(climdiff,sensitivity_array,'x')
    ax1.set_xlabel('Plio_Core - PI_CTL SAT anomaly',fontsize=15)
    ax1.set_ylabel('ECS', fontsize=15)
    ax1.set_xlim(np.floor(np.min(climdiff)),np.ceil(np.max(climdiff)))
    
    # do a linear regression
    print(modelnames)
    print(climdiff)
    print(sensitivity_array)
   
    slope, intercept, r_value, p_value, std_err = sp.stats.linregress(climdiff, sensitivity_array)
    xarray=np.arange(0,10,1)
    yarray=intercept+(slope*xarray)
    ax1.plot(xarray,yarray)
    ax1.tick_params(axis='x',  labelsize=15)
    ax1.tick_params(axis='y',  labelsize=15)
    #plt.title("R-squared: " + np.str(np.around((r_value**2.), 2)) 
    #        + ",  p-value: " + np.str(np.around(p_value, 2)) , fontsize=15) 
    
    print('julia',slope, intercept, r_value, p_value)
   # sys.exit(0)
    figtext = 'a)'
    plt.figtext(0.02, 0.97,figtext,
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0+(0.1*box.height), box.width, box.height*0.9])
  
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_globalanom.png')
    
    #plt.tight_layout()
    plt.savefig(fileout)
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_globalanom.pdf')
    plt.savefig(fileout)
    plt.close()

    rsq_std=r_value**2.
    
    f1 = openx=open(datatext,'w')
    f1.write('Data for Figure 7a\n')
    f1.write('model name, Plio_core - PI_Cntl, ECS\n')
    for i in range(0,len(modelnames)):
        f1.write(modelnames[i] + ',' + np.str(np.round(climdiff[i],2)) + ',' + 
                 np.str(np.round(sensitivity_array[i],2)) + '\n')
   

    
    ########################################################
    # plot the correlation between climate sensitivity vs the monthly mean
    rvals=np.zeros(12)
    pvals=np.zeros(12)
    for mon in range(0,12):
        slope,intercept, r_value, p_value, std_err = (
                sp.stats.linregress(climdiffmon[:,mon], sensitivity_array))
        rvals[mon]=r_value**2.
        pvals[mon]=p_value
    
    fig, ax1 = plt.subplots(figsize=(7.0, 4.0))
    labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    color = 'tab:red'
    #ax1.set_xlabel('month', fontsize=15)
    ax1.set_ylabel('Rsq', color=color, fontsize=12)
    ax1.plot(labels, rvals, color=color)
    ax1.tick_params(axis='y', labelcolor= color, labelsize=12)
    ax1.tick_params(axis='x',  labelsize=12)
    #ax1.plot([0,13],[rsq_std,rsq_std],color='black',linestyle='dashed',linewidth=2)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0+(0.1*box.height), box.width, box.height*0.9])

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axi
    ax2.set_position([box.x0, box.y0+(0.1*box.height), box.width, box.height*0.9])

    color = 'tab:blue'
    ax2.set_ylabel('p-value', color=color, fontsize=12)  # we already handled the x-label with ax1
    ax2.plot(labels, pvals, color=color)
    ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
    #plt.title('ECS vs Plio_Core - PI_CTL by month', fontsize=15)
    
    
    figtext = 'b)'
    plt.figtext(0.02, 0.97,figtext,
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
   
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_monanom.png')
    plt.savefig(fileout)
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_monanom.pdf')
    plt.savefig(fileout)
    plt.close()
    
    
    
     ########################################################
    # plot the correlation between climate sensitivity vs the latitudinal mean
    rvals=np.zeros(180)
    pvals=np.zeros(180)
    for lat in range(0,180):
        slope,intercept, r_value, p_value, std_err = (
                sp.stats.linregress(climdifflat[:,lat], sensitivity_array))
        rvals[lat]=r_value**2.
        pvals[lat]=p_value
    
    fig, ax1 = plt.subplots(figsize=(7.0, 4.0))

    color = 'tab:red'
    ax1.set_xlabel('latitude', fontsize=12)
    ax1.set_ylabel('Rsq', color=color, fontsize=12)
    ax1.plot(alllats, rvals, color=color)
    ax1.tick_params(axis='y', labelcolor= color, labelsize=12)
    ax1.tick_params(axis='x', labelsize=12)
    #ax1.plot([-90,90],[rsq_std,rsq_std],color='black',linestyle='dashed',linewidth=2)
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0+(0.1*box.height), box.width, box.height*0.9])

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.set_position([box.x0, box.y0+(0.1*box.height), box.width, box.height*0.9])
    color = 'tab:blue'
    ax2.set_ylabel('p-value', color=color, fontsize=12)  # we already handled the x-label with ax1
    ax2.plot(alllats, pvals, color=color)
    #for i, lat in enumerate(alllats):
    #    print(lat, pvals[i], rvals[i])
   
    ax2.tick_params(axis='y', labelcolor=color, labelsize=12)
    #ax2.set_ylim(0,0.1)
    #plt.title('ECS vs Plio_Core - PI_CTL by latitude', fontsize=15)
    
    figtext = 'b)'
    plt.figtext(0.02, 0.97,figtext,
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    
    
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_latanom.png')
    plt.savefig(fileout)
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_latanom.pdf')
    plt.savefig(fileout)
    
    f1.write('\n')
    f1.write('Data for Figure 7b\n')
    f1.write('latitude, Pvalue, Rsqvalue\n')
    for i in range(0,len(alllats)):
        f1.write(np.str(np.round(alllats[i],1)) + ',' + np.str(np.round(pvals[i],3)) + ',' + 
                 np.str(np.round(rvals[i],2)) + '\n')
    f1.close()
    
    #############################################################################
    # now get the global data and do a correlation
    
    cubelist = iris.cube.CubeList([])
    for mod in range(0,len(modelnames)):
        # get average anomaly
        
        fileexpt=filestart+modelnames[mod]+'/'+exptname+'.'+fieldname+'.allmean.nc'
        exptcube=iris.load_cube(fileexpt)
        filecntl=filestart+modelnames[mod]+'/'+cntlname+'.'+fieldname+'.allmean.nc'
        cntlcube=iris.load_cube(filecntl)
        
        ny,nx=np.shape(exptcube.data)
        if mod==0:
            anommap=np.zeros((len(modelnames),ny,nx))
           
       
        anommap[mod,:,:]=exptcube.data-cntlcube.data
    
        
    
    rsqmap=np.zeros((ny,nx))
    pvalmap = np.zeros((ny, nx))
    slopemap = np.zeros((ny,nx))
    interceptmap = np.zeros((ny,nx))
    for j in range(0,ny):
         for i in range(0,nx):
             slope,intercept, r_value, p_value, std_err = (
                sp.stats.linregress(anommap[:,j,i], sensitivity_array))
             rsqmap[j,i]=r_value**2.
             pvalmap[j,i] = p_value
             slopemap[j,i] = slope
             interceptmap[j,i] = intercept
           
             
    rsqmapcube=exptcube.copy(data=rsqmap) 
    rsqmapcube.units=None
    rsqmapcube.long_name = 'Rsq'
    rsqmapcube.standard_name = None
    rsqmapcube.var_name = 'Rsq'

    
    slopecube=exptcube.copy(data=slopemap) 
    slopecube.units=None
    
    interceptcube=exptcube.copy(data=interceptmap) 
    interceptcube.units=None
    
    temparr = np.where(pvalmap < 0.05, 1, 0) 
    significance_cube = rsqmapcube.copy(data=temparr)
    pval_cube = rsqmapcube.copy(data=pvalmap)
    pval_cube.units = None
    pval_cube.long_name = 'pvalue'
    pval_cube.standard_name = None
    pval_cube.var_name = 'pvalue'
    
    cubelist.append(rsqmapcube)
    cubelist.append(pval_cube)
    print(cubelist)
    iris.save(cubelist, netcdfout, netcdf_format='NETCDF3_CLASSIC')

    
    # plot the map with Rsq and the significance
    
    fig = plt.subplots(figsize=(7.0, 5.0))
    ax = plt.axes(projection = ccrs.PlateCarree())
    V=np.arange(0.0,1,0.05)
    
    qplt.contourf(rsqmapcube,V,cmap='YlGnBu')
    iplt.contourf(significance_cube, 1, hatches=[None, '///'], colors='none')
    iplt.contourf(significance_cube, 1, hatches=[None, '\\\''], colors='none')
    #titlename=modeluse+' '+exptname+': '+field
    #bar=plt.colorbar(cs,orientation="horizontal")
    #plt.title('correlation between climate sensitivity and \n mPWP-PI anomaly '
    #          'at this point')
    plt.title('')
    gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-', draw_labels='true')
    gl.xlabels_top = False
    gl.ylabels_left = False
    gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    
    #plt.title(titlename,fontsize=8)
    plt.gca().coastlines()
 
    figtext = 'c)'
    plt.figtext(0.02, 0.97,figtext,
     horizontalalignment='left',
     verticalalignment='top',
     fontsize=20)
    
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_globe.png')
    plt.savefig(fileout)
    fileout=(filestart + 'allplots/' + 
             fieldname + '/climate_sensitivity_vs_globe.pdf')
    plt.savefig(fileout)  
    plt.close()    
    
        
    # print out the slope and the intercept at each gridpoint
    plt.subplot(1,1,1)
    qplt.contourf(slopecube)
    plt.title('slope between climate sensitivity and \n mPWP-PI anomaly '
              'at this point')
    plt.gca().coastlines()
    plt.savefig(filestart + 'allplots/' + 
                fieldname + '/climate_sensitivity_map_slope.eps')  
 
    plt.close()
    
    plt.subplot(1,1,1)
    
    qplt.contourf(interceptcube)
    plt.title('intercept between climate sensitivity and \n mPWP-PI anomaly '
              'at this point')
    plt.gca().coastlines()
    plt.savefig(filestart + 'allplots/' + 
                fieldname + '/climate_sensitivity_map_intercept.eps')  
 
    plt.close()
    
    outfile = (filestart + 'allplots/' + 
              fieldname + '/climate_sensitivity_relationships.txt')
    
    txtfile1 = open(outfile,"w+") 
    txtfile1.write("longitude latitude rsq pvalue intercept slope \n")
    lons = interceptcube.coord('longitude').points
    lats = interceptcube.coord('latitude').points
    
    for j in range(0,ny):
         for i in range(0,nx):
             writestring = (np.str(np.around(lons[i],2)) + ',' + 
                            np.str(np.around(lats[j],2)) + ',' + 
                            np.str(np.around(rsqmap[j, i],2)) + ',' +
                            np.str(np.around(pvalmap[j, i],2)) + ',' +
                            np.str(np.around(interceptmap[j, i],2)) + ',' +
                            np.str(np.around(slopemap[j,i],2)) + '\n')
             txtfile1.write(writestring)
    txtfile1.close

##########################################################
# main program
        
filename=' '
linux_win='l'
#modelnames=['MIROC4m','COSMOS']   # MIROC4m  COSMOS UofT EC-Earth3.1

modelnames=['CCSM4-Utr','COSMOS', 'CCSM4',
            'EC-Earth3.3', 'CESM1.2','CESM2',
            'GISS2.1G','HadCM3',
            'IPSLCM6A','IPSLCM5A2','IPSLCM5A',
            'MIROC4m','MRI2.3',
            'NorESM-L','NorESM1-F',
            'CCSM4-UoT'
            ]

#fieldnames=['TotalPrecipitation']
#units=['mm/day']
fieldnames=['NearSurfaceTemperature']
#fieldnames = ['SST']
units=['degC']
exptname='EOI400'
cntlname='E280'

for field in range(0,len(fieldnames)):
    climate_sensitivity_analysis(modelnames,fieldnames[field],exptname,cntlname,linux_win,units[field])

#sys.exit(0)
