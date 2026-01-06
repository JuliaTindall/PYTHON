# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 13:37:24 2019

@author: julia
This program will estimate the climate sensitivity from the proxy data as follows:
    1.  read in proxy data
    2.  read in the file from the model which see's whether there is a significant relationship
        between Plio(Tanom) and ECS at each gridbox
    3.  For each proxy point
        a) check if there should be a significant relationship
        b) if so estimate the climate sensitivity using the slope and the intercept
        c) plot a map of all the climate sensitivities
        d) print out the range of all the climate sensitivities
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import sys
import os

#os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid




def main():
    """ 
    This program will estimate the climate sensitivity from the proxy data as follows:
    1.  read in proxy data
    2.  read in the file from the model which see's whether there is a significant relationship
        between Plio(Tanom) and ECS at each gridbox
    3.  For each proxy point
        a) check if there should be a significant relationship
        b) if so estimate the climate sensitivity using the slope and the intercept
        c) print out the range of all the climate sensitivities
    4.  Read in the cube showing data from figure 7d.  Which shows the 
        p value at each gridcell
    
    5. plot a map of regions where there is a significant relationship.
    6. as 5. but with overplot the climate sensitivties derived from each point
    """
   
    
    #1. read in proxy data
    proxylat, proxylon, proxy_sst_anom = readproxy()
    gridlat, gridlon,  pval, intercept, slope = readfile()
    
    nproxies = len(proxylat)
    ngrids = len(gridlat)
    clim_sens = np.zeros(nproxies)
    
    # 2. 3. check significance and estimate climate sensitivity
    for i in range(0,nproxies):
        # get the subscript from the model relationship file
        grid_ss, griddiff = get_subscript(proxylat[i], proxylon[i], gridlat, gridlon, ngrids)
        # see if it is significant (p < 0.05
        if pval[grid_ss] < 0.05:
           # if significant CS = intercept + (proxy_sst_anom) * slope
            clim_sens[i] = intercept[grid_ss] + (slope[grid_ss] * proxy_sst_anom[i])
            print('ind sens',proxylat[i],proxylon[i],clim_sens[i])
        else:
            clim_sens[i] = np.nan
        # print cs
      
    # 4. find regions that there is a significant relationship
    pval_cube = iris.load_cube(SIGNIFICANCE_FILE, SIGNIFICANCE_NAME)
    sign_data = np.where(pval_cube.data < 0.05, 1.0, 0.0)
    sign_cube=pval_cube.copy(data=sign_data)
    
    
    # put into a reduced array and plot
    nvals = np.count_nonzero(~np.isnan(clim_sens))
    print(nvals)
    count=0
    latredu = np.zeros(nvals)
    lonredu = np.zeros(nvals)
    sstanomredu = np.zeros(nvals)
    clim_sens_redu = np.zeros(nvals)
    
    for i in range(0,nproxies):
        if np.isfinite(clim_sens[i]):
            latredu[count] = proxylat[i]
            if proxylon[i] > 180:
                lonredu[count] = proxylon[i]-360.
            else:
                lonredu[count] = proxylon[i]
            clim_sens_redu[count] = clim_sens[i]
            sstanomredu[count] = proxy_sst_anom[i]
            count=count+1
            
    plotdata(latredu,lonredu,clim_sens_redu,nvals, FILEOUT)
    
    # remove sites which we are not sure about.
    # this is where the datapoint is not within 1 deg of the modelled range
    # ie the data does not even nearly agree with any of the models
    # ie where data and model do not agree at all
    
    (new_latredu, new_lonredu, 
     new_nvals, new_clim_sens_redu) = redu_sites(latredu, 
                                                 lonredu, 
                                                 sstanomredu,
                                                 clim_sens_redu,
                                                 proxy_sst_anom,
                                                 proxylat, proxylon)  
                                                 
    plotdata(new_latredu,new_lonredu,new_clim_sens_redu,new_nvals, FILEOUT_R)
    print('new clim_sens_redu',new_clim_sens_redu)
    sys.exit(0)
    
    
    # now plpot data but overplot where there is a significant relationship
    plotdata_overplot(new_latredu,new_lonredu,
                      new_clim_sens_redu,new_nvals, 
                      pval_cube, sign_cube,
                      FILEOUT_S)
    
    
def readfile():
    """
    reads data from the file
    returns numpy arrays of, lat, lon, pval, intercept, slope
    """
    f1 = open(FILECS,'r') # to count lines
    count=0
    for line in f1.readlines():
        count = count + 1
    f1.close()
     
    nvals = count 
    lats = np.zeros(nvals)
    lons = np.zeros(nvals)
    intercepts = np.zeros(nvals)
    pvals = np.zeros(nvals)
    slopes = np.zeros(nvals)
    
    f1 = open(FILECS,'r') # to read
    count=0
    for line in f1.readlines():
        if line[0:4] == 'long': # titleline ignore
            print('titleline is',line)
            pass
        else:
            vals = line.split(',')
            lons[count] = vals[0]
            lats[count] = vals[1]
            pvals[count] = vals[3]
            intercepts[count] = vals[4]
            slopes[count] = vals[5]
            
        count = count + 1
    f1.close()
    
    return lats, lons, pvals, intercepts, slopes
   
def readproxy():
    """
    reads in the excel spreadsheet of the proxy dataset
    returns arrays of the latitude and longitude and sstanom of the proxy dataset
    """
    
    lats = np.zeros(37)
    lons = np.zeros(37)
    sstanom  = np.zeros(37)
    f1 = open(FILEPROXY,'r') # to read
    count=0
    for line in f1.readlines():
        if count >= 37:
            pass
        else:
            if line[0:4] == 'Loca': # titleline ignore
                print('titleline is',line)
                pass
            else:
                vals = line.split(',')
             
                if np.float(vals[2]) > 0:
                    lons[count] = np.float(vals[2])
                else:
                    lons[count] = np.float(vals[2]) + 360.
                lats[count] = np.float(vals[1])
                sstanom[count] = np.float(vals[15]) # this is plio - noaa
                count = count + 1
          
            
       
    f1.close()
   
    return lats, lons, sstanom

def readmodel():
    """
    reads in the excel spreadsheet of the model dataset
    returns arrays of the latitude and longitude, the minimum modelled ssta and
    the maximum modelled ssta
    """
    
    lats = np.zeros(37)
    lons = np.zeros(37)
    sstanom_min  = np.zeros(37)
    sstanom_max = np.zeros(37)
    allsstanom = np.zeros((37, 17)) # the 17th is the MMM
    f1 = open(FILEMODEL,'r') # to read
    count=0
    lines = f1.readlines()
    line = lines[0]
    titles = line.split(',')
    modnames = titles[3:20]
    
    # assume data starts at row index 1
    for i in range(1,len(lines)):
      line = lines[i]
      vals = line.split(',')
      lats[i-1] = vals[1]
      lons[i-1] = vals[2]
      print(i, vals[3:19])
      sstanom_min[i-1] = np.min(np.asarray(vals[3:19], dtype = float))
      sstanom_max[i-1] = np.max(np.asarray(vals[3:19], dtype = float))
      allsstanom[i-1, :] = np.asarray(vals[3:20], dtype = float)
    

   
    return lats, lons, sstanom_min, sstanom_max, allsstanom, modnames

def print_rmse(mod_allsst, proxy_allsst, modlats, modlons, 
               proxylats, proxylons, modnames):
    """
    just prints out the rmse for all the models

    """
    print('shape proxy', np.shape(proxy_allsst))
    print('shape data', np.shape(mod_allsst))
    npoints, nmods = np.shape(mod_allsst)
    
    f1 = open('ind_model_dmc.txt','w+')
    f1.write('model        rmse      bias   within 2deg/ 1deg/ 0.5deg \n')
    for i in range(0, nmods):
        sumsq = 0.0
        avger = 0.0
        count = 0.0
        within_1deg = 0
        within_2deg = 0
        within_05deg = 0
        for j in range(0, npoints):
            sumsq = sumsq + ((proxy_allsst[j] - mod_allsst[j, i])**2.0)
            avger = avger + (proxy_allsst[j] - mod_allsst[j, i])
            if (np.abs(proxy_allsst[j] - mod_allsst[j, i]) < 1.0):
                within_1deg = within_1deg + 1
            if (np.abs(proxy_allsst[j] - mod_allsst[j, i]) < 2.0):
                within_2deg = within_2deg + 1
            if (np.abs(proxy_allsst[j] - mod_allsst[j, i]) < 0.5):
                within_05deg = within_05deg + 1
            count = count + 1.0
            
        rmse = np.sqrt(sumsq / count)
        avger = avger / count
      
        f1.write(modnames[i].ljust(12) + ',' +  np.str(np.around(rmse,2))
                 + ',   ' +  np.str(np.around(avger,2)) + ',   ' 
                 + np.str(within_2deg) +  ',' +  np.str(within_1deg)
                 +  ',' +  np.str(within_05deg) + '\n')
        print(modnames[i].ljust(12),',',np.around(rmse,2),
                 ',   ',np.around(avger,2),',   ',within_2deg,  
               '   ',within_1deg,',   ',within_05deg)
        
    f1.close()
   

def get_subscript(latreq, lonreq, gridlat_arr, gridlon_arr, ngrid):
    """
    this program is passed a latitude and longitude (latreq, lonreq)
    and also two array containing (ngrid) values.  The arrays each contain
    latitudes and longitudes
    we want to find the subscript of the array that most closely matches the
    required values and return it
    """
    
    diffvals = 100.
    subscript = 0
    
    for i in range(0, ngrid):
        thisdiff = np.abs(gridlat_arr[i] - latreq) + np.abs(gridlon_arr[i] - lonreq)
        if thisdiff < diffvals:
            diffvals = thisdiff
            subscript = [i]
    

    return subscript, diffvals

def plotdata(lat,lon,clim_sens,nproxies,fileout):
    """
    plots the cliate sensitivity on a map
    """
    m=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,
              urcrnrlat=90.0,projection='cyl',resolution='c')
    m.drawmapboundary
    m.drawcoastlines()
    parallels=np.arange(-90.,90.,50.)
    m.drawparallels(parallels,labels=[True,True,True,True],fontsize=10) # labels right
    meridians=np.arange(-180.,180.,60.)
    m.drawmeridians(meridians,labels=[True,True,True,True],fontsize=10)
    
   
    x1,y1=m(lon,lat)
    
    print(clim_sens)
    #m.scatter(x1,y1,s=sizes,c=cols,marker="o",cmap=cm.cool,alpha=0.7)
    cs = m.scatter(x1,y1,s=60,c=clim_sens,marker="o",cmap='rainbow')
    cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    #cbar.set_label('climate sensitivity (degC)',labelpad=-40,size=15)
    cbar.set_label('climate sensitivity (deg C)')
    #plt.show()
    print('saving figure as',fileout)
    plt.savefig(fileout)
    plt.close()
    
    txtfile1 = open(TEXTFILE,"w+")
   
    txtfile1.write('lon, lat, est_ECS \n')
    for i, lon in enumerate(x1):
        txtfile1.write((np.str(np.around(lon,2)) + 
                       ',' + np.str(np.around(y1[i],2)) + 
                       ',' + np.str(np.around(clim_sens[i],2)) +  '\n'))
        
    txtfile1.write('\n')
  
def plotdata_overplot(lat,lon,clim_sens,nproxies,pcube, signcube, fileout):
    """
    this will plot on a single map all the gridpoints where there is 
    a significant relationship between pliocene temp anomaly and climate sensitivity
    on top of this it will plot the climate sensitivities obtained from the data
    """
    
  
    m=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,
              urcrnrlat=90.0,projection='cyl',resolution='c')
    m.drawmapboundary
    m.drawcoastlines()
    parallels=np.arange(-90.,90.,50.)
    m.drawparallels(parallels,labels=[True,True,True,True],fontsize=10) # labels right
    meridians=np.arange(-180.,180.,60.)
    m.drawmeridians(meridians,labels=[True,True,True,True],fontsize=10)
    
    iplt.contourf(signcube, 1, colors=[[0.8, 0.8, 0.8], [1, 1, 1]])
    #iplt.contourf(signcube, 1, hatches=[None, '///'], colors='none')
    #iplt.contourf(signcube, 1, hatches=[None, '\\\''], colors='none')
   
    x1,y1=m(lon,lat)
    
    print(clim_sens)
    #m.scatter(x1,y1,s=sizes,c=cols,marker="o",cmap=cm.cool,alpha=0.7)
    #cs = m.scatter(x1,y1,s=65,marker="o","black")
    cs = m.scatter(x1, y1, s=60, c=clim_sens, marker="o",
                   cmap='rainbow', edgecolors='black')
    cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    #cbar.set_label('climate sensitivity (degC)',labelpad=-40,size=15)
    cbar.set_label('climate sensitivity (deg C)')
    m.drawmapboundary
    m.drawcoastlines()
    parallels=np.arange(-90.,90.,50.)
    m.drawparallels(parallels,labels=[True,True,True,True],fontsize=10) # labels right
    meridians=np.arange(-180.,180.,60.)
    m.drawmeridians(meridians,labels=[True,True,True,True],fontsize=10)
    plt.savefig(fileout)
    plt.close()
    
def redu_sites(lats, lons, proxysst, clim_sens, 
               proxysst_full, proxy_fulllat, proxy_fulllon):
    """
    # remove sites which we are not sure abou
    # this is where the datapoint is not within 1 deg of the modelled range
    # ie the data does not even nearly agree with any of the models
    # ie where data and model do not agree at all   
    
    input latitude longitude and original climate sensitivity
    output new latitude longitude climate sensitivitiy and nsites
    """
    print('in redu sites')
    (mod_lat, mod_lon, mod_minsst, mod_maxsst, 
    mod_allsst, modnames) =  readmodel()
    print_rmse(mod_allsst, proxysst_full, mod_lat, mod_lon, proxy_fulllat, 
               proxy_fulllon, modnames)
                                               
    nvals = 0
    new_latredu = []
    new_lonredu = []
    new_clim_sens_redu = []

    for i, lat_i in enumerate(lats):
        lon_i = np.around(lons[i],2)
        j = np.where(mod_lat == lat_i)
        
        if lon_i == mod_lon[j]:
            print('here', lat_i, lon_i, proxysst[i], mod_minsst[j], mod_maxsst[j])
            if mod_minsst[j] - 1.0 < proxysst[i] < np.min([mod_maxsst[j] + 1.0, 9.0]) :
                new_latredu.append(lat_i)
                new_lonredu.append(lon_i)
                new_clim_sens_redu.append(clim_sens[i])
                nvals = nvals + 1
         
        else:
            print('j is', j)
            print('latlon mismatch',i, lat_i, mod_lat[j],lon_i,mod_lon[j])
            sys.exit(0)
   
            
    return (np.asarray(new_latredu, dtype=float), 
            np.asarray(new_lonredu,dtype=float), 
            nvals, np.asarray(new_clim_sens_redu, dtype = float)) 
                     
    
##############################################################
LINUX_WIN = 'l'
if LINUX_WIN == 'l':
    FILECS = ('/nfs/hera1/earjcti/regridded/allplots/NearSurfaceTemperature/climate_sensitivity_relationships.txt')
#    FILECS = ('/nfs/hera1/earjcti/regridded/allplots/NearSurfaceTemperature/climate_sensitivity_relationships.txt')
    FILEPROXY = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/cs_mp_sst_data_30k_plusNOAA.csv')
    FILEOUT = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/climate_sensitivity.pdf')
    FILEOUT_R = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/climate_sensitivity_redu.pdf')
    FILEOUT_S = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/climate_sensitivity_redu_significant.pdf')
    FILEMODEL = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/modeloutput_CSCD_nearsites.csv')
    SIGNIFICANCE_FILE = ('/nfs/hera1/earjcti/regridded/alldata/data_for_7d.nc')
    SIGNIFICANCE_NAME = 'pvalue'
    TEXTFILE = '/nfs/hera1/earjcti/regridded/alldata/figure9.txt'
else:
    FILECS = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\allplots\\NearSurfaceTemperature\\climate_sensitivity_relationships.txt'
    FILEPROXY = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\proxydata\\cs_mp_sst_data_30k_plusNOAA.csv')
    FILEMODEL = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\proxydata\\modeloutput_CSCD_nearsites.csv')
    FILEOUT = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\allplots\\NearSurfaceTemperature\\climate_sensitivity.pdf')
    FILEOUT_R = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\allplots\\NearSurfaceTemperature\\climate_sensitivity_redu.pdf')
    FILEOUT_S = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\allplots\\NearSurfaceTemperature\\climate_sensitivity_redu_significant.pdf')
    SIGNIFICANCE_FILE = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\alldata\\data_for_7d.nc')
    SIGNIFICANCE_NAME = 'pvalue'
    TEXTFILE = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\alldata/figure9.txt'

main()
