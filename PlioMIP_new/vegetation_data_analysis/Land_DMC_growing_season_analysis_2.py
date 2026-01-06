#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on September 2020
# 
#  This program produces additional of figures to assess
#  the DMC in various ways.  
#
#  The ultimate aim is to do a better DMC by comparing the model
#  to the times when the proxy represents.  
#
#
#
#import os
import numpy as np
import pandas as pd
#import scipy as sp
#import cf
import iris
#import iris.util
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import netCDF4
#from mpl_toolkits.basemap import Basemap, shiftgrid
#from netCDF4 import Dataset, MFDataset
#import iris.analysis.cartography
#import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
#import cf_units as unit
#from iris.experimental.equalise_cubes import equalise_attributes
import cartopy.crs as ccrs
#import matplotlib.ticker as mticker
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from mpl_toolkits.basemap import Basemap

import sys

class GetMeteodata:
    """
    This class is everything to do with 
    getting the meteorological data.

    """
    def __init__(self):
        
        # sel.nearsites = [sitename, sitelat, sitelon]
        if SITE_REQUIRED == 'Lake Baikal':
            self.nearsites = [['Krasnojarsk', 56.0, 92.3],
                              ['Nizneudinsk', 54.9, 99.0],
                              ['Zigalovo', 54.8, 105.2],
                              ['Nizneangarsk', 55.8, 109.6],
                              ['Kalakan', 55.1, 116.8]
                              ]
            self.nearest = 'Nizneangarsk'
            
        if SITE_REQUIRED == 'Lake Baikal 19501970':
            self.nearsites = [['Zigalovo_1950_1970', 54.8, 105.2],
                              ['Nizneangarsk_1950_1970', 55.8, 109.6],
                              ['Kalakan_1950_1970', 55.1, 116.8]
                              ]
            self.nearest = 'Nizneangarsk_1950_1970'
            
    def get_meteoanncycle(self):
        
        anncycle = {
            "Krasnojarsk" : [-16.7575, -14.9292, -7.16471, 1.79832,    
                             9.57983, 16.3615, 19.0419, 15.8119,
                             9.34746, 1.61453, -8.17542, -15.3414], 
            "Nizneudinsk": [-21.3235, -19.0240, -9.20600, 0.672549,
                            7.95306, 14.5939, 17.0184, 14.1551,
                            7.68367, -0.253061, -10.9122, -19.0490],  
            "Nizneangarsk": [-22.2907, -21.3667,  -12.8852, -3.03148,
                             4.41731, 11.7189, 15.5906, 14.5269, 
                             7.85741, -1.11481, -11.3630, -17.8389],
            "Kalakan" : [-34.8708, -27.9969, -16.5879, -3.28750, 
                         6.24923, 13.7091, 16.6369, 13.5492,  
                         5.68438,  -5.90909, -22.3375, -32.9889],
            "Zigalovo" : [-27.9258, -24.3152, -12.6881, -0.409091,  
                          7.56562, 14.8478,  17.5433, 14.2940, 
                          6.71912, -2.14118, -15.0515, -24.7530],
            "Kalakan_1950_1970" : [-35.7381, -29.6000,   -17.4048,
                                   -4.49500,   5.98571,   13.4190, 
                                   16.4429,    13.5048,  5.23333,
                                    -6.30952,   -23.4200, -33.8190],     
            "Nizneangarsk_1950_1970" :[-22.3842,  -21.9474,  -13.5684,
                                       -3.92632,   3.93684,  10.6421,
                                        14.9789, 14.2000,  7.21000,
                                        -0.935000,  -12.5700,  -18.4150],
            "Zigalovo_1950_1970" : [-28.3476,  -25.6952,  -14.0571, 
                                    -1.68571,  7.35714,  14.5143,
                                    17.5952,  14.4571,  6.22857,
                                     -2.30476, -16.6238,  -25.5000]
                }
        
        monthnames = ['ja','fb','mr','ar','my','jn',
                       'jl','ag','sp','ot','nv','dc']
        allseascyc = []
        allseascyc_anom = []
         
        fig = plt.figure(figsize=(11.0, 8.5))
        ax=plt.subplot(2,1,1)
         
        for siteinfo in self.nearsites:
            sitename = siteinfo[0]
            sitelon = np.str(np.round(siteinfo[2])) + ' deg E'
            site_anncycle = anncycle.get(sitename)
            print(sitename, np.mean(site_anncycle))
            allseascyc.append(site_anncycle)
            if sitename == self.nearest:
                ann_cycle_nearest = site_anncycle
            ax.plot(monthnames, site_anncycle,label=sitename + ' ' + sitelon)
            if (sitename[0:7]) == 'Zigalov':
                ann_cycle_zigalov = np.array(site_anncycle)
                lon_zigalov = siteinfo[2]
            if (sitename[0:7]) == 'Kalakan':
                ann_cycle_kalakan = np.array(site_anncycle)
                lon_kalakan = siteinfo[2]
                anom = ((ann_cycle_zigalov - ann_cycle_kalakan) 
                                         * (109.0 - lon_kalakan)
                                         / (lon_zigalov - lon_kalakan))
                est_ann_cycle_Baikal = ann_cycle_kalakan + anom
                print('est cycle',est_ann_cycle_Baikal)
                print('est mean baikal',np.mean(est_ann_cycle_Baikal))
               
        print('anomalous annual cycle', est_ann_cycle_Baikal - ann_cycle_nearest)   
    
      
       
        box = ax.get_position()
        ax.set_position([box.x0, box.y0+(0.3*box.height), box.width * 0.8, (0.8*box.height)])
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title('annual cycle for sites near ' + SITE_REQUIRED)
        
        ax=plt.subplot(2,1,2)
        for siteinfo in self.nearsites:
            sitename = siteinfo[0]
            sitelon = np.str(np.round(siteinfo[2])) + ' deg E'

            anncyc_anom = np.array(anncycle.get(sitename)) - np.array(ann_cycle_nearest)
            allseascyc_anom.append(anncyc_anom)
            if sitename != self.nearest:
                ax.plot(monthnames, anncyc_anom, label= sitename + ' ' + sitelon)
              
        ax.axhline(y=0)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0+(0.3*box.height), box.width * 0.8, (0.8*box.height)])
        ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.title('annual cycle anomaly from ' + self.nearest)
        
        
        plt.savefig(FILESTART + '/Growing_seas/Meteorological_data/seas_cyc_near_' + SITE_REQUIRED + PLOTTYPE)
        plt.close()
        
        return (self.nearsites, allseascyc, allseascyc_anom,
                est_ann_cycle_Baikal)
        





def get_model_data(period, modelname, lat, lon):
    """
    gets the model data for each month of the year for the site.
    """
    
    lonreq=lon
    if lonreq < 0.:
        lonreq=lonreq + 360.
   
    filename = (FILESTART + '/regridded/' + modelname + '/'+ period + '.'
                + FIELD + '.mean_month.nc')
    
    
    fieldname = FIELD
   
    if modelname == 'GISS2.1G' or modelname == 'IPSLCM6A':
        if fieldname == 'NearSurfaceTemperature':
            fieldname = 'air_temperature'
        if fieldname == 'TotalPrecipitation':
            fieldname = 'precipitation_flux'
    
  
    field_cube = iris.load_cube(filename,
                              fieldname)
    
    
    lat_ix = (np.abs(field_cube.coord('latitude').points - lat)).argmin()
    lon_ix = (np.abs(field_cube.coord('longitude').points - lonreq)).argmin()
    
    seas_field = field_cube.data[:, lat_ix, lon_ix]
   
   
    return seas_field.data



def plot_seascyc_vs_data(proxyT, models_Tseas, siteinfo, period):
    """
    plots the temperature from the data vs the modelled seasonal cycle.

    """
    
    ylab = {'NearSurfaceTemperature': 'degC',
            'TotalPrecipitation': 'mm/day'}
    
    ax = plt.subplot(111)
    
    print(siteinfo[0])
    for i, mod_data in enumerate(models_Tseas):
        if i > 6:
            ax.plot(mod_data, linestyle='dashed', label=MODELNAMES[i])
        else:
            ax.plot(mod_data, label=MODELNAMES[i])
        MAT = np.mean(mod_data)
        ax.plot([-1.0,0], [MAT, MAT], color='blue')
        print(MODELNAMES[i],
              np.around(np.sqrt(((mod_data - proxyT) ** 2).mean())))
        
        
       # print('mean=',np.around(np.mean(mod_data),2), 
       #       ' median=', np.around(np.median(mod_data),2),
       #       ' avg wam/cold=', np.around((np.max(mod_data) + np.min(mod_data)) / 2.0, 2))
        
    
    ax.plot(proxyT, label=siteinfo[0], color='black', linewidth=2)
    title = (siteinfo[0] + ' lat=' + np.str(np.around(siteinfo[1])), 
        ' lon='+np.str(np.around(siteinfo[2])))
    plt.title(title)
    plt.xlabel('month')
    plt.ylabel(ylab.get(FIELD))
    
    
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    
    plt.savefig(FILESTART + '/Growing_seas/Meteorological_data/seas_cyc_' + siteinfo[0] +  '_' + period + PLOTTYPE)
    plt.close()


def plot_model_near_site(lons, all_month_pi, all_month_plio, monthname):
    """
    plots the modelled temperature at the site and for some nearby
    gridboxes
    lons = list of longitudes
    all_seascyc = seasonal cycle at (list of longs (list of models))
    """
    
    fig = plt.figure(figsize=(11.0, 11.0))
    fig.suptitle(monthname + 'temperature gradient: pi (blue), plio (red)', fontsize=20)
    fig.text(0, 0.5, 'temperature deg C', rotation=90, fontsize=20, verticalalignment = 'center')
    fig.text(0.6, 0.0, 'longitude deg E', fontsize=20, horizontalalignment = 'center')
    for i, model in enumerate(MODELNAMES):
        ax = plt.subplot(4,4,i+1)
        plt.title(model)
        monthtemp = all_month_pi[:, i] - all_month_pi[0, i]
        monthtemp_plio = all_month_plio[:, i] - all_month_plio[0, i]
        ax.plot(lons, monthtemp, color='blue')
        ax.plot(lons, monthtemp_plio, color='red')
    
       
    plt.subplots_adjust(left=0.05, bottom=0.05, right=1.00, top=0.9,
                        wspace=0.3, hspace=0.3)
       
    plt.savefig(FILESTART + '/Growing_seas/' + monthname + '_near_' + SITE_REQUIRED + PLOTTYPE, bbox_inches='tight')
    plt.close()
        




  
def main():
    """
    calling structure
    a) get modern meteorological data seasonal cycle
    b) get modern data at the location of interest
    """

#    #########################################
#    # get data for modern sites.  
    
    siteinfo = GetMeteodata() # get data for t1 timeslice
    (sitenames, seascyc, seascyc_anom,
     est_ann_cycle_Baikal_nolake)  = siteinfo.get_meteoanncycle()
       
    sitenames.append(['Estimated Lake Baikal (without lake)', 
                      55.0, 108.0])
    seascyc.append(est_ann_cycle_Baikal_nolake)
    
    #################################################
    # get PI model data for each month of the year.
    for siteno, sitedetails in enumerate(sitenames):
        sitelat = sitedetails[1]
        sitelon = sitedetails[2]
        allmodel_TseascycPI = []
   
        for i, model in enumerate(MODELNAMES):
            T_seascyc = get_model_data('E280', model, sitelat, sitelon)
            allmodel_TseascycPI.append(T_seascyc)
      
        plot_seascyc_vs_data(seascyc[siteno], allmodel_TseascycPI,
                          sitedetails, 'E280')
        
    #############################################
    # look and see if the modelled seasonal cycle
    # at the location is different from nearby gridboxes
        
    lons_list = []
    PI_Jan = []
    Plio_Jan = []
    PI_Jul = []
    Plio_Jul = []
   
    for lon in range(np.int(np.floor(sitelon)) - 5, 
                     np.int(np.ceil(sitelon)) + 5):
        lons_list.append(lon)
        allmodel_Jan = []
        allmodel_Jan_plio = []
        allmodel_Jul = []
        allmodel_Jul_plio = []
        for model in MODELNAMES:
            T_seascyc = get_model_data('E280', model, sitelat, lon)
            T_seascyc_plio = get_model_data('EOI400', model, sitelat, lon)
            allmodel_Jan.append(T_seascyc[0])
            allmodel_Jan_plio.append(T_seascyc_plio[0])
            allmodel_Jul.append(T_seascyc[6])
            allmodel_Jul_plio.append(T_seascyc_plio[6])
            
            
        PI_Jan.append(allmodel_Jan)
        Plio_Jan.append(allmodel_Jan_plio)
        PI_Jul.append(allmodel_Jul)
        Plio_Jul.append(allmodel_Jul_plio)
        
    plot_model_near_site(lons_list, np.asarray(PI_Jan),
                         np.asarray(Plio_Jan),
                         'January')
    plot_model_near_site(lons_list, np.asarray(PI_Jul),
                         np.asarray(Plio_Jul),
                         'July')
            
  
    
 
       
##########################################################
# main program


LINUX_WIN = 'w'


if LINUX_WIN == 'l':
    FILESTART = '/nfs/hera1/earjcti/'   
    LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')
    PLOTTYPE = '.eps'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
    LAND_DATAFILE = (FILESTART + '/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')
    PLOTTYPE = '.png'
 

MODELNAMES = ['CESM2', 
              'IPSLCM6A', 
              'COSMOS', 
              'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
              'MIROC4m', 'IPSLCM5A2', 'HadCM3',
              'GISS2.1G', 'CCSM4', 
              'CCSM4-Utr', 'CCSM4-UoT', 
              'NorESM-L', 'MRI2.3', 'NorESM1-F'
              ]

#MODELNAMES = ['CESM2','IPSLCM6A']

ABSOLUTE_ANOMALY_IND = 'absolute' # do we want to compare the absolute value or the anomaly value.
FIELD = 'NearSurfaceTemperature'
#FIELD = 'TotalPrecipitation'

SITE_REQUIRED = 'Lake Baikal'
#SITE_REQUIRED = 'Lake Baikal 19501970'
#SITE_REQUIRED = 'Lake Elgygytgyn'
#SITE_REQUIRED = 'Meighen Island'
#SITE_REQUIRED = 'Lost Chicken Mine'
#SITE_REQUIRED = 'James Bay Lowland'
#SITE_REQUIRED = 'Pula Maar'
#SITE_REQUIRED = 'Alpes-Maritimes'
#SITE_REQUIRED = 'Tarragona'
#SITE_REQUIRED = 'Rio Maior'
#SITE_REQUIRED = 'Yallalie, Perth'

main()
