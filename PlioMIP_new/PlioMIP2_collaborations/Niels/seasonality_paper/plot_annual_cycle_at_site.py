#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#Created on Fri Jul  5 15:11:26 2019
#Updated for JH on 17th May 2021
#
#@author: earjcti
#"""
#
#   This program will obtain the SST data from the pliovar site locations and process
#
#
# This program has been ammended from 
#PlioMIP2/large_scale_features/extract_data_locations.py
#
#

import pandas as pd
import csv
import sys
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import iris
import xlwt
from xlwt import Workbook
import os
import matplotlib.cm as cm
from matplotlib.colors import Normalize



############################################################################
class Getinitialdata:
# get all of the initial data, including filenames and the lons and lats where
# we require model output
    def __init__(self, linuxwin_, datafile_):
        self.linuxwin = linuxwin_
        if self.linuxwin == 'l':
            filename = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/'
        else:
            filename = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\proxydata\\'

        if datafile_ == 'Jonathan':
            self.filenamein = filename + 'metadata_JH_based_on_pliovar.csv'
            self.latcolumn = 3
            self.loncolumn = 2
            self.sitecolumn = 1
            self.outend = 'modeloutput_JH.xls'
            self.sitesreq = ['ODP662','ODP999','DSDP606','DSDP607','U1313',
                             'U1308','DSDP552','ODP982','ODP642']
        if datafile_ == 'Other':
            self.filenamein = filename + 'one_locality.csv'
            self.latcolumn = 1
            self.loncolumn = 2
            self.sitecolumn = 0
            self.outend = 'test_localities.xls'

        if self.linuxwin == 'l':
            self.filestart = '/nfs/hera1/earjcti/regridded/'
            if pliomip1 == 'y':
                self.fileout = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/pliomip1_'+self.outend
            else:
                self.fileout = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/'+self.outend
        else:
            self.filestart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
            self.fileout = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\proxydata\\'+self.outend
        


    def read_file(self):
        count = 0
        latlist = []
        lonlist = []
        lonlist_alt = []
        sitelist = self.sitesreq

        df = pd.read_csv(self.filenamein, delimiter=',')
        for site in sitelist:
            print('site is',site)
            dfrow = df.loc[lambda df: df['name'] == site]
            print(dfrow)
            lat = dfrow.iloc[:, self.latcolumn]
            lon = dfrow.iloc[:, self.loncolumn]
            latlist.append(pd.DataFrame.to_numpy(lat)[0])
            lonarr = pd.DataFrame.to_numpy(lon)[0]
            lonlist_alt.append(lonarr)
            if lonarr < 0.:# longitude goes from 0-360 in models
                lonlist.append(lonarr+360.)
            else:
                lonlist.append(lonarr)
            lonlist.append
           

        for i, site in enumerate(sitelist):
            print(site, lonlist[i], lonlist_alt[i], latlist[i])

        returndata = self.filestart,self.fileout,lonlist, latlist, lonlist_alt, sitelist
        return returndata

# end of class getinitdata
###############################################################################

class Getmodeldata:
    # get all of the data from the model at the required gridpoints
    def __init__(self, test, modelstart_, field, latlist, lonlist, period):

        fieldunits = {
            "SST" : "degC",
            "TotalPrecipitation" : "mm/day"
                        }


        self.fieldnames = field
        self.latlist = latlist
        self.lonlist = lonlist
        self.modelstart = modelstart_ # the start of the filename for the model
        self.period = period

       
        if period == 'E280':
            self.modelnames.append('HadISST')
            self.modelnames.append('NOAAERSST5')

        self.units = fieldunits.get(fieldnames)


        #
    def extract_model_points(self, filenameuse):
        """
        will extract the data at each point from 'filenameuse'

        calls: get_near_data
        """

        cube = iris.load_cube(filenameuse)

        cubelats = cube.coord('latitude').points
        cubelons = cube.coord('longitude').points

        model_data = np.zeros(len(self.lonlist))
        model_data_near = np.zeros(len(self.lonlist)) # values near the point
        near_distance = np.zeros(len(self.lonlist)) # how far away we have to look to get data
        ngbox_avg = np.zeros(len(self.lonlist)) # how many gridboxes we are averaging over to get data

        for i in range(0,len(self.lonlist)):
            # find nearest latitude and lontiude to the value
            latix = (np.abs(cubelats-self.latlist[i])).argmin()
            lonix = (np.abs(cubelons-self.lonlist[i])).argmin()

            print(self.latlist[i], self.lonlist[i])
            # get data from this location
            data_slice  =  cube.extract(iris.Constraint(
                        latitude = cubelats[latix],longitude = cubelons[lonix]))
            
            data_slice_data = data_slice.data
            if -100. < data_slice_data < 100.: 
                model_data[i] = data_slice_data
            else:
                model_data[i] = float('NaN')
        
            model_data_near[i] = model_data[i]

            count_near_gb = 0 # how many gridboxes away are we looking for data
            ngboxes = 1 # number of gridboxes we are averaging over when looking at 'near points'

            # while value is unknown gradually expand the region to look for near gridboxes
            while np.isnan(model_data_near[i]):
                # get nearest neighbours within 'count_near_gb' gridboxes
                count_near_gb = count_near_gb+1
                print(count_near_gb)
                neardata,ngboxes = self.get_near_data(cube,lonix,latix,cubelons,cubelats,count_near_gb)
                model_data_near[i] = neardata

            near_distance[i] = count_near_gb # how far away are we looking for data
            ngbox_avg[i] = ngboxes


        returndata = [model_data,model_data_near,near_distance,ngbox_avg]
        return returndata



    def get_near_data(self, cube, lonix, latix, cubelons, cubelats, npt):
    # if there is no data at the given gridpoint get the data near the gridpoint

        count_finite = 0
        count_nan = 0
        totdata = 0.
        nlons = len(cubelons)
        for i2 in range(lonix-npt,lonix+npt+1):
            i3 = i2
            if i2 >=  nlons:
                i3 = i2-nlons
            for j2 in range(latix-npt,latix+npt+1):
                data_slice_new = cube.extract(iris.Constraint(
                     latitude = cubelats[j2],longitude = cubelons[i3]))
                data2 = data_slice_new.data
                if np.ma.is_masked(data2):
                    count_nan = count_nan+1
                else:
                    count_finite = count_finite+1
                    totdata = totdata+data2
        if count_finite > 0:
            data_near = totdata/count_finite
        else:
            data_near = float('NaN') # if no data near set to nan

        return data_near,count_finite


    def extract_all(self):
        """
        extract points from all models for timeperiod.
        timeperiod is likely to be 'E280' or 'EOI400'

        returns
        modelnames (strarr) modelnames used for this period
        sitevals (np.arr): the values at the sites
        sitenear (np.arr): an average of the values nearest to the sites
        sitenear_dist (np.arr): how far away the values presented are
        sitenear_ngbox_avg) (np.arr): the number of gridboxes averaged where the
                      values near to the sites are used

        """

        npoints = len(self.lonlist)
        nmodels = len(self.modelnames)
        sitenear = np.zeros((nmodels, npoints)) # data near point
        sitevals = np.zeros((nmodels, npoints)) # data at point
        sitenear_dist = np.zeros((nmodels, npoints)) # how far away we have to look
        sitenear_ngbox_avg = np.zeros((nmodels, npoints)) # how many gridboxes we are averaging over


        for model in range(0, len(self.modelnames)):
            print(self.modelnames[model])
            filename = (self.modelstart + self.modelnames[model] + '/' +
                               self.period + '.SST.allmean.nc')

            # get model points and how far away they are from data
            (sitevals[model,:],sitenear[model,:],sitenear_dist[model,:],
                sitenear_ngbox_avg[model,:]) = self.extract_model_points(filename)

        return [self.modelnames,sitevals,
                sitenear,sitenear_dist,sitenear_ngbox_avg]

# end of class Getmodeldata

###############################################################################

class Getmodeldata_p1:
    # get all of the data from the pliomip1 models at the required gridpoints
    def __init__(self, test, modelstart_, field, latlist, lonlist, period):

        fieldunits = {
            "SST" : "degC",
            "TotalPrecipitation" : "mm/day"
                        }


        self.fieldnames = field
        self.latlist = latlist
        self.lonlist = lonlist
        self.modelstart = modelstart_ # the start of the filename for the model
        self.period = period

        if test == 'y':
            self.modelnames = ['NOR']
        else:
            self.modelnames = ['COSMOS', 'Had', 'CCSM',
                               'IPSL', 'MIROC', 'MRI', 'NOR']

        if period == 'E280':
            self.modelnames.append('HadISST')
            self.modelnames.append('NOAAERSST5')

        self.units = fieldunits.get(fieldnames)


    def extract_obs(self, filenameuse):
        """
        will extract the data at each point from NOAAERSST5 or HadISST
        """

        cube = iris.load_cube(filenameuse)

        cubelats = cube.coord('latitude').points
        cubelons = cube.coord('longitude').points

        model_data = np.zeros(len(self.lonlist))
        for i in range(0,len(self.lonlist)):
            # find nearest latitude and lontiude to the value
            latix = (np.abs(cubelats-self.latlist[i])).argmin()
            lonix = (np.abs(cubelons-self.lonlist[i])).argmin()

            # get data from this location
            data_slice  =  cube.extract(iris.Constraint(
                        latitude = cubelats[latix],longitude = cubelons[lonix]))
            
            data_slice_data = data_slice.data
            if -100. < data_slice_data < 100.: 
                model_data[i] = data_slice_data
            else:
                model_data[i] = float('NaN')
        
         
        return model_data

    def extract_model_points_p1(self, model):
        """
        will extract the data at each point from 'filenameuse'

        calls: get_near_data
        """
        filename = '/nfs/hera1/earjcti/PLIOMIP/PlioMIP1_regridded.nc'
        print(self.period, model)
        if self.period == 'E280':
            field = model + '_ctrl_sst'
        if self.period == 'EOI400':
            field = model + '_plio_sst'

        cubeall = iris.load(filename)
        for cube_temp in cubeall:
            var = cube_temp.var_name
            if field.lower() in var.lower():
               cube=cube_temp

        print(field, cube)
       
        cubelats = cube.coord('latitude').points
        cubelons = cube.coord('longitude').points

        model_data = np.zeros(len(self.lonlist))
       
        for i in range(0,len(self.lonlist)):
            # find nearest latitude and lontiude to the value
            if self.lonlist[i] > 180.:
                lonreq = self.lonlist[i] - 360.
            else:
                lonreq = self.lonlist[i]
            
            latix = (np.abs(cubelats-self.latlist[i])).argmin()
            lonix = (np.abs(cubelons-lonreq)).argmin()
            # get data from this location
            data_slice  =  cube.extract(iris.Constraint(
                        latitude = cubelats[latix],longitude = cubelons[lonix]))
            
            if model == 'NOR':
                data_slice_data = data_slice.data - 273.15
            else:
                data_slice_data = data_slice.data

            if -100. < data_slice_data < 100.: 
                model_data[i] = data_slice_data
            else:
                model_data[i] = float('NaN')
       
        return model_data



   

    def extract_all_p1(self):
        """
        extract points from all models for timeperiod.
        timeperiod is likely to be 'E280' or 'EOI400'

        returns
        modelnames (strarr) modelnames used for this period
        sitevals (np.arr): the values at the sites
       
        """

        npoints = len(self.lonlist)
        nmodels = len(self.modelnames)
        sitevals = np.zeros((nmodels, npoints)) # data at point

        for modno, modelname in enumerate(self.modelnames):
           if modelname == 'HadISST' or modelname == 'NOAAERSST5':
              filename = (self.modelstart + modelname + '/' +
                          self.period + '.SST.allmean.nc')

              sitevals[modno,:] = self.extract_obs(filename)
           else:
               # get model data
               sitevals[modno,:] = self.extract_model_points_p1(modelname)

        return [self.modelnames,sitevals]
# end of class Getmodeldata

###############################################################################
def plotpoints(lonlist,latlist,datalist):
# plot the points we have got from the file



    fig,ax = plt.subplots()
    alllons = np.arange(-180,180,1)
    alllats = np.arange(-90,90,1)
    lons,lats = np.meshgrid(alllons,alllats)
    map = Basemap(llcrnrlon = -180.0,urcrnrlon = 180.0,llcrnrlat = -90.0,
                urcrnrlat = 90.0,projection = 'cyl',resolution = 'c')
    map.drawmapboundary
    x,y = map(lons,lats)
    map.drawcoastlines()

    valmin = np.nanmin(datalist)
    valmax = np.nanmax(datalist)

    norm  =  mpl.colors.Normalize(vmin = valmin, vmax = valmax)
    cmap  =  cm.brg


    xpts,ypts = map(lonlist,latlist)
    incr = (valmax-valmin+1.0)/10.
    V = np.arange(valmin,valmax,incr)
    cvals = (datalist-valmin)/(valmax-valmin) # scale cval onto same scale as colorbar
    coluse = cmap(cvals)
    cs  =  map.scatter(xpts,ypts,color = coluse,marker = 'o')

    sm  =  plt.cm.ScalarMappable(cmap = cmap, norm = norm)
    sm.set_array([])
    plt.colorbar(sm, ticks = V,#ticks = np.linspace(valmin,valmax,incr),
             orientation = "horizontal",extend = "both")

    plt.show()

def model_correct():
    # calculate eoi400 - e280 + noaa_ersstv5 and eoi400-e280 + hadisst
    noaaix = modelnames_e280.index('NOAAERSST5')
    hadix = modelnames_e280.index('HadISST')

    ny, nx = eoi400.shape
    eoi400_corr_hadiss = eoi400 - e280[0:ny,:] + e280[hadix, :]
    eoi400_corr_noaa = eoi400 - e280[0:ny,:] + e280[noaaix, :]
   
    print(eoi400_corr_noaa.shape, eoi400.shape)
    return eoi400_corr_hadiss, eoi400_corr_noaa

def write_sheet(wb, style, sheetname, modelnames, lonlist_alt, latlist, datawrite, sitename):
    sheet  =  wb.add_sheet(sheetname)
    sheet.write(0,0,'site')
    sheet.write(0,1,'lat')
    sheet.write(0,2,'lon')
    for model in range(0,nmodels):
        sheet.write(0,3+model,modelnames[model])

    for i in range(0,npoints):
        sheet.write(i+1,0,sitename[i],style)
        sheet.write(i+1,1,latlist[i],style)
        sheet.write(i+1,2,lonlist_alt[i],style)
        for model in range(0,nmodels):
            sheet.write(i+1,3+model,datawrite[model,i],style)

    sheet.write(0,3+nmodels,'MMM')
    sheet.write(0,4+nmodels,'MM-SD')
    for i in range(0,npoints):
        sheet.write(i+1,3+nmodels,np.nanmean(datawrite[0:nmodels,i]),style)
        sheet.write(i+1,4+nmodels,np.nanstd(datawrite[0:nmodels,i]),style)

    # add extra columns if we have them  this is likely to be hadisst
    if len(modelnames) > nmodels:
        print(modelnames,nmodels)
        for model in range(nmodels,len(modelnames)):
            sheet.write(0,5+model,modelnames[model])
            for i in range(0,npoints):
                sheet.write(i+1,5+model,datawrite[model,i],style)



######################################################################################
def write_to_book(fileout,lonlist,latlist,lonlist_alt, sitename):
    # write to workbook
    # calls write_sheet

    # Workbook is created
    wb  =  Workbook()

    style  =  xlwt.XFStyle()
    style.num_format_str  =  '0.00'


    # add_sheet for Eoi400 E280 and difference
    write_sheet(wb,style, 'EOI400_corr_NOAA', modelnames_eoi400, lonlist_alt, latlist, eoi400_corr_noaa, sitename)
    write_sheet(wb,style, 'EOI400_RAW', modelnames_eoi400, lonlist_alt, latlist, eoi400, sitename)
    write_sheet(wb,style, 'E280', modelnames_e280, lonlist_alt, latlist,e280, sitename)
    write_sheet(wb,style, 'EOI400-E280',modelnames_eoi400,lonlist_alt,latlist,
                eoi400[0:nmodels]-e280[0:nmodels], sitename)

    # add_sheet near Eoi400 E280 and difference
    # JULIA NOTE:  I HAVE GOT RID OF ALL THESE FOR JH BECAUSE ALL HIS POINTS
    # ARE DEFINATELY OCEAN.  
    #write_sheet(wb,style, 'EOI400near',modelnames_eoi400,lonlist_alt,latlist,eoi400_near, sitename)
    #write_sheet(wb,style, 'E280near',modelnames_e280,lonlist_alt,latlist,e280_near, sitename)
    #write_sheet(wb,style, 'EOI400-E280near',modelnames_eoi400,lonlist_alt,latlist,
    #            eoi400_near[0:nmodels]-e280_near[0:nmodels], sitename)

    # add sheet for how far away we need to look for data
    #write_sheet(wb,style, 'EOI400distance',modelnames_eoi400,lonlist_alt,latlist,eoi400_near_distance, sitename)
    #write_sheet(wb,style, 'E280distance',modelnames_e280,lonlist_alt,latlist,e280_near_distance, sitename)

    # add sheet for how many gridboxes we are averaging over
    #write_sheet(wb,style, 'EOI400nboxes',modelnames_eoi400,lonlist_alt,latlist,eoi400_ngbox_avg, sitename)
    #write_sheet(wb,style, 'E280nboxes',modelnames_e280,lonlist_alt,latlist,e280_ngbox_avg, sitename)




    # remove output file if it exists
    exists  =  os.path.isfile(fileout)
    if exists:
        os.remove(fileout)
    wb.save(fileout)

#################################################################################
#def plot_points():
    # plot all the points from eoi400_near[model,i]-3280_near[model,i] to a map
#    for model in range(0,len(modelnames)):


def plot_by_lat(SST, lat, outmid):
# do a temperature by latitude plot for all of the data
    print(modelnames_eoi400)
    print(SST.shape, np.mean(SST, axis=0))
    plt.scatter(lat, np.mean(SST, axis=0),label='multimodel mean')
    for i, model in enumerate(modelnames_eoi400):
        if i < 8:
            plt.scatter(lat, SST[i, :], s=5, marker='^', label=model)
        else:
            plt.scatter(lat, SST[i, :], s=5, marker='v', label=model)
   
    plt.errorbar(np.asarray(lat)-1.0, np.mean(SST, axis=0), yerr=np.std(SST, axis=0))
    plt.title(outmid)
    plt.xlabel('latitude')
    plt.ylabel('temp deg C')
    plt.legend(ncol=2, prop={'size':6})
    

    outstart = '/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/'
    outstart = outstart + 'Collaborators/Jonathan_Hall/'
    if pliomip1 == 'y':
        outstart + outstart + 'PlioMIP1_'
    plt.savefig(outstart + outmid + '.png')
    plt.savefig(outstart + outmid + '.eps')
    plt.close()


def plot_data(eoi400_data, e280_data):
    """
    eoi400/e280_data = np.zeros(len modelnames, nmonths)
    """

    # also write data to a text file
    textout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/Niels/'+ FIELD + '_seas_cyc_at_'+np.str(LAT_REQ)+'N_'+np.str(LON_REQ)+'E.txt')
    f1 = open(textout,'w+')
   
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:grey','tab:olive','tab:cyan','lightgray','indigo','coral','palegreen','salmon','chocolate','orchid']
    # absolute values
    ax1 = plt.subplot(2,1,1)
    mean400=np.zeros(12)
    meancount=0
    f1.write('absolute values \n')
    f1.write('modelname,')
    for month in months: f1.write(month + ', ')
    f1.write('\n')
    for i, model in enumerate(MODELNAMES):
        if (-500 < eoi400_data[i,0] < 1E10):
            ax1.plot(months, eoi400_data[i, :],color=colors[i], label=model,
                     linestyle='dashed')
            f1.write(model + ', ')
            for val in eoi400_data[i,:]: f1.write(np.str(np.round(val,2)) + ', ')
            f1.write('\n')
            mean400 = mean400 + eoi400_data[i, :]
            meancount = meancount+1
    mean400 = mean400 / meancount
    ax1.plot(months, mean400,color='black', label='MEAN')
          
    plt.title('mPWP temperature at '+np.str(LAT_REQ)+'N: '+np.str(LON_REQ)+'E')
    plt.legend(ncol=2,prop={'size':6})
    plt.ylabel('temp deg C')
  
    # anomaly
    f1.write('\nanomalies Plio-PI \n')
    f1.write('Modelname, ')
    for month in months: f1.write(month + ', ')
    f1.write('\n')
   
    ax2 = plt.subplot(2,1,2)
    meananom=np.zeros(12)
    meancount=0
    for i, model in enumerate(MODELNAMES):
        if ((-500 < eoi400_data[i,0]  < 1E10) 
            and (-500 < e280_data[i,0] < 1E10))  :
            ax2.plot(months, eoi400_data[i, :] - e280_data[i,:], 
                     color=colors[i],label=model, linestyle='dashed')
            # write to file
            f1.write(model + ', ')
            for val in eoi400_data[i,:] -  e280_data[i,:]: 
                f1.write(np.str(np.round(val,2)) + ', ')
            f1.write('\n')
            
            meananom = meananom + eoi400_data[i, :] - e280_data[i,:]
            meancount = meancount+1
    meananom = meananom / meancount
    ax2.plot(months, meananom,color='black', label='MEAN')
    print('anom is',meananom)
    f1.close()
  
    plt.title('mPWP - PI temperature at '+np.str(LAT_REQ)+'N: '+np.str(LON_REQ)+'E')
    plt.ylabel('temp deg C')
  
    plt.tight_layout()

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/Niels/'+ FIELD + '_seas_cyc_at_'+np.str(LAT_REQ)+'N_'+np.str(LON_REQ)+'E')
    plt.savefig(fileout + '.png')
    plt.savefig(fileout + '.eps')
    plt.close()
    sys.exit(0)

    
def extract_lon_lat(cube):
    """
    extracts the required longitude and latitude from the cube
    """ 

    cubelats = cube.coord('latitude').points
    cubelons = cube.coord('longitude').points

    # find nearest latitude and lontiude to the value
    latix = (np.abs(cubelats - LAT_REQ)).argmin()
    lonix = (np.abs(cubelons - LON_REQ)).argmin()

    data_slice  =  cube.extract(iris.Constraint(
                        latitude = cubelats[latix],longitude = cubelons[lonix]))
       

    return data_slice

def extract_data():
    """
    extracts the data at the required lat and long
    """
   
    eoi400_alldata = np.zeros((len(MODELNAMES),12))
    e280_alldata = np.zeros((len(MODELNAMES),12))

    for i, model in enumerate(MODELNAMES):
        filestart = '/nfs/hera1/earjcti/regridded/' + model + '/'
        eoi400_cube = iris.load_cube(filestart + 'EOI400.'+ FIELD + '.mean_month.nc')
        eoi400_cube_slice = extract_lon_lat(eoi400_cube)
        eoi400_alldata[i, :] = eoi400_cube_slice.data
    

        eoi400_cube = iris.load_cube(filestart + 'E280.'+ FIELD+ '.mean_month.nc')
        e280_cube_slice = extract_lon_lat(eoi400_cube)
        e280_alldata[i, :] = e280_cube_slice.data
    
    return eoi400_alldata, e280_alldata    
        


#################
# MAIN PROGRAM
################

###################################
# get initial data including the lats and longs we require

linuxwin = 'l'
LAT_REQ = 52.5 # 51.2N, 4.4E
LON_REQ = 3.0
test ='n'
#FIELD = 'NearSurfaceTemperature'
FIELD = 'SST'

if test == 'y':
     MODELNAMES = ['NorESM1-F', 'HadCM3']
else:
     MODELNAMES = ['CCSM4', 'CCSM4-UoT', 'CCSM4-Utr',  'CESM1.2','CESM2',
                   'COSMOS', 'EC-Earth3.3', 'GISS2.1G', 'HadCM3','HadGEM3',
                   'IPSLCM6A', 'IPSLCM5A2', 'IPSLCM5A', 'MIROC4m', 'MRI2.3',
                   'NorESM-L', 'NorESM1-F'
                               ]
eoi400_alldata, e280_alldata = extract_data()

plot_data(eoi400_alldata, e280_alldata)
