#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#Created on Fri Jul  5 15:11:26 2019
#
#@author: earjcti
#"""
#
#   This program will obtain the SST data from the pliovar site locations and process
#
#
import csv
import sys
#from mpl_toolkits.basemap import Basemap
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

        if datafile_ == 'Erin':
            self.filenamein = filename + 'pliovar_metadata_global_02102019.csv'
            self.latcolumn = 3
            self.loncolumn = 2
            self.sitecolumn = 1
            self.outend = 'modeloutput_pliovar.xls'
        if datafile_ == 'Harry':
            self.filenamein = filename + 'Copy_of_CSCD_localities.csv'
            self.latcolumn = 1
            self.loncolumn = 2
            self.sitecolumn = 0
            self.outend = 'modeloutput_CSCD_localities.xls'
        if datafile_ == 'Other':
            self.filenamein = filename + 'one_locality.csv'
            self.latcolumn = 1
            self.loncolumn = 2
            self.sitecolumn = 0
            self.outend = 'test_localities.xls'

        if self.linuxwin == 'l':
            self.filestart = '/nfs/hera1/earjcti/regridded/'
            self.fileout = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/'+self.outend
        else:
            self.filestart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
            self.fileout = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\proxydata\\'+self.outend

 
    def read_file(self):
        count = 0
        latlist = []
        lonlist = []
        lonlist_alt = []
        sitelist = []

        with open(self.filenamein) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            for row in readCSV:
                if count != 0: # not titleline
                    sitelist.append(row[self.sitecolumn])
                    latlist.append(np.float(row[self.latcolumn]))
                    lon = (np.float(row[self.loncolumn]))
                    lonlist_alt.append(lon)
                    if lon < 0.: # longitude goes from 0-360 in models
                        lonlist.append(lon+360.)
                    else:
                        lonlist.append(lon)
                count = count+1
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

        if test == 'y':
            self.modelnames = ['NorESM1-F', 'HadCM3']
        else:
            self.modelnames = ['CCSM4', 'CCSM4-UoT',
                               'CCSM4-Utr',  
                               'CESM1.2','CESM2',
                               'COSMOS', 'EC-Earth3.3', 
                               'GISS2.1G', 'HadCM3','HadCM3_new',
                               'IPSLCM6A', 'IPSLCM5A2', 'IPSLCM5A',
                               'MIROC4m', 'MRI2.3',
                               'NorESM-L', 'NorESM1-F'
                               ]

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

    #plt.show()



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
    write_sheet(wb,style, 'EOI400', modelnames_eoi400, lonlist_alt, latlist, eoi400, sitename)
    write_sheet(wb,style, 'E280', modelnames_e280, lonlist_alt, latlist,e280, sitename)
    write_sheet(wb,style, 'EOI400-E280',modelnames_eoi400,lonlist_alt,latlist,
                eoi400[0:nmodels]-e280[0:nmodels], sitename)

    # add_sheet near Eoi400 E280 and difference
    write_sheet(wb,style, 'EOI400near',modelnames_eoi400,lonlist_alt,latlist,eoi400_near, sitename)
    write_sheet(wb,style, 'E280near',modelnames_e280,lonlist_alt,latlist,e280_near, sitename)
    write_sheet(wb,style, 'EOI400-E280near',modelnames_eoi400,lonlist_alt,latlist,
                eoi400_near[0:nmodels]-e280_near[0:nmodels], sitename)

    # add sheet for how far away we need to look for data
    write_sheet(wb,style, 'EOI400distance',modelnames_eoi400,lonlist_alt,latlist,eoi400_near_distance, sitename)
    write_sheet(wb,style, 'E280distance',modelnames_e280,lonlist_alt,latlist,e280_near_distance, sitename)

    # add sheet for how many gridboxes we are averaging over
    write_sheet(wb,style, 'EOI400nboxes',modelnames_eoi400,lonlist_alt,latlist,eoi400_ngbox_avg, sitename)
    write_sheet(wb,style, 'E280nboxes',modelnames_e280,lonlist_alt,latlist,e280_ngbox_avg, sitename)




    # remove output file if it exists
    exists  =  os.path.isfile(fileout)
    if exists:
        os.remove(fileout)
    wb.save(fileout)

#################################################################################
#def plot_points():
    # plot all the points from eoi400_near[model,i]-3280_near[model,i] to a map
#    for model in range(0,len(modelnames)):


#################
# MAIN PROGRAM
################

###################################
# get initial data including the lats and longs we require

linuxwin = 'l'
datafile = 'Harry' # could have Harry or Erin or Other for test file
testdata = 'n'   # yes use one model no use full range of models
fieldnames = 'SST'

indata = Getinitialdata(linuxwin,datafile)
(modelstart, outputfile, longitudes, 
 latitudes, longitudes_alt, sitenames) = indata.read_file() # get the lats lons required and the number of sites
npoints = len(longitudes) # get the number of points


######################################
# setup a map and plot the points

#if linuxwin = ='l':
#    plotpoints(lonlist_alt,latlist,np.zeros(count))


##############################
# get the SST data from IRIS cubes

#eoi400
modeldata=Getmodeldata(testdata,modelstart,fieldnames,latitudes,longitudes,'EOI400')
(modelnames_eoi400,eoi400,eoi400_near,
 eoi400_near_distance,eoi400_ngbox_avg)=modeldata.extract_all()


modeldata=Getmodeldata(testdata,modelstart,fieldnames,latitudes,longitudes,'E280')
(modelnames_e280,e280,e280_near,
 e280_near_distance,e280_ngbox_avg)=modeldata.extract_all() # extract the data from all the models

nmodels=len(modelnames_eoi400) # we also have HadISST in e280
#######################################
# write data out to a workbook
write_to_book(outputfile,longitudes,latitudes,longitudes_alt, sitenames)

###################################
# plot model points
#for model in range(0,len(modelnames_eoi400)):
    #plotpoints(lonlist,latlist,eoi400[model,:]-e280[model,:])

#    plotpoints(longitudes_alt,latitudes,eoi400[model,:])

