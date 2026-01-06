#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#"""
#Created on Fri Jul  5 15:11:26 2019
#
#@author: earjcti
#"""
#
#   This program will obtain the SST data from the pliovar site locations 
#   from Aislings PlioMIP1 data
#
#
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

        if datafile_ == 'Erin':
            self.filenamein = filename + 'pliovar_metadata_global_alan.csv'
            self.latcolumn = 2
            self.loncolumn = 3
            self.outend = 'PlioMIP1output_pliovar.xls'
        if datafile_ == 'Harry':
            self.filenamein = filename + 'Copy_of_CSCD_localities.csv'
            self.latcolumn = 1
            self.loncolumn = 2
            self.outend = 'PlioMIP1output_CSCD_localities.xls'

        if self.linuxwin == 'l':
            self.fileout = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/'+self.outend
        else:
            self.fileout = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\proxydata\\'+self.outend


    def read_file(self):
        count = 0
        latlist = []
        lonlist = []
        lonlist_alt = []

        with open(self.filenamein) as csvfile:
            readCSV = csv.reader(csvfile, delimiter=',')
            for row in readCSV:
                if count != 0: # not titleline
                    latlist.append(np.float(row[self.latcolumn]))
                    lon = (np.float(row[self.loncolumn]))
                    lonlist_alt.append(lon)
                    if lon < 0.: # longitude goes from 0-360 in models
                        lonlist.append(lon+360.)
                    else:
                        lonlist.append(lon)
                count = count+1
        returndata = self.fileout,lonlist, latlist, lonlist_alt
        return returndata

# end of class getinitdata
###############################################################################

class Getmodeldata:
    # get all of the data from the model at the required gridpoints
    def __init__(self, test, field, latlist, lonlist, lonlist_alt, period):

        fieldunits = {
            "SST" : "degC",
            "TotalPrecipitation" : "mm/day"
                        }


        self.fieldnames = field
        self.latlist = latlist
        self.lonlist = lonlist
        self.lonlist_alt=lonlist_alt
        self.period = period

        if test == 'y':
            self.modelnames = ['COSMOS']
        else:
            self.modelnames = ['COSMOS', 'GISS', 'HAD', 
                               'IPSL', 'MIROC', 'MRI',
                               'NOR'
                               ]

        self.units = fieldunits.get(fieldnames)


        #
    def extract_model_points(self, allcube_,fieldreq_):
        """
        will extract the data at each point from 'filenameuse'

        calls: get_near_data
        """

        
        ncubes = len(allcube_)
        for i in range(0, ncubes):
            if allcube_[i].var_name == fieldreq_:
                cube = allcube_[i]

      
        print(cube,'found')
        cubelats = cube.coord('latitude').points
        cubelons = cube.coord('longitude').points

        model_data = np.zeros(len(self.lonlist))
        model_data_near = np.zeros(len(self.lonlist)) # values near the point
        near_distance = np.zeros(len(self.lonlist)) # how far away we have to look to get data
        ngbox_avg = np.zeros(len(self.lonlist)) # how many gridboxes we are averaging over to get data

        for i in range(0,len(self.lonlist)):
            # find nearest latitude and lontiude to the value
            latix = (np.abs(cubelats-self.latlist[i])).argmin()
            lonix = (np.abs(cubelons-self.lonlist_alt[i])).argmin()
            
            print(self.lonlist_alt[i],self.latlist[i],latix,lonix)

            # get data from this location
            data_slice  =  cube.extract(iris.Constraint(
                        latitude = cubelats[latix],longitude = cubelons[lonix]))
            model_data[i] = data_slice.data
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
            data_near = np.nan # if no data near set to nan

        return data_near,count_finite

    def get_fieldreq(self,model_):
        
        PeriodE280Use={
                       "COSMOS" : "Ctrl",
                       "GISS" : "Ctrl",
                       "HAD" : "ctrl",
                       "IPSL" : "ctrl",
                       "MIROC" : "ctrl",
                       "MRI" : "ctrl",
                       "NOR" : "ctrl"
                       }
        
        PeriodEoi400Use={
                         "COSMOS" : "Plio",
                         "GISS" : "Plio",
                         "HAD" : "plio",
                         "IPSL" : "plio",
                         "MIROC" : "plio",
                         "MRI" : "plio",
                         "NOR" : "plio"
                         }
        
        fieldname={
                   "COSMOS" : "SST",
                   "GISS" : "SST",
                   "HAD" : "sst",
                   "IPSL" : "sst",
                   "MIROC" : "sst",
                   "MRI" : "sst",
                   "NOR" : "sst"
                   }
         
        if self.period == 'E280':
            self.fieldreq = (model_ + '_' + 
                        PeriodE280Use.get(model_)+
                        '_' + fieldname.get(model_))
        if self.period == 'EOI400':
            self.fieldreq = (model_ + '_' + 
                        PeriodEoi400Use.get(model_)+
                        '_' + fieldname.get(model_)) 
            
            

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
          
            filename = ('/nfs/hera1/earjcti/PLIOMIP/PlioMIP1_regridded.nc')
            self.get_fieldreq(self.modelnames[model])
            print('fieldreq is',self.fieldreq)
            
            allcube = iris.load(filename)
            
            # get model points and how far away they are from data
            (sitevals[model,:],sitenear[model,:],sitenear_dist[model,:],
                sitenear_ngbox_avg[model,:]) = self.extract_model_points(allcube,self.fieldreq)

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



def write_sheet(wb, style, sheetname, modelnames, lonlist_alt, latlist, datawrite):
    sheet  =  wb.add_sheet(sheetname)
    sheet.write(0,0,'lat')
    sheet.write(0,1,'lon')
    for model in range(0,nmodels):
        sheet.write(0,2+model,modelnames[model])

    for i in range(0,npoints):
        sheet.write(i+1,0,latlist[i],style)
        sheet.write(i+1,1,lonlist_alt[i],style)
        for model in range(0,nmodels):
            sheet.write(i+1,2+model,datawrite[model,i],style)

    # add multimodel mean and multimodel standard deviation
    sheet.write(0,2+nmodels,'MMM')
    sheet.write(0,3+nmodels,'MM-SD')
    for i in range(0,npoints):
        sheet.write(i+1,2+nmodels,np.nanmean(datawrite[0:nmodels,i]),style)
        sheet.write(i+1,3+nmodels,np.nanstd(datawrite[0:nmodels,i]),style)

    # add extra columns if we have them  this is likely to be hadisst
    if len(modelnames) > nmodels:
        print(modelnames,nmodels)
        for model in range(nmodels,len(modelnames)):
            sheet.write(0,4+model,modelnames[model])
            for i in range(0,npoints):
                sheet.write(i+1,4+model,datawrite[model,i],style)



######################################################################################
def write_to_book(fileout,lonlist,latlist,lonlist_alt):
    # write to workbook
    # calls write_sheet

    # Workbook is created
    wb  =  Workbook()

    style  =  xlwt.XFStyle()
    style.num_format_str  =  '0.00'


    # add_sheet for Eoi400 E280 and difference
    write_sheet(wb,style, 'EOI400', modelnames_eoi400, lonlist_alt, latlist, eoi400)
    write_sheet(wb,style, 'E280', modelnames_e280, lonlist_alt, latlist,e280)
    write_sheet(wb,style, 'EOI400-E280',modelnames_eoi400,lonlist_alt,latlist,
                eoi400[0:nmodels]-e280[0:nmodels])

    # add_sheet near Eoi400 E280 and difference
    write_sheet(wb,style, 'EOI400near',modelnames_eoi400,lonlist_alt,latlist,eoi400_near)
    write_sheet(wb,style, 'E280near',modelnames_e280,lonlist_alt,latlist,e280_near)
    write_sheet(wb,style, 'EOI400-E280near',modelnames_eoi400,lonlist_alt,latlist,
                eoi400_near[0:nmodels]-e280_near[0:nmodels])

    # add sheet for how far away we need to look for data
    write_sheet(wb,style, 'EOI400distance',modelnames_eoi400,lonlist_alt,latlist,eoi400_near_distance)
    write_sheet(wb,style, 'E280distance',modelnames_e280,lonlist_alt,latlist,e280_near_distance)

    # add sheet for how many gridboxes we are averaging over
    write_sheet(wb,style, 'EOI400nboxes',modelnames_eoi400,lonlist_alt,latlist,eoi400_ngbox_avg)
    write_sheet(wb,style, 'E280nboxes',modelnames_e280,lonlist_alt,latlist,e280_ngbox_avg)




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
datafile = 'Harry' # could have Harry or Erin
testdata = 'n'   # yes use one model no use full range of models
fieldnames = 'SST'

indata = Getinitialdata(linuxwin,datafile)
outputfile,longitudes,latitudes,longitudes_alt = indata.read_file() # get the lats lons required and the number of sites
npoints = len(longitudes) # get the number of points


######################################
# setup a map and plot the points

#if linuxwin = ='l':
#    plotpoints(lonlist_alt,latlist,np.zeros(count))


##############################
# get the SST data from IRIS cubes

#eoi400
modeldata=Getmodeldata(testdata,fieldnames,
                       latitudes,longitudes,longitudes_alt,'EOI400')
(modelnames_eoi400,eoi400,eoi400_near,
 eoi400_near_distance,eoi400_ngbox_avg)=modeldata.extract_all()


modeldata=Getmodeldata(testdata,fieldnames,
                       latitudes,longitudes,longitudes_alt,'E280')
(modelnames_e280,e280,e280_near,
 e280_near_distance,e280_ngbox_avg)=modeldata.extract_all() # extract the data from all the models

nmodels=len(modelnames_eoi400) # we also have HadISST in e280
#######################################
# write data out to a workbook
write_to_book(outputfile,longitudes,latitudes,longitudes_alt)

###################################
# plot model points
for model in range(0,len(modelnames_eoi400)):
    #plotpoints(lonlist,latlist,eoi400[model,:]-e280[model,:])

    plotpoints(longitudes_alt,latitudes,eoi400[model,:])

