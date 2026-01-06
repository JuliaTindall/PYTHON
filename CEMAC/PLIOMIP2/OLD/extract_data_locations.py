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
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import iris
################################
def plotpoints(lonlist,latlist):
# plot the points we have got from the file

    alllons=np.arange(-180,180,1)
    alllats=np.arange(-90,90,1)
    lons,lats=np.meshgrid(alllons,alllats)
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,
                urcrnrlat=90.0,projection='cyl',resolution='c')
    map.drawmapboundary
    x,y=map(lons,lats)
    map.drawcoastlines()


    for i in range(0,len(latlist)):
        lon=lonlist[i]
        lat=latlist[i] 
        # convert to map projection coords.
        # Note that lon,lat can be scalars, lists or numpy arrays.
        xpt,ypt = map(lon,lat)
        map.plot(xpt,ypt,'bo')  # plot a blue dot there
    plt.show()

###################################################################
def extract_model_points(lonlist,latlist,filename_eoi400,filename_e280):

# will extract the mpwp and the preindustrial data at each point
    eoi400cube=iris.load_cube(filename_eoi400)
    cubelats=eoi400cube.coord('latitude').points
    cubelons=eoi400cube.coord('longitude').points    

    eoi400_data=[]
    for i in range(0,len(lonlist)):
        # find nearest latitude and lontiude to the value
        
        latix=(np.abs(cubelats-latlist[i])).argmin()
        lonix=(np.abs(cubelons-lonlist[i])).argmin()
        print(latlist[i],lonlist[i],cubelats[latix],cubelons[lonix])

        # get data from this location
        data_slice = eoi400cube.extract(iris.Constraint(
                latitude=cubelats[latix],longitude=cubelons[lonix]))
        eoi400_data.append(data_slice.data)
        
    for i in range(0,len(lonlist)):
        print(i,lonlist[i],latlist[i],eoi400_data[i])


#################
# MAIN PROGRAM
################

###################################
# get the sitenames from a file

filename='/nfs/hera1/earjcti/PLIOMIP2/proxydata/pliovar_metadata_global_alan.csv'

count=0
latlist=[]
lonlist=[]

with open(filename) as csvfile:
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        if count!=0: # not titleline
            latlist.append(np.float(row[2]))
            lonlist.append(np.float(row[3]))
        count=count+1



######################################
# setup a map and plot the points

#plotpoints(lonlist,latlist)


##############################
# get the SST data from IRIS cubes

#modelnames=['HadCM3','NorESM-L','NorESM1-F','IPSLCM6A','IPSLCM5A2','IPSLCM5A',
#            'MIROC4m','COSMOS','UofT',
#            'EC-Earth3.1','MRI-CGCM2.3',
#            ]

filestart='/nfs/hera1/earjcti/regridded/'
modelnames=['HadCM3','NorESM-L']
fieldnames=['NearSurfaceTemperature']
units=['degC']

#fieldnames=['TotalPrecipitation']
#units=['mm/day']

for model in range(0,len(modelnames)):
    filename_eoi400=filestart+modelnames[model]+'/EOI400.SST.allmean.nc'
    filename_e280=filestart+modelnames[model]+'/E280.SST.allmean.nc'
    extract_model_points(lonlist,latlist,filename_eoi400,filename_e280)
