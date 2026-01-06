#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on Monday June 3rd 2019

#@author: earjcti
#
# Manuel has asked if I can extract some fields from Kanhu's experiments so 
# that he can do some averages on them.  He would like data from the 
# xkrax experiment which has time varying vegetation, greenhouse gases varying 
#
# Because he is comparing to GNIP he would like data from years 1957-2014 
# as these are the years where we have observations.  This is j57-k14 I think
#
#
########################################################
# other notes are

import os
import numpy as np
import scipy as sp
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris.analysis.cartography
from iris.experimental.equalise_cubes import equalise_attributes
from iris.util import unify_time_units
import cf_units as unit
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans
import sys



# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,
             minval,maxval,valinc,V,uselog,cbarname,lonmin,lonmax,latmin,latmax):
    lons, lats = np.meshgrid(lon,lat)
    if fileno !=99:        
        plt.subplot(1,2,fileno+1)



    map=Basemap(llcrnrlon=lonmin,urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax,projection='cyl',resolution='c')
    x, y = map(lons, lats)
    map.drawcoastlines()
    
    # set up for drawing gridlines
    if lonmax-lonmin <= 60:
        londiff=10
    else:
        londiff=30
     
   
    if latmax-latmin <= 60:
        latdiff=10
    else:
        latdiff=30

    parallels=np.arange(-90,90,latdiff)
    meridians=np.arange(-180,360,londiff)

   
    map.drawparallels(parallels,labels=[False,True,False,False])
    map.drawmeridians(meridians,labels=[False,False,False,True])


    plotdata2=plotdata
    #plotdata=maskoceans(x,y,plotdata)
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu',extend='both')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='both')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='ra':
                    cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    print(np.shape(plotdata))
                    if uselog=='nb':  # use bluescale
                        cmapuse='Blues'
                    else:
                        cmapuse='rainbow'
                    cs = map.contourf(x,y,plotdata,V,cmap=cmapuse,extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")


    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        cbar.set_label(cbarname,labelpad=-70,size=20)
        cbar.ax.tick_params(labelsize=20)
        plt.title(titlename,loc='left',fontsize=20)
   

    plotdata=plotdata2



#end def plotdata

##########################################################
# main program

# this is regridding where all results are in a single file
# create a dictionary with the long field names in and the field names we want
# we are also using dictionaries so that we only have to change timeperiod name
# when rerunning

linuxwin='l'
            	
filein={"l":'output/data_xkrax.nc',
           "w":'C:/Users/julia/OneDrive/WORK/DATA/TEMPORARY/xkraxa@pd'}



fieldname='d18o'
cube=iris.load_cube(filein.get(linuxwin),fieldname)

if fieldname=='d18o': # also get precipitation
    fieldp='TOTAL PRECIPITATION RATE     KG/M2/S'
    cubeprecip=iris.load_cube(filein.get(linuxwin),fieldp)
    
print(cube.coord('t').points)    
datajune= cube.data[5,:,:]
datajul= cube.data[6,:,:]
dataaug=cube.data[7,:,:]

# find average without weighting
datajja=(datajune+datajul+dataaug)/3.0

lon=cube.coord('longitude').points
lat=cube.coord('latitude').points

print(datajja)

plotdata(datajja,0,lon,lat,'d18o no weighting',-25.,1.,1.,0,'n',
         'permille',270.,330.,-60.,20.)


# find and plot average with weighting

precipjune= cubeprecip.data[5,:,:]
precipjul= cubeprecip.data[6,:,:]
precipaug=cubeprecip.data[7,:,:]

# find average without weighting
datajja_w=(((datajune * precipjune)+(datajul*precipjul)+(dataaug*precipaug))
        /(precipjune+precipjul+precipaug))

plotdata(datajja_w,1,lon,lat,'d18o weighting',-25.,1.,1.,0,'n',
         'permille',270.,330.,-60.,20.)

#plotdata(datajja,99,lon,lat,'d18o no weighting',-25.,1.,1.,0,'n',
#         'permille',0.,360.,-90.,90.)
plt.show()

sys.exit(0)
