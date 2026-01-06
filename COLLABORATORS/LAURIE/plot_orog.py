#!/usr/bin/env python2.7
#NAME
#    PLOT_OROG
#PURPOSE
#    This program will plot the orography for the eocene simulation
# search for 'main program' to find end of functions
# Julia 13/09/2017



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans


#functions are:
#  def plotdata
#  def annmean
#  def seasmean

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    if fileno !=99:
        plt.subplot(2,2,fileno+1)



   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
    x, y = map(lons, lats)
#    map.drawcoastlines()

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
                    cs = map.contourf(x,y,plotdata,V,extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")


    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        cbar.set_label(cbarname,labelpad=-70,size=20)
        cbar.ax.tick_params(labelsize=20)
        plt.title(titlename,loc='left',fontsize=20)
   

    plotdata=plotdata2

    #map.drawmapboundary

#end def plotdata



################################
# main program

# get orography
f=Dataset('/nfs/hera2/apps/metadata/experiments/tbomh/ancil/tbomh.qrparm.orog.nc','r')
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]
height=f.variables['ht'][:]
height=np.squeeze(height)
f.close()


# get land sea mask

f=Dataset('/nfs/hera2/apps/metadata/experiments/tbomh/ancil/tbomh.qrparm.mask.nc','r')
latm = f.variables['latitude'][:]
lonm = f.variables['longitude'][:]
mask=f.variables['lsm'][:]
mask=np.squeeze(mask)
f.close()
# replace 1's with zeros and zeros with 1's
mask=(mask - 1.0)*1.0

masked_height=np.ma.masked_array(height,mask=mask)

# plot it
masked_height,lon = shiftgrid(180.,masked_height,lon,start=False)
plotdata(masked_height,99,lon,lat,'HadCM3 Eocene Orography',0,3000,250,0.0,'n','height (m)')
plt.show()




sys.exit(0)

####

