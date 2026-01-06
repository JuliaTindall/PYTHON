#!/usr/bin/env python2.7
#NAME
#    PLOT_D18O
#PURPOSE 
#    PLOT D18O FROM LOUISE/MAX TIMESLICE EXPERIMENTS from a single file
#
#
# Julia 20.09.2016


# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid



# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,cbartitle,minval,maxval,diffval):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(2,2,fileno+1)
   # map=Basemap(projection='robin',resolution='l')
    map=Basemap(llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=80,projection='mill')
   # map.drawmapboundary(fill_color='aqua')
    x, y = map(lons, lats)
    map.drawcoastlines()
    #V=np.arange(np.amin(plotdata),np.amax(plotdata),np.amin(plotdata)/10)
    V=np.arange(minval,maxval,diffval)
   # cs = map.contourf(x,y,plotdata,V)
    cs=map.contourf(x,y,plotdata,V)
    plt.title(titlename)
  
    cbar=map.colorbar(cs,location='bottom',pad="5%")
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_title(cbartitle)


#  functions end here


# 1.  Set up details.  Gridbox required and filename

#Jordan
latreq=37.5
longreq=30.0

#os.chdir("/home/users/jctindall/mholloway/xlubb/pcpd/")
os.chdir("/home/users/jctindall/maxoutput/")


#2. Set up filename and extract stash code 338
#   print d18o and dD for that file

#filename='iso_xlubba@pcq77.nc'
filename='iso_xlubba@pcr20_pcr49.nc'

f=Dataset(filename,mode='r')
f.dimensions
f.variables

lon = f.variables['longitude_1'][:]
lat = f.variables['latitude_1'][:]
times = f.variables['t'][:]
xsize=len(lon)
ysize=len(lat)
tsize=len(times)
dD=f.variables['dD'][:]
dD=np.squeeze(dD)
d18o=f.variables['dO18'][:]
d18o=np.squeeze(d18o)
h2o=f.variables['h2o'][:]
h2o=np.squeeze(h2o)


f.close()

# we need to get annual average dD and d18o by weighting by precipitation amount

dD_weight=dD*h2o
dD_weightsum=np.sum(dD_weight,axis=0)
d18o_weight=d18o*h2o
d18o_weightsum=np.sum(d18o_weight,axis=0)
tot_h2o=np.sum(h2o,axis=0)
dD_weightavg=dD_weightsum/tot_h2o
d18o_weightavg=d18o_weightsum/tot_h2o
#for t in range(0:tsize-1)
print("shape tpt h2o",tot_h2o.shape)
print("shape mean dDweight",dD_weightsum.shape)

cbartitle='mm'
titlename='avg'
plotdata(tot_h2o*60.*60.*24.*30/tsize,0,lon,lat,titlename,cbartitle,0,200,10)

print('amax',np.amax(dD))



#weighted dD
titlename='dD'
cbartitle='permil '
plotdata(dD_weightavg,1,lon,lat,titlename,cbartitle,-200,100,20)

#weighted d18o
titlename='d18o'
cbartitle='permil '
plotdata(d18o_weightavg,2,lon,lat,titlename,cbartitle,-50,10,2)

#weighted dxs
titlename='dxs'
cbartitle='permil '
dxs=dD_weightavg - (8.0 * d18o_weightavg)
plotdata(dxs,3,lon,lat,titlename,cbartitle,0,35,2)


plt.show()

# z end program here
