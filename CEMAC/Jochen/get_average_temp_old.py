#!/usr/bin/env python2.7
#NAME
#    get_average_temp.py
#PURPOSE
#    This program will get the long term annual average temperature from
#    both individual and database files.
#
#    It was originally written for Jochen who wanted to verify some values
#    in a spreadsheet of Fergus
#
# search for 'main program' to find end of functions
# Julia 3/8/2017



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid


#functions are:
#  def plotdata  (not currently used)
#  def annmean_ind
#  def annmean_database

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    if fileno != 99:
        plt.subplot(2,2,fileno+1)

   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe

    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary

    x, y = map(lons, lats)

    map.drawcoastlines()
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal")
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu_r')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='max')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='i': #increasing
                    cs = map.contourf(x,y,plotdata,V,norm=mp.colors.LogNorm(vmin=0,vmax=32),cmap='Reds')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    cs = map.contourf(x,y,plotdata,V,extend='both',cmap='rainbow')
                    cbar = plt.colorbar(cs,orientation="horizontal")


    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        cbar.set_label(cbarname,labelpad=-70,size=20)
        cbar.ax.tick_params(labelsize=20)
        plt.title(titlename,loc='left',fontsize=20)
   


#end def plotdata

def annmean_ind(exptname,extra,startyear,nyears):

    # if nyears is a negative number read in all the files in the directory
    if nyears < 0:
        f=MFDataset('/nfs/hera2/apps/metadata/experiments/'+exptname+'/netcdf/'+exptname+'a@pd*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]

        atemp=f.variables['temp'][:] # this is temp at 1.5m
        atemp=np.squeeze(atemp)
        ntimes,ny,nx=np.shape(atemp)
    
        btemp=f.variables['temp_1'][:] # surface temperature after timestep
        btemp=np.squeeze(btemp)
        f.close()

    else:
        print('you need to set up for a subset of the directory')
        sys.exit()
        
    
    #average across the time dimension
    temp_1point5=np.mean(atemp,axis=0)
    temp_surf=np.mean(btemp,axis=0)

    # create weighting array
    weightarr=np.zeros(np.shape(temp_surf))
    for i in range(0,len(lon)):
        weightarr[:,i]=np.cos(np.deg2rad(lat))

    # find weighted mean

    avg_temp_surf=np.average(temp_surf,weights=weightarr)
    avg_temp_1point5=np.average(temp_1point5,weights=weightarr)

    print( )
    print('number of files is',ntimes,'for expt',exptname)
    print('Raw files average surface temperature=',avg_temp_surf-273.15)
    print('Raw files average 1.5m temperature=',avg_temp_1point5-273.15)

#end def annmean


def annmean_database(exptname,extra,startyear,nyears):

    f=Dataset('/nfs/hera2/apps/metadata/experiments/'+exptname+'/timeseries/'+exptname+'a@pd_SurfaceTemperature.nc')
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]

    atemp=f.variables['temp'][:] # this is temp at 1.5m
    atemp=np.squeeze(atemp)
    ntimes,ny,nx=np.shape(atemp)

    starttime=ntimes-(nyears*12)
    surface_temperature=atemp[starttime:ntimes,:,:]
       
    #average across the time dimension
    temp_surf=np.mean(surface_temperature,axis=0)

    # create weighting array
    weightarr=np.zeros(np.shape(temp_surf))
    for i in range(0,len(lon)):
        weightarr[:,i]=np.cos(np.deg2rad(lat))

    # find weighted mean

    avg_temp_surf=np.average(temp_surf,weights=weightarr)
   
    print('Database files average surface temperature ',exptname,'=',avg_temp_surf-273.15)

    return(avg_temp_surf-273.15)


#end def annmean



################################
# main program

# annual mean from individual files

exptname='xhckb'
extra='z'
startyear=70
nyears_ind=-10  # if negative number use all of them
nyears_db=30  # if negative number use all of them

annmean_ind(exptname,extra,startyear,nyears_ind)
database_tsurf=annmean_database(exptname,extra,startyear,nyears_db)


exptname='xiomv'
annmean_ind(exptname,extra,startyear,nyears_ind)
database_tsurf=annmean_database(exptname,extra,startyear,nyears_db)


exptname='xhckf'
annmean_ind(exptname,extra,startyear,nyears_ind)
database_tsurf=annmean_database(exptname,extra,startyear,nyears_db)


# experiments needed
# tdcza
# tdhzt
# tdlqa
# tdlqb
# tdlqc
# tdlqd
# tdlqe
# tdlqf
# tdlqg
# tdlqh
# tdlqi
# tdlqj
# tdlxf
# tdlxi
# xgraf
# xgrag
# xgrah
# xgrai
# xgygi
# xgygj
# xgygk
# xgygl
# xgygm
# xgygn
# xgygo
# xgygp
# xgygq
# xgygz
# xhckb
# xhckd
# xhckf
# xhckh
# xhcki
# xhckj
# xhckk
# xhckl
# xhckm
# xhckn
# xhcko
# xhckp
# xhckq
# xhckz
# xhgfk4
# xhsop4
# xhsor4
# xilug
# xiomv
# xjpli
# xhgfa
# tdipf
# tdkgg
# tdkgi
# tdkgf





sys.exit(0)

####

