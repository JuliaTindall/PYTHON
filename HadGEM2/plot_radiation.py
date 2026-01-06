#!/usr/bin/env python2.7
#NAME
#    PLOT_RADIATION
#PURPOSE
#    This program will plot the radiation budget for the pliocene simulations
#
# search for 'main program' to find end of functions
# Julia 22/11/2016



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess


#functions are:
#  def plotdata
#  def annmean
#  def seasmean

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
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
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.),extend="both")
        cbar = plt.colorbar(cs,orientation="horizontal")
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu_r',extend="both")
            cbar = plt.colorbar(cs,orientation="horizontal")

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend="both")
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                print(np.shape(plotdata))
                cs = map.contourf(x,y,plotdata,V,extend="both")
                cbar = plt.colorbar(cs,orientation="horizontal",)

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
#end def plotdata

def annmean(switch,HadCM3,expt,extra):
    # switch is a dummy variable to allow the program to be called

    # we will plot a) incoming sw ra flux (toa) field200
    #              b) incoming sw ra flux (toa) field201
    #              c) outgoing lw rad flux (toa) olr
    # other things if the budgets don't balance

    if HadCM3 == 'y':
        datasetname='/nfs/hera1/earjcti/um/netcdf/'+expt+'_netcdf/'+expt+'a@pd'+extra+'[7-9]*.nc'
        print(datasetname)
        f=MFDataset(datasetname)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        insw=f.variables['field200'][:]
        outsw=f.variables['field201'][:]
        outlw=f.variables['olr'][:]
        netlw= (-1.0) *outlw
        netsw=insw-outsw


        titlePI='PI-TAnn_HadCM3'
        titlediff='Plio - PI Tanom_HadCM3'
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt+'/netcdf/'
        filename=expt+'a@pd'+extra+'[7-9]*.nc'
        os.system('ls '+dirname+filename)
        os.system('ls '+dirname+filename+' | wc -l')
        nval=os.system('ls '+dirname+filename+' | wc -l')
        allfiles=subprocess.check_output('ls '+dirname+filename+' | wc -l',shell=True)
        print('number of files=',allfiles)
        f=MFDataset(dirname+filename)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        # net downward sw flux at the trop
        netsw=f.variables['solar_2'][:]
        # net downward lw flux at the trop
        netlw=f.variables['longwave_2'][:]

        # alternative 1
        insw=f.variables['field200'][:]
        outsw=f.variables['field201'][:]
        outlw=f.variables['olr'][:]
        netlw= (-1.0) *outlw
        netsw=insw-outsw


    netsw=np.squeeze(netsw)
    netlw=np.squeeze(netlw)
    ntimes,ny,nx=np.shape(netsw)
    print(ntimes,ny,nx)
    
#average across the time dimension
    netsw_ann=np.mean(netsw,axis=0)
    netlw_ann=np.mean(netlw,axis=0)
    
    
    plt.figure(0)
    wm2='w/m2'
    lontemp=lon
    titlename=expt+'net SW ann'
    netsw_ann,lon = shiftgrid(180.,netsw_ann,lon,start=False)    
    plotdata(netsw_ann,0,lon,lat,titlename,0,400,40.0,0.0,'n',wm2)


    lon=lontemp
    netlw_ann,lon = shiftgrid(180.,netlw_ann,lon,start=False)    
    plotdata(netlw_ann * (-1.0),1,lon,lat,' net lw ann',100,300,40.0,0.0,'n',wm2)

    # calculate residual and mean residual weighted by cos latitude

    residual=netsw_ann+ netlw_ann
  
    weights=np.cos(np.radians(lat))
    print('len weights',len(weights))
    resid_zon=np.average(residual,axis=0,weights=weights)
    average_residual=np.average(resid_zon)

    print(average_residual)
    titlename='residual, avg='+str(average_residual)+'(w/m2)'
    plotdata(residual,2,lon,lat,titlename,-100,100,20.0,0.0,'n',wm2)
    


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_radiation/annmean_'+expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


#end def annmean



################################
# main program

# annual mean
figureno=0


#exptname='xkvjg'  #xkvje xkvjf xkvjg
#extra='n'
#HadCM3='n'
#exptname='ximut'
#extra='l'
#plt.figure(figureno)
#annmean('y',HadCM3,exptname,extra)
#figureno=figureno+1


HadCM3='y'
exptname='xhcph'
extra='j'
plt.figure(figureno)
annmean('y',HadCM3,exptname,extra)
figureno=figureno+1





#djf mean
#plt.figure(figureno)
#seasmean('dc','ja','fb',figureno,'djf',HadCM3)
#figureno=figureno+1

#jja mean
#plt.figure(figureno)
#seasmean('jn','jl','ag',figureno,'jja',HadCM3)
#figureno=figureno+1


sys.exit(0)

####

