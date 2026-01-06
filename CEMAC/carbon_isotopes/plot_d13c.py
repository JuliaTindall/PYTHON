#!/usr/bin/env python2.7
#NAME
#    PLOT_d13c.py
#PURPOSE
#    This program will plot d13c from ocean tracer 15 in one of jennies pg files
#
# search for 'main program' to find end of functions
# Julia 28/7/2016



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
#  def plotdata
#  def surf_lat_lon

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
   


#end def plotdata

def surf_lat_lon(filein):
    # switch is a dummy variable to allow the program to be called
    #==============
    # preindustrial


    # read in data from multiple files
    print(filein)
    f=Dataset(filein)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['otracer14_2'][:]
    atemp=np.squeeze(atemp)
    nz,ny,nx=np.shape(atemp)
    print(nz,ny,nx)
    
#average across the surface layer
    c13_modelunits=atemp[0,:,:]
    c13_modelunits_10=atemp[15,:,:]

    d13C=(c13_modelunits-100.)*10. 
    d13C_10=(c13_modelunits_10-100.)*10. 

    
    plt.figure(0)
    lontemp=lon
    d13C,lon = shiftgrid(180.,d13C,lon,start=False)
    lon=lontemp
    d13C_10,lon = shiftgrid(180.,d13C_10,lon,start=False)
    
    print('about to plot')

    plotdata(d13C,0,lon,lat,'d13C',0,3,0.2,0.0,'n','permille')
    plotdata(d13C_10,1,lon,lat,'d13C_10',0,3,0.2,0.0,'n','permille')
    plt.show()
    print('plotted first')
    f.close()

#    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surftemp/MAT_anom_only_'+pliop2_expt+'.eps' 
#    plt.savefig(fileout, bbox_inches='tight')  

#    plt.close()

 

#end def annmean





#end def annmean

################################
# main program

# annual mean
figureno=0
# biology only
#filein='/nfs/see-fs-02_users/earjcti/temporary/jennie_data/xnbqpo#pg000002073c1+.nc'
# biology and air sea
#filein='/nfs/see-fs-02_users/earjcti/temporary/jennie_data/xnbqoo#pg000002038c1+.nc'
# biology and air sea after one year
filein='/nfs/see-fs-02_users/earjcti/temporary/jennie_data/xnbqoo#pg000000001c1+.nc'
surf_lat_lon(filein)


sys.exit(0)

####

