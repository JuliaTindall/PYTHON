#!/usr/bin/env python2.7
#NAME
#    PLOT_SALINITY
#PURPOSE
#    This program will do all the plots to do with salinity
#    we will begin by plotting salinity by depth across the Atlantic
#    However we will move on to other things to do with satlinity
# search for 'main program' to find end of functions
# Julia 19/1/2017



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
#  def annmean
#  def seasmean

# functions start here
def plotmap(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
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
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal",extend='max')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu_r')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='max')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                print(np.shape(plotdata))
                cs = map.contourf(x,y,plotdata,V)
                cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
 
#end def plotdata
#####################################
def plot_lat_dep(plotdata,fileno,lat,dep,titlename,minval,maxval,valinc,cbarname):
    lats, deps  = np.meshgrid(lat,dep)

    V=np.arange(minval,maxval,valinc)
    
    print(np.shape(plotdata))
    cs = plt.contourf(lats,deps,plotdata,V,extend="both")
    plt.gca().invert_yaxis()

    cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
 
#end def plotdata

#  to check if a character is numeric
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#============================
def Atlantic_salinity_depth_plot(expt_name):



#  get data from files

    f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pg2/'+expt_name+'o@pg*.nc')
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    dep = f.variables['depth_1'][:]
    salin=f.variables['salinity'][:]
    salin=np.squeeze(salin)
    f.close()

    ntimes,ndep,nlat,nlon=np.shape(salin)
    print(ntimes,ndep,nlat,nlon)

    mean_salin=np.mean(salin,axis=0)


    print(np.shape(mean_salin))

#  get mask of Oceans
    filemask='/nfs/see-fs-02_users/earjcti/MOC/merid/basin_hadgom_216'
    f1=open(filemask,'r')
    # discard 3 title lines
    textline=f1.readline()
    textline=f1.readline()
    textline=f1.readline()
    # next line is details from basin
    nrows=nlat
    Indian_d=np.zeros((nrows,4))
    Pacific_d=np.zeros((nrows,4))
    Atlantic_d=np.zeros((nrows,4))
    Combined_d=np.zeros((nrows,4))
    for line in f1:
        linesplit=line.split()  # split line by space
        if is_number(linesplit[0]):
            rowno=int(linesplit[0])
            Indian_d[rowno-1,:]=int(linesplit[1]),int(linesplit[2]),\
                int(linesplit[3]),int(linesplit[4])
            Pacific_d[rowno-1,:]=int(linesplit[5]),int(linesplit[6]),\
                int(linesplit[7]),int(linesplit[8])
            Atlantic_d[rowno-1,:]=int(linesplit[9]),int(linesplit[10]),\
                int(linesplit[11]),int(linesplit[12])
            Combined_d[rowno-1,:]=int(linesplit[13]),int(linesplit[14]),\
                int(linesplit[15]),int(linesplit[16])

    f1.close()

    Indian_mask=np.ones((nlat,nlon),dtype=bool)
    Pacific_mask=np.ones((nlat,nlon),dtype=bool)
    Atlantic_mask=np.ones((nlat,nlon),dtype=bool)
    Combined_mask=np.ones((nlat,nlon),dtype=bool)
    for j in range(0,nrows):
        Indian_mask[j,Indian_d[j,0]:Indian_d[j,1]]=0
        Indian_mask[j,Indian_d[j,2]:Indian_d[j,3]]=0
        Pacific_mask[j,Pacific_d[j,0]:Pacific_d[j,1]]=0
        Pacific_mask[j,Pacific_d[j,2]:Pacific_d[j,3]]=0
        Atlantic_mask[j,Atlantic_d[j,0]:Atlantic_d[j,1]]=0
        Atlantic_mask[j,Atlantic_d[j,2]:Atlantic_d[j,3]]=0
        Combined_mask[j,Combined_d[j,0]:Combined_d[j,1]]=0
        Combined_mask[j,Combined_d[j,2]:Combined_d[j,3]]=0



    # get mean salinity over indian ocean
    Indian_salin=np.ma.masked_array(mean_salin,mask=np.tile(Indian_mask,(mean_salin.shape[0],1)))
    Pacific_salin=np.ma.masked_array(mean_salin,mask=np.tile(Pacific_mask,(mean_salin.shape[0],1)))
    Atlantic_salin=np.ma.masked_array(mean_salin,mask=np.tile(Atlantic_mask,(mean_salin.shape[0],1)))
    Combined_salin=np.ma.masked_array(mean_salin,mask=np.tile(Combined_mask,(mean_salin.shape[0],1)))

    lontemp=lon
    Indian_salin,lon = shiftgrid(180.,Indian_salin,lon,start=False)
    lon=lontemp
    Pacific_salin,lon = shiftgrid(180.,Pacific_salin,lon,start=False)
    lon=lontemp
    Atlantic_salin,lon = shiftgrid(180.,Atlantic_salin,lon,start=False)
    lon=lontemp
    Combined_salin,lon = shiftgrid(180.,Combined_salin,lon,start=False)


   # plotmap(Indian_salin[0,:,:]*1000.,0,lon,lat,'test Ind',-4,4,0.5,0.0,'n','psu')
   # plotmap(Pacific_salin[0,:,:]*1000.,1,lon,lat,'test Pac',-4,4,0.5,0.0,'n','psu')
   # plotmap(Atlantic_salin[0,:,:]*1000.,2,lon,lat,'test Atl',-4,4,0.5,0.0,'n','psu')
   # plotmap(Combined_salin[0,:,:]*1000.,3,lon,lat,'test Com',-4,4,0.5,0.0,'n','psu')
   # plt.show()


    # julia note - atlantic basin is fine - not sure about other ones


    # now average the salinity across the latitude over the Atlantic basin

    print('Atl Salin',np.shape(Atlantic_salin))
    Atlantic_sal_avg=np.mean(Atlantic_salin,axis=2)
    print('Atl salin mean',np.shape(Atlantic_sal_avg),np.shape(lat),np.shape(dep))
    

    titlename='Atlantic Salinity: '+expt_name
    plot_lat_dep((Atlantic_sal_avg*1000.)+35.0,0,lat,dep,titlename,33,36,0.2,'psu')


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_salin/'+expt_name+'_salinity.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

    return(lat,dep,Atlantic_sal_avg)
    

#end def Atlantic_salinity_depth_plot
#end def annmean



################################
# main program

# annual mean salinity by depth

retdata=Atlantic_salinity_depth_plot('xkvjg')
lat=retdata[0]
dep=retdata[1]
xkvjg_salinity=retdata[2]

retdata=Atlantic_salinity_depth_plot('xkvjf')
xkvjf_salinity=retdata[2]

retdata=Atlantic_salinity_depth_plot('xkvje')
xkvje_salinity=retdata[2]


# plot lat depthanomalies
titlename='Atlantic Salinity xkvjg - xkvje'
sal_anom=xkvjg_salinity-xkvje_salinity
plot_lat_dep((sal_anom*1000.),0,lat,dep,titlename,-0.5,0.6,0.1,'psu')
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_salin/xkvjg-xkvje_salinity_anom.eps' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

titlename='Atlantic Salinity xkvjf - xkvje'
sal_anom=xkvjf_salinity-xkvje_salinity
plot_lat_dep((sal_anom*1000.),0,lat,dep,titlename,-0.5,0.6,0.1,'psu')
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_salin/xkvjf-xkvje_salinity_anom.eps' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()



sys.exit(0)

####

