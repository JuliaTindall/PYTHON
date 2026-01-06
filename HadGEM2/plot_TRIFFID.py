#!/usr/bin/env python2.7
#NAME
#    PLOT_TRIFFID
#PURPOSE
#    This program will plot the vegetation types
# search for 'main program' to find end of functions
# Julia 31/1/2019



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid


#functions are:
#  def plotdata
#  def annmean
#  def seasmean

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,ygrid,xgrid,yuse,xuse,yspan,xspan):


    if fileno > 3:
        print('ERROR NOT ENOUGH SPACE ON PAGE ',fileno)
        sys.exit()
    #plt.subplot2grid((10,12),(fileno*2,0),colspan=9,rowspan=2)

    print(ygrid,xgrid,yuse,xuse,yspan,xspan)
    plt.subplot2grid((ygrid,xgrid),(yuse,xuse),colspan=yspan,rowspan=xspan)

    lons, lats = np.meshgrid(lon,lat)
   
    northlat=90.0
    southlat=-90.0

    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=southlat,urcrnrlat=northlat,projection='cyl',resolution='c',fix_aspect=False)
    x, y = map(lons, lats)
    map.drawcoastlines()

    plotdata2=plotdata


    # plot limits
    if V == 0:
        V=np.arange(minval,maxval+valinc,valinc)

    # plot map
        
    if uselog == 'n':
        mycmap=mp.cm.get_cmap('jet',len(V+2))
        newcolors=mycmap(np.linspace(0,1,len(V+2)))
        white=([1,1,1,1])
        newcolors[0,:]=white
        mycmap=ListedColormap(newcolors)
        cs = map.contourf(x,y,plotdata,V,cmap=mycmap)
      
        #cs = map.contourf(x,y,plotdata,V,cmap='YlGnBu')
      
        #cs = map.contourf(x,y,plotdata,V,cmap='rainbow')
        if fileno == 0:
            plt.title('mPWP')
        if fileno == 1:
            plt.title('PI')
    else:
        mycmap=mp.cm.get_cmap('RdBu_r',len(V+2))
        newcolors=mycmap(np.linspace(0,1,len(V+2)))
        white=([1,1,1,1])
        newcolors[(len(V)/2)-2:(len(V)/2)+3,:]=white
        mycmap=ListedColormap(newcolors)
        cs = map.contourf(x,y,plotdata,V,cmap=mycmap,extend='none')
       
  
    fontsize=10
    if fileno == 0 or fileno ==2:
        plt.text(-180.0-6,northlat-fontsize-1,titlename,fontsize=fontsize,ha="right",bbox=dict(boxstyle="square,pad=0.1",color="white"))
 

    # colorbar
    if fileno==0:

        #plt.subplot2grid((ygrid,xgrid),(yuse,xuse),colspan=yspan,rowspan=xspan)

        plt.subplot2grid((ygrid,xgrid),(ygrid-1,xuse),colspan=4,rowspan=xspan)

        plt.gca().set_visible(False)
        cbar = plt.colorbar(cs,orientation="horizontal",fraction=1.0,format='%0.1f')
        cbar.set_label(cbarname)
        cbar.ax.tick_params(labelsize=7)
         

    plotdata=plotdata2
    
    # plot map boundary
    map.drawmapboundary



#end def plotdata



def annmean(preind_expt,plio_expt,HadCM3):

    if HadCM3 == 'n':
        filestart='/nfs/hera1/earjcti/um/HadGEM_data/'
        filemid='/netcdf/pifiles/'
        modname='HadGEM2'


    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'
        filemid='_netcdf/'
        modname='HadCM3'

    if HadCM3 == 'HGpeter':
        filestart='/nfs/hera1/earjcti/um/HadGEM_data/'
        filemid='/netcdf/pifiles/'
        modname='HadGEM2_peter'
    

    #==============
    # preindustrial

    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    vegdata=0.
    count=0

    for mon in range (0,len(monthnames)):
       
        f=Dataset(filestart+preind_expt+filemid+preind_expt+'a@pin99'+monthnames[mon]+'.nc')
    
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        ntypes=f.variables['pseudo_1'][:]
        atemp=f.variables['field1391'][:]
        atemp_mean=np.mean(atemp,axis=0)
        atemp_mean=np.squeeze(atemp_mean)

        vegdata=vegdata+atemp_mean
        count=count+1

    vegdata_pi=vegdata/count # get mean of vegdata

     #==============
     # Pliocene

    vegdata=0.
    count=0
    for mon in range (0,len(monthnames)):
        print(filestart+plio_expt+filemid+plio_expt+'a@piy99'+monthnames[mon]+'.nc')
        if HadCM3 =='HGpeter':
            f=Dataset('/nfs/hera1/earjcti/um/ximup/datam/ximupa@pik99'+monthnames[mon]+'.nc')
        else:
            f=Dataset(filestart+plio_expt+filemid+plio_expt+'a@pio00'+monthnames[mon]+'.nc')
    
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        ntypes=f.variables['pseudo_1'][:]
        atemp=f.variables['field1391'][:]
        atemp_mean=np.mean(atemp,axis=0)
        atemp_mean=np.squeeze(atemp_mean)

        vegdata=vegdata+atemp_mean
        count=count+1

    vegdata_plio=vegdata/count # get mean of vegdata
     
    broadleaf_pi=vegdata_pi[0,:,:]
    broadleaf_plio=vegdata_plio[0,:,:]
    needleleaf_pi=vegdata_pi[1,:,:]
    needleleaf_plio=vegdata_plio[1,:,:]
    grasses_pi=vegdata_pi[2,:,:]+vegdata_pi[3,:,:]+vegdata_pi[4,:,:]
    grasses_plio=vegdata_plio[2,:,:]+vegdata_plio[3,:,:]+vegdata_plio[4,:,:]
    soil_pi=vegdata_pi[7,:,:]
    soil_plio=vegdata_plio[7,:,:]

    # shiftgrid

    lontemp=lon
    broadleaf_pi,lon = shiftgrid(180.,broadleaf_pi,lon,start=False)
   
    lon=lontemp
    broadleaf_plio,lon = shiftgrid(180.,broadleaf_plio,lon,start=False)
   
    lon=lontemp
    needleleaf_pi,lon = shiftgrid(180.,needleleaf_pi,lon,start=False)
   
    lon=lontemp
    needleleaf_plio,lon = shiftgrid(180.,needleleaf_plio,lon,start=False)
   

    lon=lontemp
    grasses_pi,lon = shiftgrid(180.,grasses_pi,lon,start=False)
   
    lon=lontemp
    grasses_plio,lon = shiftgrid(180.,grasses_plio,lon,start=False)
   
    lon=lontemp
    soil_pi,lon = shiftgrid(180.,soil_pi,lon,start=False)
   
    lon=lontemp
    soil_plio,lon = shiftgrid(180.,soil_plio,lon,start=False)
   

   
    plotdata(broadleaf_plio,0,lon,lat,'BLT',0.,1,0.05,0,'n','fractional coverage',17,7,0,0,2,4)
    plotdata(needleleaf_plio,2,lon,lat,'NLT',0.,1,0.05,0,'n','fractional coverage',17,7,4,0,2,4)
    plotdata(grasses_plio,2,lon,lat,'Grasses',0.,1,0.05,0,'n','fractional coverage',17,7,8,0,2,4)
    plotdata(soil_plio,2,lon,lat,'Soil',0.,1,0.05,0,'n','fractional coverage',17,7,12,0,2,4)
   
    plotdata(broadleaf_pi,1,lon,lat,'BLT',0.,1,0.05,0,'n','fractional coverage',17,7,0,2,2,4)
    plotdata(needleleaf_pi,3,lon,lat,'NLT',0.,1,0.05,0,'n','fractional coverage',17,7,4,2,2,4)
    plotdata(grasses_pi,3,lon,lat,'Grasses',0.,1,0.05,0,'n','fractional coverage',17,7,8,2,2,4)
    plotdata(soil_pi,3,lon,lat,'Soil',0.,1,0.05,0,'n','fractional coverage',17,7,12,2,2,4)
   
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_TRIFFID/'+modname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_TRIFFID/'+modname+'.png' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()


 

#end def annmean




################################
# main program

# annual mean

HadCM3='n'
preind_expt='xkvje'
pliop2_expt='xkvjg'

annmean(preind_expt,pliop2_expt,HadCM3)

#HadCM3='y'
#preind_expt='xiboi'
#plio_expt='xibol'

#annmean(preind_expt,plio_expt,HadCM3)

#HadCM3='n'
#preind_expt='xkvje'
#pliop2_expt='ximup'
#HadCM3='HGpeter' # check against peters vegetation#

#annmean(preind_expt,pliop2_expt,HadCM3)




sys.exit(0)

####

