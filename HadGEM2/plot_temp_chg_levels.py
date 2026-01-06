#!/usr/bin/env python2.7
#NAME
#    PLOT_temp_chg_levels.py
#PURPOSE
#    This program will plot the temperature change between the mPWP
#    and PI by level.  We will normalise it by the surface temperature
#    change
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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid


#functions are:
#  def plotdata
#  def annmean
#  def seasmean

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(2,2,fileno+1)

    print(fileno)

   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe
    map=Basemap(llcrnrlon=0.0,urcrnrlon=360.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='l')
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
                mycmap=mp.cm.get_cmap('bwr',len(V+2))
                newcolors=mycmap(np.linspace(0,1,len(V+2)))
                white=([1,1,1,1])
                newcolors[(len(V)/2)-1:(len(V)/2)+2,:]=white
                mycmap=ListedColormap(newcolors)
                print(mycmap)
                
                cs = map.contourf(x,y,plotdata,V,cmap=mycmap,extend="both")
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                print(np.shape(plotdata))
                cs = map.contourf(x,y,plotdata,V,cmap='rainbow',extend="both")
                cbar = plt.colorbar(cs,orientation="horizontal",)

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
#end def plotdata

def annmean(preind_expt,plio_expt,pliop2_expt,extra,HadCM3):
   

    if HadCM3 == 'y':
        filepi='/nfs/hera1/earjcti/um/'+preind_expt+'/netcdf/'+preind_expt+'a@pc'+extra+'[7-9]*.nc'
        fileplio='/nfs/hera1/earjcti/um/'+pliop2_expt+'/netcdf/'+pliop2_expt+'a@pc'+extra+'[7-9]*.nc'
        latname='latitude_1'
        lonname='longitude_1'
    else:
        filepi='/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/netcdf/'+preind_expt+'a@pc'+extra+'[7-9]*.nc'
        fileplio='/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/netcdf/'+pliop2_expt+'a@pc'+extra+'[7-9]*.nc'
        latname='latitude'
        lonname='longitude'

    # get preindustrial data
    f=MFDataset(filepi)
    lat = f.variables[latname][:]
    lon = f.variables[lonname][:]
    pressure = f.variables['p'][:]
    atemp=f.variables['temp'][:]
  
    temp_pi=np.mean(atemp,axis=0)
    temp_pi=np.squeeze(temp_pi)
    nz,ny,nx=np.shape(temp_pi)
  

   # get pliocene data
    f=MFDataset(fileplio)
    atemp=f.variables['temp'][:]
  
    temp_plio=np.mean(atemp,axis=0)
    temp_plio=np.squeeze(temp_plio)
   

    # get pliocene-preindustrial data
    temp_diff=temp_plio-temp_pi
    # normalize temp_diff by surface
    temp_diff_norm=temp_diff / temp_diff[0,:,:]

    nplots=2 # 2 plots per page
    npages=np.int(np.floor(nz/nplots))
    for figure in range(0,npages):
       for i in range(0,nplots):
          titlename='tdiff at pressure '+np.str(pressure[(figure *nplots)+i])
          titlename2='normalized at pressure '+np.str(pressure[(figure *nplots)+i])
          print('nums',(figure*nplots),((figure*nplots)+i),2.*((figure*nplots)+i))
          plotdata(temp_diff[(figure*nplots)+i,:,:],2*i,lon,lat,titlename,0,5.5,0.5,0.0,'n','degC')
          plotdata(temp_diff_norm[(figure*nplots)+i,:,:],2*i+1,lon,lat,titlename2,0,2.1,0.1,0.0,'a','fraction')

       fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_temp_chg_levels/'+pliop2_expt+'-'+preind_expt+'_'+np.str(pressure[(figure *nplots)+i])+'.eps' 
       plt.savefig(fileout, bbox_inches='tight')
       plt.close()
          
   
#end def annmean


def seasmean(monthnames,seasonname,preind_expt,plio_expt,pliop2_expt,extra,HadCM3):
   

    if HadCM3 == 'y':
        filepi='/nfs/hera1/earjcti/um/'+preind_expt+'/netcdf/'+preind_expt+'a@pc'+extra+'[7-9]*'
        fileplio='/nfs/hera1/earjcti/um/'+pliop2_expt+'/netcdf/'+pliop2_expt+'a@pc'+extra+'[7-9]*'
        latname='latitude_1'
        lonname='longitude_1'
    else:
        filepi='/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/netcdf/'+preind_expt+'a@pc'+extra+'[7-9]*'
        fileplio='/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/netcdf/'+pliop2_expt+'a@pc'+extra+'[7-9]*'
        latname='latitude'
        lonname='longitude'

   
    # get preindustrial data
    tottemp=0
    for i in range(0,len(monthnames)):
        print(filepi+monthnames[i]+'.nc')
        f=MFDataset(filepi+monthnames[i]+'.nc')
        lat = f.variables[latname][:]
        lon = f.variables[lonname][:]
        pressure = f.variables['p'][:]
        atemp=f.variables['temp'][:]
        tottemp=tottemp+atemp
  
    tottemp=tottemp / len(monthnames)
    temp_pi=np.mean(tottemp,axis=0)
    temp_pi=np.squeeze(temp_pi)
    nz,ny,nx=np.shape(temp_pi)
  

   # get pliocene data
    tottemp=0
    for i in range(0,len(monthnames)):
        f=MFDataset(fileplio+monthnames[i]+'.nc')
        atemp=f.variables['temp'][:]

    tottemp=tottemp+atemp
    temp_plio=np.mean(atemp,axis=0)
    temp_plio=np.squeeze(temp_plio)
   

    # get pliocene-preindustrial data
    temp_diff=temp_plio-temp_pi
    # normalize temp_diff by surface
    temp_diff_norm=temp_diff / temp_diff[0,:,:]

    nplots=2 # 2 plots per page
    npages=np.int(np.floor(nz/nplots))
    for figure in range(0,npages):
       for i in range(0,nplots):
          titlename='tdiff '+seasonname+' at press'+np.str(pressure[(figure *nplots)+i])
          titlename2='normalized at pressure '+np.str(pressure[(figure *nplots)+i])
          print('nums',(figure*nplots),((figure*nplots)+i),2.*((figure*nplots)+i))
          plotdata(temp_diff[(figure*nplots)+i,:,:],2*i,lon,lat,titlename,0,5.5,0.5,0.0,'n','degC')
          plotdata(temp_diff_norm[(figure*nplots)+i,:,:],2*i+1,lon,lat,titlename2,0,2.1,0.1,0.0,'a','fraction')

       fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_temp_chg_levels/'+pliop2_expt+'-'+preind_expt+'_'+np.str(pressure[(figure *nplots)+i])+'_'+seasonname+'.eps' 
       plt.savefig(fileout, bbox_inches='tight')
       plt.close()
          
   
#end def seasname


################################
# main program

# annual mean
#preind_expt='xiboi'
#plio_expt='xibol'
#pliop2_expt='xibol'
#extra='y'
#HadCM3='y'


preind_expt='xkvje'
plio_expt='xkvjf'
pliop2_expt='xkvjg'
extra='n'
HadCM3='n'
#annmean(preind_expt,plio_expt,pliop2_expt,extra,HadCM3)
seasmean(['dc','ja','fb'],'djf',preind_expt,plio_expt,pliop2_expt,extra,HadCM3)
seasmean(['jn','jl','ag'],'jja',preind_expt,plio_expt,pliop2_expt,extra,HadCM3)




sys.exit(0)

####

