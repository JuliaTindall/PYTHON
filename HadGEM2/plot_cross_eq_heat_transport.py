#!/usr/bin/env python2.7
#NAME
#    PLOT_cross_eq_heat_transport.py
#PURPOSE
#    This program is wanting to find out more about the movement of the 
#    ITCZ and so initially wants to find out whether the NH is warming more
#    than the southern hemisphere, or vice versa.  It will then look at the 
#    cross equatorial heat transport
#
# search for 'main program' to find end of functions
# Julia 1/11/2018



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
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(2,2,fileno+1)

   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='l')
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

#######################################################
def interhem_tdiff(mPWP_expt,pi_expt,extra,HadCM3,regionname,lonmin,lonmax):
# this will get the interhemispheric temperature differences by month
# for a given longitudinal range

    #==============
    # get all data

    if HadCM3 == 'y':
        filestart_pi='/nfs/hera1/earjcti/um/netcdf/'+pi_expt+'_netcdf/'+pi_expt+'a@pd'+extra+'[7-9]*'
        filestart_plio='/nfs/hera1/earjcti/um/netcdf/'+mPWP_expt+'_netcdf/'+mPWP_expt+'a@pd'+extra+'[7-9]*'
        fieldname='temp'
        suffix='.nc'
    else:
        filestart_pi='/nfs/hera1/earjcti/um/HadGEM_data/'+pi_expt+'/temp_data/'+pi_expt+'a@pd'+extra+'[7-9]*'
        filestart_plio='/nfs/hera1/earjcti/um/HadGEM_data/'+mPWP_expt+'/temp_data/'+mPWP_expt+'a@pd'+extra+'[7-9]*'
        fieldname='temp_1'
        suffix='_temp.nc'

    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    nmon=len(monthnames)

    for month in range(0,nmon):
        # get preindustrial data (average monthly temperature over globe)
        f=MFDataset(filestart_pi+monthnames[month]+suffix)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        atemp=f.variables[fieldname][:]
        atemp=np.squeeze(atemp)
        temp_avg=np.mean(atemp,axis=0)
        ny,nx=np.shape(temp_avg)

        # shiftgrid if required
        if lonmin < 0 and lonmax > 0:
            temp_avg,lon=shiftgrid(180.,temp_avg,lon,start=False)

        if month==0:
            # find the number of longitudes
            nx_redu=0
            lonstart=999
            for x in range(0,len(lon)):
                if lon[x] >=lonmin and lon[x] <=lonmax:
                    nx_redu=nx_redu+1
                    if lonstart > 900:
                        lonstart=x
            temp_pi=np.zeros((nmon,ny,nx_redu))
            lonredu=np.zeros(nx_redu)

            
        for x in range(0,nx_redu):
            temp_pi[month,:,x]=temp_avg[:,x+lonstart]
            lonredu[x]=lon[x+lonstart]

        if month == 0:
            print(lonmin,lonmax)
            plotdata(temp_pi[0,:,:],0,lonredu,lat,'region',270,300,5.0,0.0,'n',regionname)
            plt.show()

    
        # get mPWP data (average monthly temperature over globe)
        f=MFDataset(filestart_plio+monthnames[month]+suffix)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        atemp=f.variables[fieldname][:]
        atemp=np.squeeze(atemp)
        temp_avg=np.mean(atemp,axis=0)
        ny,nx=np.shape(temp_avg)

        # shiftgrid if required
        if lonmin < 0 and lonmax > 0:
            temp_avg,lon=shiftgrid(180.,temp_avg,lon,start=False)

        if month==0:
            temp_plio=np.zeros((nmon,ny,nx_redu))
        for x in range(0,nx_redu):
            temp_plio[month,:,x]=temp_avg[:,x+lonstart]
            lonredu[x]=lon[x+lonstart]
    

    ###########################################################
    # get the average interhemispheric temperature difference

    # create weighting array for the NH and the SH
    nh_weights=np.zeros((ny,nx_redu))
    sh_weights=np.zeros((ny,nx_redu))
   
    for i in range(0,nx_redu):
        nh_weights[:,i]=np.cos(np.deg2rad(lat))
        sh_weights[:,i]=np.cos(np.deg2rad(lat))

    for j in range(0,ny):
        if lat[j] < 0:
            nh_weights[j,:]=0.0
        if lat[j] > 0:
            sh_weights[j,:]=0.0

    nh_plio=np.zeros(nmon)
    nh_pi=np.zeros(nmon)
    sh_plio=np.zeros(nmon)
    sh_pi=np.zeros(nmon)

    for mon in range(0,nmon):
        nh_plio[mon]=np.average(temp_plio[mon,:,:],weights=nh_weights)
        nh_pi[mon]=np.average(temp_pi[mon,:,:],weights=nh_weights)
        sh_plio[mon]=np.average(temp_plio[mon,:,:],weights=sh_weights)
        sh_pi[mon]=np.average(temp_pi[mon,:,:],weights=sh_weights)
                        

    plt.subplot(3,1,1)
    plt.plot(nh_plio,label='nh_plio')
    plt.plot(sh_plio,label='sh_plio')
    plt.plot(nh_pi,label='nh_pi')
    plt.plot(sh_pi,label='sh_pi')
    plt.title(regionname)
    plt.legend()
    
    plt.subplot(3,1,2)
    plt.plot(nh_plio-sh_plio,label='NH-SH plio')
    plt.plot(nh_pi-sh_pi,label='NH-SH pi')
    plt.legend()

    plt.subplot(3,1,3)
    plt.plot((nh_plio-sh_plio)-(nh_pi-sh_pi),label='(pliocene - pi) nh-sh contrast')
    plt.legend()

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/cross_eq_heat_transport/NH-SH_'+mPWP_expt+'_'+regionname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/cross_eq_heat_transport/NH-SH_'+mPWP_expt+'_'+regionname+'.tiff' 
    plt.savefig(fileout, bbox_inches='tight')  
   
    plt.close()
    

    ###########################################################
    # get the average interhemispheric temperature difference
    # for the tropics only

    # create weighting array for the NH and the SH
    nh_weights_tropics=nh_weights
    sh_weights_tropics=sh_weights
   
    for j in range(0,ny):
        if lat[j] < -45:
            sh_weights[j,:]=0.0
        if lat[j] > 45:
            nh_weights[j,:]=0.0


    #plotdata(sh_weights_tropics,0,lon,lat,'weights',0,1,0.1,0.0,'n','weights')
    #plt.show()

    nh_tropics_plio=np.zeros(nmon)
    nh_tropics_pi=np.zeros(nmon)
    sh_tropics_plio=np.zeros(nmon)
    sh_tropics_pi=np.zeros(nmon)

    for mon in range(0,nmon):
        nh_tropics_plio[mon]=np.average(temp_plio[mon,:,:],weights=nh_weights_tropics)
        nh_tropics_pi[mon]=np.average(temp_pi[mon,:,:],weights=nh_weights_tropics)
        sh_tropics_plio[mon]=np.average(temp_plio[mon,:,:],weights=sh_weights_tropics)
        sh_tropics_pi[mon]=np.average(temp_pi[mon,:,:],weights=sh_weights_tropics)
                        

    plt.subplot(3,1,1)
    plt.plot(nh_tropics_plio,label='nh_tropics_plio')
    plt.plot(sh_tropics_plio,label='sh_tropics_plio')
    plt.plot(nh_tropics_pi,label='nh_tropics_pi')
    plt.plot(sh_tropics_pi,label='sh_tropics_pi')
    plt.title(regionname)
    plt.legend()
    
    plt.subplot(3,1,2)
    plt.plot(nh_tropics_plio-sh_tropics_plio,label='NH-SH tropics plio')
    plt.plot(nh_tropics_pi-sh_tropics_pi,label='NH-SH tropicspi')
    plt.legend()

    plt.subplot(3,1,3)
    plt.plot((nh_tropics_plio-sh_tropics_plio)-(nh_tropics_pi-sh_tropics_pi),label='Tropics (pliocene - pi) nh-sh contrast')
    plt.legend()

   
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/cross_eq_heat_transport/NH-SH_tropics_'+mPWP_expt+'_'+regionname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/cross_eq_heat_transport/NH-SH_tropics_'+mPWP_expt+'_'+regionname+'.tiff' 
    plt.savefig(fileout, bbox_inches='tight')  
   
    plt.close()
  
#end def interhem_tdiff



#######################################################
def heat_transport(mPWP_expt,pi_expt,extra,HadCM3,regionname,lonmin,lonmax):
# this will get the interhemispheric temperature differences by month
# for a given longitudinal range

  
#end def heat_transport


################################
# main program

# annual mean
figureno=0

HadCM3='n'
mPWP_expt='xkvjg'
pi_expt='xkvje'
extra='n'
#plt.figure(figureno)
#interhem_tdiff(mPWP_expt,pi_expt,extra,HadCM3,'globe',0,360)
interhem_tdiff(mPWP_expt,pi_expt,extra,HadCM3,'Atlantic',-60,20)
#interhem_tdiff(mPWP_expt,pi_expt,extra,HadCM3,'Indian',30,120)
#interhem_tdiff(mPWP_expt,pi_expt,extra,HadCM3,'W_Pacific',150,200)
#interhem_tdiff(mPWP_expt,pi_expt,extra,HadCM3,'E_Pacific',240,270)

#figureno=figureno+1


sys.exit(0)

####

