#!/usr/bin/env python2.7
#NAME
#    PLOT_SEAS_CYCLE.py
#PURPOSE
#    This program will plot the seasonal cycle over a given region for the
#    HadCM3 data.
#    It was originally set up for the isotopes el nino paper to look at how
#    seasonal differences had changed between the Pliocene and the preindustrial
#
# search for 'main program' to find end of functions
# Julia 2/4/2017



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
#  def seasmean

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)

   # this is good for a tropical region
    map=Basemap(llcrnrlon=90.0,urcrnrlon=300.0,llcrnrlat=-45.0,urcrnrlat=45.0,projection='cyl',resolution='c')
   # this is good for the globe
   # map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
    map.drawmapboundary

    # if it comes from pf file we want to fill the land in white
    nlons=len(lons)
    if  nlons > 120:
        map.drawmapboundary(fill_color='white')

    x, y = map(lons, lats)
    map.drawcoastlines()

    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='both')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='ar':
                    cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    print(np.shape(plotdata))
                    cs = map.contourf(x,y,plotdata,V)
                    cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename,loc='left',fontsize=25)
    cbar.ax.tick_params(labelsize=15)
    #cbar.set_label(cbarname,fontsize=15)
    cbar.ax.set_title(cbarname,fontsize=20)
#end def plotdata


#==================================================
def seasmean(m1,m2,m3,seasname,exptname,fieldname):
    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean
    #==============
    # preindustrial
    
    if fieldname=='precip':
        filetype1='a'
        filetype2='d'
        varname='precip'
    
    if fieldname =='SST':
        filetype1='o'
        filetype2='f'
        varname='temp'

    if fieldname =='d18op':
        filetype1='a'
        filetype2='d'
        varname='QCL'

    if fieldname =='d18osw':
        filetype1='o'
        filetype2='f'
        varname='otracer1'

   
    # read in data from multiple files
    datasetname='/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+filetype1+'@p'+filetype2+'[w-y]*'+m1+'.nc'
    #datasetname='/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+filetype1+'@p'+filetype2+'y[7-9]*'+m1+'.nc'
    print(datasetname)
    fa=MFDataset(datasetname)
    lat = fa.variables['latitude'][:]
    lon = fa.variables['longitude'][:]
    avar=fa.variables[varname][:]
    fa.close()


    datasetname='/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+filetype1+'@p'+filetype2+'[w-y]*'+m2+'.nc'
    #datasetname='/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+filetype1+'@p'+filetype2+'y[7-9]*'+m2+'.nc'
    fb=MFDataset(datasetname)
    bvar=fb.variables[varname][:]
    fb.close()


    datasetname='/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+filetype1+'@p'+filetype2+'[w-y]*'+m3+'.nc'
    #datasetname='/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+filetype1+'@p'+filetype2+'y[7-9]*'+m3+'.nc'
    fc=MFDataset(datasetname)
    cvar=fc.variables[varname][:]
    fc.close()

    avar=np.squeeze(avar)
    bvar=np.squeeze(bvar)
    cvar=np.squeeze(cvar) 

    if fieldname=='d18op' or fieldname =='d18osw':
        ntimes,nz,ny,nx=np.shape(avar)
        ntimesb,nz,ny,nx=np.shape(bvar)
        ntimesc,nz,ny,nx=np.shape(cvar)
    else:
        ntimes,ny,nx=np.shape(avar)
        ntimesb,ny,nx=np.shape(bvar)
        ntimesc,ny,nx=np.shape(cvar)

    # abort if we are getting a different number of files from each month
    if ntimes != ntimesb or ntimes !=ntimesc:
        print('you are not getting the same number of times for each month')
        print(ntimes,ntimesb,ntimesc)
        sys.exit()

    
#average across the time dimension
    avar_avg=np.mean(avar,axis=0)
    bvar_avg=np.mean(bvar,axis=0)
    cvar_avg=np.mean(cvar,axis=0)
    
    seasvar=np.mean((avar_avg,bvar_avg,cvar_avg),axis=0)
    if fieldname=='precip':
        seasvar=seasvar * 60. * 60. * 30. * 24.

    if fieldname=='d18op':
        temporarydata=seasvar
        seasvar=0
        tot16o=temporarydata[0,:,:]+temporarydata[3,:,:]+ \
            temporarydata[6,:,:]+temporarydata[9,:,:]
        tot18o=temporarydata[1,:,:]+temporarydata[4,:,:]+ \
            temporarydata[7,:,:]+temporarydata[10,:,:]
        seasvar=((tot18o/tot16o)-2005.2E-6)/2005.2E-9
    
    if fieldname=='d18osw': #we just need level 0
        temporarydata=seasvar
        seasvar=0
        seasvar=(temporarydata[0,:,:]-2005.2E-6)/2005.2E-9
    
    allretdata=[lon,lat,seasvar]
    print('size of allretdata',np.shape(allretdata))
    return allretdata


#end def seasmean

################################
# main program

precipanom='y'
tempanom='y'
d18opanom='y'
d18oswanom='y'


if precipanom=='y':
    # figure 1d PI panom
    #djf mean for xiboi

    retdata=seasmean('dc','ja','fb','djf','xiboi','precip')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_djf_mean=retdata[2]

    #jja mean for xiboi
    retdata=seasmean('jn','jl','ag','jja','xiboi','precip')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_jja_mean=retdata[2]

    # djf - jja  mean for xiboi
    djf_jja_xiboi_precip=xiboi_djf_mean-xiboi_jja_mean
    plotdata(djf_jja_xiboi_precip,0,xiboi_lon,xiboi_lat,'d) Preindustrial DJF-JJA precipitation',-400,405,5,0.0,'a','mm/month')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1d_precip_pi_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    ##########################
    # figure 1e Pliocene panom
    #djf mean for xibol
    retdata=seasmean('dc','ja','fb','djf','xibol','precip')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_djf_mean=retdata[2]

    #jja mean for xibol
    retdata=seasmean('jn','jl','ag','jja','xibol','precip')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_jja_mean=retdata[2]

    # djf - jja  mean for xibol
    djf_jja_xibol_precip=xibol_djf_mean-xibol_jja_mean
    plotdata(djf_jja_xibol_precip,0,xibol_lon,xibol_lat,'e) mPWP DJF-JJA precipitation',-400,405,5,0.0,'a','mm/month')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1e_precip_plio_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    # djf - jja  mean for xibol-xiboi
    seas_precip_anom_xibol_xiboi=djf_jja_xibol_precip - djf_jja_xiboi_precip
    plotdata(seas_precip_anom_xibol_xiboi,0,xibol_lon,xibol_lat,'f) mPWP-PI,  DJF-JJA precipitation',-100,120,20,0.0,'a','mm/month')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1f_precip_plio_pi_djf_jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()


if tempanom=='y':
    ##########################
    # temperature anomaly fig 1a
    #djf mean for xiboi

    retdata=seasmean('dc','ja','fb','djf','xiboi','SST')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_djf_mean=retdata[2]
    
    #jja mean for xiboi
    retdata=seasmean('jn','jl','ag','jja','xiboi','SST')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_jja_mean=retdata[2]
    
    # djf - jja  mean for xiboi
    degC=u'\N{DEGREE SIGN}'+'C'
    djf_jja_xiboi_temp=xiboi_djf_mean-xiboi_jja_mean
    plotdata(djf_jja_xiboi_temp,0,xiboi_lon,xiboi_lat,'a) Preindustrial DJF-JJA SST',-15,15.5,0.5,0,'ar',degC)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1a_SST_pi_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    ##########################
    # figure 1b Pliocene Tanom
    #djf mean for xibol
    retdata=seasmean('dc','ja','fb','djf','xibol','SST')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_djf_mean=retdata[2]
    
    #jja mean for xibol
    retdata=seasmean('jn','jl','ag','jja','xibol','SST')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_jja_mean=retdata[2]

    # djf - jja  mean for xibol
    djf_jja_xibol_temp=xibol_djf_mean-xibol_jja_mean
    plotdata(djf_jja_xibol_temp,0,xibol_lon,xibol_lat,'b) mPWP DJF-JJA SST',-15,15.5,0.5,0,'ar',degC)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1b_SST_plio_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    # djf - jja  mean for xibol-xiboi
    seas_temp_anom_xibol_xiboi=djf_jja_xibol_temp - djf_jja_xiboi_temp
    plotdata(seas_temp_anom_xibol_xiboi,0,xibol_lon,xibol_lat,'c) mPWP-PI,  DJF-JJA SST',-1.0,1.2,0.2,0,'ar',degC)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1c_SST_plio_pi_djf_jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()



if d18opanom=='y':
    ##########################
    # d18op anomaly fig 1g
    #djf mean for xiboi

    retdata=seasmean('dc','ja','fb','djf','xiboi','d18op')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_djf_mean=retdata[2]
    
    #jja mean for xiboi
    retdata=seasmean('jn','jl','ag','jja','xiboi','d18op')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_jja_mean=retdata[2]
    
    # djf - jja  mean for xiboi
    djf_jja_xiboi_d18op=xiboi_djf_mean-xiboi_jja_mean
    plotdata(djf_jja_xiboi_d18op,0,xiboi_lon,xiboi_lat,'g) Preindustrial DJF-JJA d18Op',-10,10.5,0.5,0,'ar',u'\u2030')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1g_d18op_pi_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    ##########################
    # figure 1h Pliocene d18opanom
    #djf mean for xibol
    retdata=seasmean('dc','ja','fb','djf','xibol','d18op')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_djf_mean=retdata[2]
    
    #jja mean for xibol
    retdata=seasmean('jn','jl','ag','jja','xibol','d18op')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_jja_mean=retdata[2]

    # djf - jja  mean for xibol
    djf_jja_xibol_d18op=xibol_djf_mean-xibol_jja_mean
    plotdata(djf_jja_xibol_d18op,0,xibol_lon,xibol_lat,'h) mPWP DJF-JJA d18Op',-10,10.5,0.5,0,'ar',u'\u2030')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1h_d18op_plio_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    # djf - jja  mean for xibol-xiboi
    seas_d18op_anom_xibol_xiboi=djf_jja_xibol_d18op - djf_jja_xiboi_d18op
    plotdata(seas_d18op_anom_xibol_xiboi,0,xibol_lon,xibol_lat,'i) mPWP-PI,  DJF-JJA d18Op',-3.0,3.5,0.5,0,'ar',u'\u2030')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1i_d18op_plio_pi_djf_jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()




if d18oswanom=='y':
    ##########################
    # d18osw anomaly fig 1j
    #djf mean for xiboi

    retdata=seasmean('dc','ja','fb','djf','xiboi','d18osw')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_djf_mean=retdata[2]
    
    #jja mean for xiboi
    retdata=seasmean('jn','jl','ag','jja','xiboi','d18osw')
    xiboi_lon=retdata[0]
    xiboi_lat=retdata[1]
    xiboi_jja_mean=retdata[2]
    
    # djf - jja  mean for xiboi
    djf_jja_xiboi_d18osw=xiboi_djf_mean-xiboi_jja_mean
    plotdata(djf_jja_xiboi_d18osw,0,xiboi_lon,xiboi_lat,'j) Preindustrial DJF-JJA d18osw',-0.15,0.18,0.03,0,'ar',u'\u2030')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1j_d18osw_pi_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    ##########################
    # figure 1k Pliocene d18oswanom
    #djf mean for xibol
    retdata=seasmean('dc','ja','fb','djf','xibol','d18osw')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_djf_mean=retdata[2]
    
    #jja mean for xibol
    retdata=seasmean('jn','jl','ag','jja','xibol','d18osw')
    xibol_lon=retdata[0]
    xibol_lat=retdata[1]
    xibol_jja_mean=retdata[2]

    # djf - jja  mean for xibol
    djf_jja_xibol_d18osw=xibol_djf_mean-xibol_jja_mean
    plotdata(djf_jja_xibol_d18osw,0,xibol_lon,xibol_lat,'k) mPWP DJF-JJA d18Osw',-0.15,0.18,0.03,0,'ar',u'\u2030')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1k_d18osw_plio_djf-jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()

    # djf - jja  mean for xibol-xiboi
    seas_d18osw_anom_xibol_xiboi=djf_jja_xibol_d18osw - djf_jja_xiboi_d18osw
    plotdata(seas_d18osw_anom_xibol_xiboi,0,xibol_lon,xibol_lat,'l) mPWP-PI,  DJF-JJA d18osw',-0.15,0.18,0.03,0,'ar',u'\u2030')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/seas_cycle/fig_s1l_d18osw_plio_pi_djf_jja.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()








sys.exit(0)

####

