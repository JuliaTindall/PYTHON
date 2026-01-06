#!/usr/bin/env python2.7
#NAME
#    PLOT_d18o_timeslice_seasonal
#PURPOSE
#    This program will plot d18o diagnostics from timeslice experiments
#    it will plot the difference between two given experiments
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



# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,latmin,latmax,lonmin,lonmax):
    lons, lats = np.meshgrid(lon,lat)
    if fileno != 99:
        plt.subplot(2,2,fileno+1)

 
    map=Basemap(llcrnrlon=lonmin,urcrnrlon=lonmax,llcrnrlat=latmin,urcrnrlat=latmax,projection='cyl',resolution='c')
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary
    map.drawcountries()

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
                    print(V)
                    cs = map.contourf(x,y,plotdata,V,norm=mp.colors.LogNorm(vmin=0,vmax=32),cmap='Reds')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    cs = map.contourf(x,y,plotdata,V,extend='both',cmap='rainbow')
                    cbar = plt.colorbar(cs,orientation="horizontal")


    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        #cbar.set_label(cbarname,labelpad=-70,size=15)
        cbar.set_label(cbarname,size=15)
        cbar.ax.tick_params(labelsize=15)
        plt.title(titlename,loc='left',fontsize=15)
   


#end def plotdata


def get_init_data(exptname):
    
# timeslices are xiboi=preindustrial, xibol=3205 - km5c'
#xiboi/l use extra=y and dec=7-9

#  xjplc=3205-km5c, xjlpd=3060 (K1), xjple=2950 (G17), xjplf=3155 (KM3)
                #for xjplc-xjplf extra is 6 and decs are from 7-9
# continued timeslices (with smc output are) xogzs=PI, xogzt=3205-km5c xogzu=3060-k1,
#                                            xogzv=g17-2950,xogzw=km3-3155
                #for xogzt-xogzw suggest using extra 7 and dec from 1-3
                #for xogzs suggest using extra z and dec 2-4
   

    # only use this to define pliocene timeslices pi timeslices are defined elsewhere
    timeplio = {
                "xjplc" : "KM5c",
                "xjpld" : "K1",
                "xjple" : "G17",
                "xjplf" : "KM3",
                "xogzt" : "KM5c",
                "xogzu" : "K1",
                "xogzv" : "G17",
                "xogzw" : "KM3",
                }

    #print('exptname',exptname,exptname[0:3])

    if exptname[0:4]=='xibo':
        extra='y'
        startdec='7'
        enddec='9'
    
    if exptname[0:4]=='xjpl':
        extra='6'
        startdec='7'
        enddec='9'
    
    if exptname[0:4]=='xogz':
        if exptname[4]=='s':
            extra='z'
            startdec='0'
            enddec='3'
        else:
            extra='7'
            startdec='0'
            enddec='3'


    if exptname=='xiboi' or exptname=='xogzs':
        timeslice='PI'
        icevolcorr=0.0
    else:
        timeslice=timeplio.get(exptname)
        icevolcorr=-0.3
   
    #print(extra)
    
    return extra,startdec,enddec,timeslice,icevolcorr



def d18o_precip(monthnames,expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2,seasonname):

    atemp_all=0
    for mon in range(0,len(monthnames)):
        # read in data from expt1 files
        filenames='/nfs/hera1/earjcti/um/'+expt1+'/netcdf/'+expt1+'a@pd'+extra1+'['+startdec1+'-'+enddec1+']*'+monthnames[mon]+'.nc'
        print(filenames)
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        atemp=f.variables['QCL'][:] # get d18o
        atemp_all=atemp_all+atemp
        f.close()
    atemp_all=atemp_all / len(monthnames)
    atemp_avg=np.mean(atemp,axis=0)  # average by time
    tot_18o=atemp_avg[1,:,:]+atemp_avg[4,:,:]+atemp_avg[7,:,:]+atemp_avg[10,:,:]
    tot_16o=atemp_avg[0,:,:]+atemp_avg[3,:,:]+atemp_avg[6,:,:]+atemp_avg[9,:,:]
    

    # check if tot_16o is zero at poles.  If so set it to vsmow
    for j in range(0,len(lat)):
        for i in range(0,len(lon)):
            if tot_16o[j,i]==0:
                if lat[j]==90. or lat[j]==-90.:
                    tot_16o[j,i]=1.0
                    tot_18o[j,i]=2005.2E-6
  
   
    d18o_expt1=((tot_18o/tot_16o)-2005.2E-6)/2005.2E-9
    d18o_expt1,lon1 = shiftgrid(180.,d18o_expt1,lon,start=False)


 

    latmin=0.
    #latmax=40.
    latmax=60.
    lonmin=60.
    lonmax=140.


    titlename='d18O_p ('+timeslice1+')'
    plotdata(d18o_expt1,99,lon1,lat,titlename,-20,2.0,1.0,0.0,'n','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()


    atemp_all=0
    for mon in range(0,len(monthnames)):
        # read in data from expt2 files
        filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'a@pd'+extra2+'['+startdec2+'-'+enddec2+']*'+monthnames[mon]+'.nc'
        print(filenames,monthnames[mon])
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        atemp=f.variables['QCL'][:] # get d18o
        atemp_all=atemp_all+atemp
        f.close()
    atemp_all=atemp_all / len(monthnames)
    atemp_avg=np.mean(atemp,axis=0)  # average by time
    tot_18o=atemp_avg[1,:,:]+atemp_avg[4,:,:]+atemp_avg[7,:,:]+atemp_avg[10,:,:]
    tot_16o=atemp_avg[0,:,:]+atemp_avg[3,:,:]+atemp_avg[6,:,:]+atemp_avg[9,:,:]
    
    d18o_expt2=((tot_18o/tot_16o)-2005.2E-6)/2005.2E-9
    d18o_expt2,lon = shiftgrid(180.,d18o_expt2,lon,start=False)

    # plot out timeslice 2.
    titlename='d18O_p ('+timeslice2+')'
    plotdata(d18o_expt2,99,lon1,lat,titlename,-20,2.0,1.0,0.0,'n','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt2+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt2+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
    titlename='d18O_p ('+timeslice2+'-'+timeslice1+')'
    plotdata(d18o_expt2-d18o_expt1,99,lon1,lat,titlename,-5,5.5,0.5,0.0,'a','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt2+'-'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt2+'-'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
    
 

#end def d18o_precip



def d18o_runoff(monthnames,expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2,seasonname):

    totrunoff=0
    for mon in range(0,len(monthnames)):                    
        # read in data from expt1 files
        filenames='/nfs/hera1/earjcti/um/'+expt1+'/netcdf/'+expt1+'a@pd'+extra1+'['+startdec1+'-'+enddec1+']*'+monthnames[mon]+'.nc'
        print(filenames)
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        #sruoff_18o=f.variables['slowrunoff_1'][:] 
        #fruoff_18o=f.variables['fastrunoff_1'][:]
        sruoff_16o=f.variables['slowrunoff'][:] 
        fruoff_16o=f.variables['fastrunoff'][:]
        f.close()
        totrunoff=sruoff_16o+fruoff_16o+totrunoff

    totrunoff=totrunoff/len(monthnames)
           
    #runoff_18o=np.mean(sruoff_18o,axis=0)+np.mean(fruoff_18o,axis=0)
    runoff_16o=np.mean(totrunoff,axis=0)

    # calculate runoff in mm/day (currently mm/30minutes)    
    runoff_16o_expt1=np.squeeze(runoff_16o)*2. *24.
    runoff_16o_expt1,lon1 = shiftgrid(180.,runoff_16o_expt1,lon,start=False)
    #d18o_expt1=((runoff_18o/runoff_16o)-2005.2E-6)/2005.2E-9
    #d18o_expt1=np.squeeze(d18o_expt1)
    #d18o_expt1,lon1 = shiftgrid(180.,d18o_expt1,lon,start=False)

    latmin=0.
    latmax=40.
    lonmin=50.
    lonmax=130.


    #titlename='d18O_runoff ('+timeslice1+')'
    #plotdata(d18o_expt1,99,lon1,lat,titlename,-20,2.1,0.1,0.0,'n','permille',latmin,latmax,lonmin,lonmax)
    #fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_runoff_'+expt1+'.eps'
    #plt.savefig(fileout,bbox_inches='tight')
    #fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_runoff_'+expt1+'.png'
    #plt.savefig(fileout,bbox_inches='tight')
    #plt.close()


  
    # read in data from expt2 files
    totrunoff=0
    for mon in range(0,len(monthnames)):                    
        filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'a@pd'+extra2+'['+startdec2+'-'+enddec2+']*'+monthnames[mon]+'.nc'
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        #sruoff_18o=f.variables['slowrunoff_1'][:] 
        #fruoff_18o=f.variables['fastrunoff_1'][:]
        sruoff_16o=f.variables['slowrunoff'][:] 
        fruoff_16o=f.variables['fastrunoff'][:]
        f.close()
        totrunoff=sruoff_16o+fruoff_16o+totrunoff

    totrunoff=totrunoff/len(monthnames)
           
    #runoff_18o=np.mean(sruoff_18o,axis=0)+np.mean(fruoff_18o,axis=0)
    runoff_16o=np.mean(totrunoff,axis=0)
    
    # calculate runoff in mm/day (currently mm/30minutes)
    runoff_16o_expt2=np.squeeze(runoff_16o)*2. *24.
    runoff_16o_expt2,lon2 = shiftgrid(180.,runoff_16o_expt2,lon,start=False)
    #d18o_expt2=((runoff_18o/runoff_16o)-2005.2E-6)/2005.2E-9
    #d18o_expt2=np.squeeze(d18o_expt2)
    #d18o_expt2,lon = shiftgrid(180.,d18o_expt2,lon,start=False)




    latmin=0.
    latmax=40.
    lonmin=50.
    lonmax=130.

    #latmin=-90
    #latmax=90
    #lonmin=-180.
    #lonmax=180.
    
    # plot d18o runoff
    #titlename='d18O_runoff ('+timeslice2+'-'+timeslice1+')'
    #plotdata(d18o_expt2-d18o_expt1,99,lon1,lat,titlename,-10,11,1.0,0.0,'a','permille',latmin,latmax,lonmin,lonmax)
    #fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18orunoff_'+expt2+'-'+expt1+'.eps'
    #plt.savefig(fileout,bbox_inches='tight')
    #fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18orunoff_'+expt2+'-'+expt1+'.png'
    #plt.savefig(fileout,bbox_inches='tight')
    #plt.close()

    # plot runoff amount
    titlename='Total Runoff ('+timeslice2+'-'+timeslice1+')'
    plotdata(runoff_16o_expt2-runoff_16o_expt1,99,lon1,lat,titlename,-5,5.5,0.5,0.0,'a','mm/day',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/runoff_'+expt2+'-'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/runoff_'+expt2+'-'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
  
    #titlename='d18O_runoff ('+timeslice2+'-'+timeslice1+')'  
    #plotdata(d18o_expt2-d18o_expt1,0,lon1,lat,titlename,-10,11,1.0,0.0,'a','permille',latmin,latmax,lonmin,lonmax)
    #titlename='Total Runoff ('+timeslice2+'-'+timeslice1+')'
    #plotdata(runoff_16o_expt2-runoff_16o_expt1,1,lon1,lat,titlename,-0.5,0.6,0.1,0.0,'a','mm/day',latmin,latmax,lonmin,lonmax)
    #titlename='d18o Runoff ('+timeslice1+')'
    #plotdata(d18o_expt1,2,lon1,lat,titlename,-30.0,0.0,1.0,0.0,'n','mm/day',latmin,latmax,lonmin,lonmax)
    #titlename='Total Runoff ('+timeslice1+')'
    #plotdata(runoff_16o_expt1,3,lon1,lat,titlename,0,4.1,0.1,0.0,'n','mm/day',latmin,latmax,lonmin,lonmax)
   

#end def d18o_runoff
    
    
    
def d18o_smc(monthnames,expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2,seasonname):

    # read data from experiment 1 files
    for mon in range(0,len(monthnames)):                    
        # read in data from expt1 files
        filenames='/nfs/hera1/earjcti/um/'+expt1+'/netcdf/'+expt1+'a@pd'+extra1+'['+startdec1+'-'+enddec1+']*'+monthnames[mon]+'.nc'
        print(filenames)
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        smc_18o=f.variables['sm_2'][:] 
        smc_16o=f.variables['sm'][:] 
        f.close()
     
    smc_16o_1=np.squeeze(np.mean(smc_16o,axis=0)) # get smc
    
    d18o_smc_1=((np.mean(smc_18o,axis=0)/np.mean(smc_16o,axis=0))-2005.2E-6)/2005.2E-9
    d18o_smc_1=np.squeeze(d18o_smc_1)
    # mask out where smc is low in preindustrial (ie smc_16o < 5)
    d18o_smc_1_mask=np.ma.masked_where(smc_16o_1 < 1.0, d18o_smc_1)
  
    # read in data from expt2 file
    for mon in range(0,len(monthnames)):                    
        filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'a@pd'+extra2+'['+startdec2+'-'+enddec2+']*'+monthnames[mon]+'.nc'
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        smc_18o=f.variables['sm_2'][:] 
        smc_16o=f.variables['sm'][:] 
        f.close()
    
    d18o_smc_2=((np.mean(smc_18o,axis=0)/np.mean(smc_16o,axis=0))-2005.2E-6)/2005.2E-9
    d18o_smc_2=np.squeeze(d18o_smc_2)
    # mask where smc is low in preindustrial
    d18o_smc_2_mask=np.ma.masked_where(smc_16o_1 < 1.0, d18o_smc_2)
  
           
   
    latmin=0.
    latmax=60.
    lonmin=60.
    lonmax=140.
 
    # plot d18o smc
   
    titlename='SMC ('+timeslice2+'-'+timeslice1+')'
    plotdata(d18o_smc_2_mask-d18o_smc_1_mask,99,lon,lat,titlename,-5,5.5,0.5,0.0,'a','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/smc_'+expt2+'-'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/smc_'+expt2+'-'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
   
    # plot smc16o in control so that we can see 'dry'areas'
    #titlename='SMC ('+timeslice2+'-'+timeslice1+')'
    #plotdata(smc_16o_1,99,lon,lat,titlename,-0,10.5,0.5,0.0,'n','kg/m2',latmin,latmax,lonmin,lonmax)
    #plt.show()
    #sys.exit(0)
   

#end def d18o_smc



def d18o_ocean(monthnames,expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2,seasname):

    latmin=0.
    latmax=40.
    lonmin=50.
    lonmax=130.

    icevolcorr_diff=icevolcorr2-icevolcorr1
   
    # read in data from expt1 files
    atemp_tot=0
    for mon in range(0,len(monthnames)):
        filenames='/nfs/hera1/earjcti/um/'+expt1+'/netcdf/'+expt1+'o@pf'+extra1+'['+startdec1+'-'+enddec1+']*'+monthnames[mon]+'.nc'
        print(filenames)
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        atemp=f.variables['otracer1'][:] # get d18o
        f.close()
        atemp_tot=atemp_tot+atemp
    atemp_tot=atemp_tot / len(monthnames)
    atemp_avg=np.mean(atemp_tot,axis=0)  # average by time
    atemp_top=atemp_avg[0,:,:]
   
    d18o_expt1=((atemp_top - 2005.2E-6)/2005.2E-9)+icevolcorr1
    d18o_expt1,lon1 = shiftgrid(180.,d18o_expt1,lon,start=False)

    
  
    # read in data from expt2 files

    atemp_tot=0.
    for mon in range(0,len(monthnames)):
        filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'o@pf'+extra2+'['+startdec2+'-'+enddec2+']*'+monthnames[mon]+'.nc'
        print(filenames)
        f=MFDataset(filenames)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        atemp=f.variables['otracer1'][:] # get d18o
        f.close()
        atemp_tot=atemp_tot+atemp
    atemp_tot=atemp_tot / len(monthnames)
    atemp_avg=np.mean(atemp_tot,axis=0)  # average by time
    atemp_top=atemp_avg[0,:,:]
   
    d18o_expt2=((atemp_top - 2005.2E-6)/2005.2E-9)+icevolcorr2
    d18o_expt2,lon = shiftgrid(180.,d18o_expt2,lon,start=False)

    latmin=0.
    latmax=40.
    lonmin=50.
    lonmax=130.

    titlename='d18O_sw ('+timeslice2+'-'+timeslice1+')'
    plotdata(d18o_expt2-d18o_expt1,99,lon1,lat,titlename,-0.5+icevolcorr_diff,0.55+icevolcorr_diff,0.05,0.0,'a','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'-'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'-'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
  
    titlename='d18O_sw ('+timeslice2+')'
    plotdata(d18o_expt2,0,lon1,lat,titlename,-1.0,1.1,0.1,0.0,'n','permille',latmin,latmax,lonmin,lonmax)
    titlename='d18O_sw ('+timeslice1+')'
    plotdata(d18o_expt1,1,lon1,lat,titlename,-1.0,1.1,0.1,0.0,'n','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'_and_'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'_and_'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
    
 

#end def d18o_ocean


def salinity_ocean(monthnames,expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2,seasname):

    latmin=0.
    latmax=40.
    lonmin=50.
    lonmax=130.


    latmin=-90
    latmax=90
    lonmin=-180.
    lonmax=180.


    # read in data from expt1 files
    filenames='/nfs/hera1/earjcti/um/'+expt1+'/netcdf/'+expt1+'o@pg'+extra1+'['+startdec1+'-'+enddec1+']*.nc'
    print(filenames)
    f=MFDataset(filenames)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['salinity'][:] # get d18o
    f.close()
    atemp_avg=np.mean(atemp,axis=0)  # average by time
    atemp_top=atemp_avg[0,:,:]
   
    sal_expt1=(atemp_top * 1000.)+35.0
    sal_expt1,lon1 = shiftgrid(180.,sal_expt1,lon,start=False)

    
  
    # read in data from expt2 files

    filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'o@pg'+extra2+'['+startdec2+'-'+enddec2+']*.nc'
    print(filenames)
    f=MFDataset(filenames)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['salinity'][:] # get sal
    f.close()
    atemp_avg=np.mean(atemp,axis=0)  # average by time
    atemp_top=atemp_avg[0,:,:]
   
    sal_expt2=(atemp_top *1000.)+35.0
    sal_expt2,lon = shiftgrid(180.,sal_expt2,lon,start=False)

   
    titlename='sal_sw ('+timeslice2+'-'+timeslice1+')'
    plotdata(sal_expt2-sal_expt1,99,lon1,lat,titlename,-3.0,3.0,0.1,0.0,'a','psu',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'-'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'-'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
  
    titlename='sal_sw ('+timeslice2+')'
    plotdata(sal_expt2,0,lon1,lat,titlename,30.0,40.0,1.0,0.0,'n','psu',latmin,latmax,lonmin,lonmax)
    titlename='sal_sw ('+timeslice1+')'
    plotdata(sal_expt1,1,lon1,lat,titlename,30.0,40.0,1.0,0.0,'n','psu',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'_and_'+expt1+'_'+seasonname+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'_and_'+expt1+'_'+seasonname+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
    
 

#end def salinity_ocean

################################
# main program

# get all the data that we need to run the program
figureno=0
cntlexpt='xogzs' # control experiment is normally pi or km5c
expt='xogzt'
#season=['jn','jl','ag']
#season=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
season=['ag']


extracntl,startdeccntl,enddeccntl,timeslicecntl,icevolcorrcntl=get_init_data(cntlexpt)
extra,startdec,enddec,timeslice,icevolcorr=get_init_data(expt)

# get seasonname by cobbling together the first character of each 
seasonname=''
for i in range(0,len(season)):
    seasonname=seasonname+season[i][0]
 

#######################################################################################
# choose what we want to output by uncommenting as appropriate



d18o_precip(season,cntlexpt,extracntl,startdeccntl,enddeccntl,timeslicecntl,expt,extra,startdec,enddec,timeslice,icevolcorrcntl,icevolcorr,seasonname)

#d18o_ocean(season,cntlexpt,extracntl,startdeccntl,enddeccntl,timeslicecntl,expt,extra,startd#ec,enddec,timeslice,icevolcorrcntl,icevolcorr,seasonname)

#salinity_ocean(cntlexpt,extracntl,startdeccntl,enddeccntl,timeslicecntl,expt,extra,startdec,enddec,timeslice,icevolcorrcntl,icevolcorr)

#d18o_runoff(season,cntlexpt,extracntl,startdeccntl,enddeccntl,timeslicecntl,expt,extra,start#dec,enddec,timeslice,icevolcorrcntl,icevolcorr,seasonname)

d18o_smc(season,cntlexpt,extracntl,startdeccntl,enddeccntl,timeslicecntl,expt,extra,startdec,enddec,timeslice,icevolcorrcntl,icevolcorr,seasonname)



#sys.exit(0)

####

