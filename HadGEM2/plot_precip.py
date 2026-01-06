#!/usr/bin/env python2.7
#NAME
#    PLOT_SURFPRECIP
#PURPOSE
#    This program will plot the precipitation (annual and seasonal) and
#    the precipitation anomaly (annual and seasonal) for the final 30 years
#    of the HadGEM2 simulations
#
# search for 'main program' to find end of functions
# Julia 22/11/2016
# Julia 26/08/2018 added a subprogram to focus on plotting precipitation
#                  in ITCZ region



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans


#functions are:
#  def plotdata
#  def plotmap_itcz # a copy of plotdata for doing itcz
#  def annmean
#  def seasmean
#  def focusitcz

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,mask_ind):
    lons, lats = np.meshgrid(lon,lat)
    if fileno !=99:
        plt.subplot(2,2,fileno+1)


    if mask_ind == 'l': # ;land mask
        plotnew=maskoceans(lons,lats,plotdata)
        plotdata=plotnew
        if cbarname=='mm/day':
            minval=minval/2.
            maxval=maxval/2.
            valinc=valinc/2.

    if mask_ind =='t': # tropics mask
        northlat=30.0
        southlat=-30.0
    else:
        northlat=90.0
        southlat=-90.0

   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=southlat,urcrnrlat=northlat,projection='cyl',resolution='c')
    x, y = map(lons, lats)
    map.drawcoastlines()

    plotdata2=plotdata
    #plotdata=maskoceans(x,y,plotdata)
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu',extend='both')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='both')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='ra':
                    cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
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
   

    plotdata=plotdata2

    if mask_ind == 'l':
        map.drawmapboundary(fill_color='white')
    else:
        map.drawmapboundary

#end def plotdata


def plotmap_itcz(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,mask_ind,ygrid,xgrid,yuse,xuse,yspan,xspan):


    if fileno > 3:
        print('ERROR NOT ENOUGH SPACE ON PAGE ',fileno)
        sys.exit()
    #plt.subplot2grid((10,12),(fileno*2,0),colspan=9,rowspan=2)

    print(ygrid,xgrid,yuse,xuse,yspan,xspan)
    plt.subplot2grid((ygrid,xgrid),(yuse,xuse),colspan=yspan,rowspan=xspan)

    lons, lats = np.meshgrid(lon,lat)
   
    if mask_ind =='t': # tropics mask
        northlat=30.0
        southlat=-30.0
    elif mask_ind == 'talt': # alternative tropical region:
        northlat=45.0
        southlat=-30.0
    else:
        northlat=90.0
        southlat=-90.0

    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=southlat,urcrnrlat=northlat,projection='cyl',resolution='c',fix_aspect=False)
    x, y = map(lons, lats)
    map.drawcoastlines()

    plotdata2=plotdata


    # plot limits
    if V == 0:
        V=np.arange(minval,maxval,valinc)

    # plot map
        
    if uselog == 'n':
        cs = map.contourf(x,y,plotdata,V,cmap='YlGnBu',extend='both')
    else:
        cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
   
    #parallels=np.arange(-90.,90.,15.)
    parallels=(-15,0,15,30)
    map.drawparallels(parallels,labels=[False,True,False,False]) # labels right
    meridians=np.arange(-180.,180.,60.)
    if fileno !=3:
        map.drawmeridians(meridians,labels=[False,False,False,False]) # nolabels
    else:
        map.drawmeridians(meridians,labels=[False,False,False,True]) # labels bottom
   
    fontsize=10
    plt.text(-180.0-6,northlat-fontsize-1,titlename,fontsize=fontsize,ha="right",bbox=dict(boxstyle="square,pad=0.1",color="white"))
 
    # colorbar
    if fileno==0:
        plt.subplot2grid((10,12),(9,0),colspan=9,rowspan=1)
        plt.gca().set_visible(False)
        #cax=plt.axes([0.85,0.1,0.075,0.8])
        #cbar = plt.colorbar(cax=cax)
        #cbar = plt.colorbar(cs,orientation="horizontal",cax=cax)
        cbar = plt.colorbar(cs,orientation="horizontal",fraction=1.0)         
        cbar.set_label(cbarname)
        cbar.ax.tick_params(labelsize=10)
        
   

    plotdata=plotdata2
    
    # plot map boundary
    map.drawmapboundary

#end def plotmap_itcz

def annmean(figureno,preind_expt,plio_expt,pliop2_expt,extra,mask_ind):
    #==============
    # preindustrial


    # read in data from multiple files
    f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/precip_data/'+preind_expt+'a@pd'+extra+'[5-9]*.nc')
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    aprecip=f.variables['precip_1'][:]
    aprecip=np.squeeze(aprecip)
    ntimes,ny,nx=np.shape(aprecip)
    print(ntimes,ny,nx)
    f.close()

#average across the time dimension
    pi_precip_ann=np.mean(aprecip,axis=0)
    print('new shape',np.shape(pi_precip_ann))
    
    pi_precip_ann=pi_precip_ann * 60. * 60. * 24. * 30.
    
    plt.figure(0)
    lonprecip=lon
    pi_precip_ann,lon = shiftgrid(180.,pi_precip_ann,lon,start=False)
    
    plotdata(pi_precip_ann,0,lon,lat,'PI-Ann_HadGEM2',0,275,25.0,0.0,'n','mm/month',mask_ind)


     #==============
     # Pliocene


    f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+plio_expt+'/precip_data/'+plio_expt+'a@pd'+extra+'[5-9]*.nc')
    aprecip=f.variables['precip_1'][:]
    aprecip=np.squeeze(aprecip)
    ntimes,ny,nx=np.shape(aprecip)
    print(ntimes,ny,nx)
    f.close()


    plio_precip_ann=np.mean(aprecip,axis=0)
    plio_precip_ann=plio_precip_ann * 60. * 60. * 24. * 30.

    lon=lonprecip
    plio_precip_ann,lon = shiftgrid(180.,plio_precip_ann,lon,start=False)

    plotdata(plio_precip_ann,1,lon,lat,'Plio-PAnn_HG2',0,275,25,0.0,'n','mm/month',mask_ind)


     #==============
     # Pliocene+2


     # read in data from multiple files
    f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/precip_data/'+pliop2_expt+'a@pd'+extra+'[5-9]*.nc')
    aprecip=f.variables['precip_1'][:]
    aprecip=np.squeeze(aprecip)
    ntimes,ny,nx=np.shape(aprecip)
    print(ntimes,ny,nx)
    f.close()

    
    #average across the time dimension
    plio_precipp2_ann=np.mean(aprecip,axis=0)
    plio_precipp2_ann=plio_precipp2_ann * 60. * 60. * 24. * 30.
    lon=lonprecip
    plio_precipp2_ann,lon = shiftgrid(180.,plio_precipp2_ann,lon,start=False)




    # Pliocene - preindustrial

    plio_anom=plio_precip_ann-pi_precip_ann

    print('plio_precip_ann',plio_precip_ann[20,:])
    print('pi_precip_ann',pi_precip_ann[20,:])
    print('plio_anom',plio_anom[20,:])


    V=[-128,-64,-32,-16,-8,-4,0,4,8,16,32,64,128]
    plotdata(plio_anom,2,lon,lat,'Plio - PI Panom_HG2',0,275,25,V,'la','mm/month',mask_ind)
    
    # Pliocene+2 - preindustrial

    pliop2_anom=plio_precipp2_ann-pi_precip_ann
    V=[-128,-64,-32,-16,-8,-4,0,4,8,16,32,64,128]

    plotdata(pliop2_anom,3,lon,lat,'PlioP2 - PI Panom_HG2',0,275,25,V,'la','mm/month',mask_ind)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_anom_'+pliop2_expt+'_'+plio_expt+'_'+preind_expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

    ###############################################
    # plot anomaly in mm/day for paper.
    # also calculate mean changes
    #V=[-4,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,4]
    #plotdata(pliop2_anom/30.,99,lon,lat,'b) mPWP precipitation anomaly',0,275,25,V,'la','mm/day')
    #plt.show()

    V=0
    plotdata(pliop2_anom/30.,99,lon,lat,'b) mPWP precipitation anomaly',-2.0,2.2,0.2,V,'a','mm/day',mask_ind)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_anom_only_'+pliop2_expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

  # get land mask and put on correct grid

    fm=Dataset('/nfs/hera1/earjcti/um/HadGEM_ancils/qrparm.mask.nc')
    lsmlon=fm.variables['longitude'][:]
    lsmlat=fm.variables['latitude'][:]
    lsm=fm.variables['lsm'][:]
    lsm=np.squeeze(lsm)
    lsm,lsmlon = shiftgrid(180.,lsm,lsmlon,start=False)
    fm.close()

    if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
        anom_land=pliop2_anom / (30. * lsm)
        anom_sea=pliop2_anom / (30. * np.abs(lsm-1.0))
    else:
        print('error lon/lat of land sea mask dont match')
        anom_land=pliop2_anom * lsm
        plotdata(anom_land,99,lon,lat,'a) mPWP temperature anomaly',0,10,1.0,V,'i',degC,mask_ind)
        plt.show()
        sys.exit()


    # create weighting array
    weightarr=np.zeros(np.shape(anom_sea))
    for i in range(0,len(lon)):
        weightarr[:,i]=np.cos(np.deg2rad(lat))

    print('mean anom_sea',np.average(pliop2_anom/30.,weights=weightarr * np.abs(lsm-1.0)))
    print('mean anom_land',np.average(pliop2_anom/30.,weights=weightarr*lsm))
    print('allmean',np.average(pliop2_anom/30.,weights=weightarr))



    print('non plus 2')
    print('mean anom_sea',np.average(plio_anom/30.,weights=weightarr * np.abs(lsm-1.0)))
    print('mean anom_land',np.average(plio_anom/30.,weights=weightarr*lsm))
    print('allmean',np.average(plio_anom/30.,weights=weightarr))

    

 


    


    # Pliocene - preindustrial percentage change

    plio_peranom=((plio_precip_ann-pi_precip_ann)/pi_precip_ann)*100.
    plotdata(plio_peranom,0,lon,lat,'Plio - PI Panom_HG2 %',-50,60,5,0,'a','%',mask_ind)

    plio_peranomp2=((plio_precipp2_ann-pi_precip_ann)/pi_precip_ann)*100.
    plotdata(plio_peranomp2,1,lon,lat,'Plio - PI Panom_HG2+2 %',-50,60,5,0,'a','%',mask_ind)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_anom_percent_'+pliop2_expt+'_'+plio_expt+'_'+preind_expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

    # plot percentage change on a single plot

    plotdata(plio_peranomp2,99,lon,lat,'c) mPWP precipitation anomaly',-70,78,5,0,'a','%',mask_ind)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_anom_only_pcent'+pliop2_expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()



#end def annmean


def seasmean(m1,m2,m3,figureno,seasname,preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind):

    land_ocn_ext=''
    if mask_ind == 'l':
        land_ocn_ext='_land'
    if mask_ind == 'o':
        land_ocn_ext='_ocean'

    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean
    #==============
    # preindustrial

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/precip_data/'
    fileend='_precip.nc'
    fieldreq='precip_1'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+preind_expt+'_netcdf/'
        fileend='.nc'
        fieldreq='precip'
    
    filestart=filestart+preind_expt+'a@pd'+extra+'*'

   
    # read in data from multiple files
    print(filestart+m1+fileend)
    fa=MFDataset(filestart+m1+fileend)
    fb=MFDataset(filestart+m2+fileend)
    fc=MFDataset(filestart+m3+fileend)
    lat = fa.variables['latitude'][:]
    lon = fa.variables['longitude'][:]
    aprecip=fa.variables[fieldreq][:]
    bprecip=fb.variables[fieldreq][:]
    cprecip=fc.variables[fieldreq][:]
    aprecip=np.squeeze(aprecip)
    bprecip=np.squeeze(bprecip)
    cprecip=np.squeeze(cprecip)
    ntimes,ny,nx=np.shape(aprecip)
    print(ntimes,ny,nx)

    fa.close()
    fb.close()
    fc.close()

    
#average across the time dimension
    pi_aprecip_avg=np.mean(aprecip,axis=0)
    pi_bprecip_avg=np.mean(bprecip,axis=0)
    pi_cprecip_avg=np.mean(cprecip,axis=0)
    

# find standard deviation across time dimension
    pi_aprecip_stdev=np.std(aprecip*60.*60.*24.*30.,axis=0)
    pi_bprecip_stdev=np.std(bprecip*60.*60.*24.*30.,axis=0)
    pi_cprecip_stdev=np.std(cprecip*60.*60.*24.*30.,axis=0)



    pi_seasprecip=np.mean((pi_aprecip_avg,pi_bprecip_avg,pi_cprecip_avg),axis=0)
    pi_seasprecip=pi_seasprecip * 60. * 60. * 30. * 24.
    
    
    lonprecip=lon
    pi_seasprecip,lon = shiftgrid(180.,pi_seasprecip,lon,start=False)
    
    plotdata(pi_seasprecip,0,lon,lat,'PI HadGEM2: '+seasname,0,275,25,0.0,'n','mm/month',mask_ind)
    
     #==============
     # Pliocene

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+plio_expt+'/precip_data/'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+plio_expt+'_netcdf/'

    filestart=filestart+plio_expt+'a@pd'+extra+'*'


    fa=MFDataset(filestart+m1+fileend)
    fb=MFDataset(filestart+m2+fileend)
    fc=MFDataset(filestart+m3+fileend)
    aprecip=fa.variables[fieldreq][:]
    bprecip=fb.variables[fieldreq][:]
    cprecip=fc.variables[fieldreq][:]
    aprecip=np.squeeze(aprecip)
    bprecip=np.squeeze(bprecip)
    cprecip=np.squeeze(cprecip)

    # average across the time dimension    
    plio_aprecip_avg=np.mean(aprecip,axis=0)
    plio_bprecip_avg=np.mean(bprecip,axis=0)
    plio_cprecip_avg=np.mean(cprecip,axis=0)
    
    fa.close()
    fb.close()
    fc.close()



    plio_seasprecip=np.mean((plio_aprecip_avg,plio_bprecip_avg,plio_cprecip_avg),axis=0)
    plio_seasprecip=plio_seasprecip * 60. * 60. * 30. * 24.

    lon=lonprecip
    plio_seasprecip,lon = shiftgrid(180.,plio_seasprecip,lon,start=False)
    
    



     #==============
     # Pliocene+2

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/precip_data/'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+pliop2_expt+'_netcdf/'

    filestart=filestart+pliop2_expt+'a@pd'+extra+'*'


    fa=MFDataset(filestart+m1+fileend)
    fb=MFDataset(filestart+m2+fileend)
    fc=MFDataset(filestart+m3+fileend)
    aprecip=fa.variables[fieldreq][:]
    bprecip=fb.variables[fieldreq][:]
    cprecip=fc.variables[fieldreq][:]
    aprecip=np.squeeze(aprecip)
    bprecip=np.squeeze(bprecip)
    cprecip=np.squeeze(cprecip)

    fa.close()
    fb.close()
    fc.close()

    
    pliop2_aprecip_avg=np.mean(aprecip,axis=0)
    pliop2_bprecip_avg=np.mean(bprecip,axis=0)
    pliop2_cprecip_avg=np.mean(cprecip,axis=0)


    # find standard deviation across time dimension
    pliop2_aprecip_stdev=np.std(aprecip*60.*60.*24.*30.,axis=0)
    pliop2_bprecip_stdev=np.std(bprecip*60.*60.*24.*30.,axis=0)
    pliop2_cprecip_stdev=np.std(cprecip*60.*60.*24.*30.,axis=0)

    
    pliop2_seasprecip=np.mean((pliop2_aprecip_avg,pliop2_bprecip_avg,pliop2_cprecip_avg),axis=0)
    pliop2_seasprecip=pliop2_seasprecip * 60. * 60. * 30. * 24.

    lon=lonprecip
    pliop2_seasprecip,lon = shiftgrid(180.,pliop2_seasprecip,lon,start=False)
    

    plotdata(plio_seasprecip,1,lon,lat,'Plio+2 HadGEM2: '+seasname,0,275,25,0.0,'n','mm/month',mask_ind)
 

    # Pliocene - preindustrial

    plio_anom=plio_seasprecip-pi_seasprecip
    V=[-128,-64,-32,-16,-8,-4,0,4,8,16,32,64,128]
   

    plotdata(plio_anom,2,lon,lat,'Plio - PI Panom_HG2',0,275,25,V,'la','mm/month',mask_ind)
    
    # Pliocene+2 - preindustrial

    pliop2_anom=pliop2_seasprecip-pi_seasprecip
    #V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]
    V=[-128,-64,-32,-16,-8,-4,0,4,8,16,32,64,128]
   

    plotdata(pliop2_anom,3,lon,lat,'PlioP2 - PI Panom_HG2',0,275,25,V,'la','mm/month',mask_ind)



    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_'+seasname+'anom_'+pliop2_expt+'_'+plio_expt+'_'+preind_expt+land_ocn_ext+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()


    ###############################################
    # plot anomaly in mm/day for paper.
    # also calculate mean changes

    V=0
    titlename='b) mPWP precipitation anomaly '+seasname
    plotdata(pliop2_anom/30.,99,lon,lat,titlename,-2.0,2.2,0.2,V,'a','mm/day',mask_ind)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_anom_only_'+pliop2_expt+'_'+seasname+land_ocn_ext+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()



    # Pliocene - preindustrial percentage change

    plio_peranom=((plio_seasprecip-pi_seasprecip)/pi_seasprecip)*100.
    titlename='Plio-PI precip %'+seasname
    plotdata(plio_peranom,0,lon,lat,titlename,-50,50,5,0,'a','%',mask_ind)

    plio_peranomp2=((pliop2_seasprecip-pi_seasprecip)/pi_seasprecip)*100.
    titlename='Plio+2-PI precip %'+seasname
    plotdata(plio_peranomp2,1,lon,lat,titlename,-50,50,5,0,'a','%',mask_ind)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_anom_percent_'+seasname+'_'+pliop2_expt+'_'+plio_expt+'_'+preind_expt+land_ocn_ext+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()


    # plot percentage change on a single plot
    titlename='c) mPWP precipitation anomaly'+seasname
    plotdata(plio_peranomp2,99,lon,lat,titlename,-70,78,5,0,'a','%',mask_ind)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/MAP_anom_only_pcent'+pliop2_expt+'_'+seasname+land_ocn_ext+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()





# plot standard deviation for plio for months

    plotdata(pliop2_aprecip_stdev,0,lon,lat,'Plio stdev m1: '+seasname,0,100,5.0,0.0,'n','mm/month',mask_ind)
    plotdata(pliop2_bprecip_stdev,1,lon,lat,'Plio stdev m2: '+seasname,0,100,5.0,0.0,'n','mm/month',mask_ind)
    plotdata(pliop2_cprecip_stdev,2,lon,lat,'Plio stdev m3: '+seasname,0,100,5.0,0.0,'n','mm/month',mask_ind)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/stdevP'+seasname+'_'+pliop2_expt+land_ocn_ext+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()
   
# plot difference in standard deviation for each month

    stdev_anom_1=pliop2_aprecip_stdev-pi_aprecip_stdev
    stdev_anom_2=pliop2_bprecip_stdev-pi_bprecip_stdev
    stdev_anom_3=pliop2_cprecip_stdev-pi_cprecip_stdev
    plotdata(stdev_anom_1,0,lon,lat,'Plio stdev m1: '+seasname,-20,25.0,5.0,0.0,'ra','mm/month',mask_ind)
    plotdata(stdev_anom_2,1,lon,lat,'Plio stdev m2: '+seasname,-20,25.0,5.0,0.0,'ra','mm/month',mask_ind)
    plotdata(stdev_anom_3,2,lon,lat,'Plio stdev m3: '+seasname,-20,25.0,5.0,0.0,'ra','mm/month',mask_ind)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/stdevP'+seasname+'anom'+pliop2_expt+land_ocn_ext+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()



     


#end def seasmean



def seas_snowrain(m1,m2,m3,seasname,preind_expt,plio_expt,extra,HadCM3,mask_ind):

    land_ocn_ext=''
    if mask_ind == 'l':
        land_ocn_ext='_land'
    if mask_ind == 'o':
        land_ocn_ext='_ocean'

    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean
    #==============
    # preindustrial

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/netcdf/pdfiles/'
    fileend='.nc'
    fieldrain='rain'
    fieldsnow='snow'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+preind_expt+'_netcdf/'
        fileend='.nc'
        fieldreq='precip'
    
    filestart=filestart+preind_expt+'a@pd'+extra+'*'

   
    # read in data from multiple files
    print(filestart+m1+fileend)
    fa=MFDataset(filestart+m1+fileend)
    fb=MFDataset(filestart+m2+fileend)
    fc=MFDataset(filestart+m3+fileend)
    lat = fa.variables['latitude'][:]
    lon = fa.variables['longitude'][:]
    arain=fa.variables[fieldrain][:]
    brain=fb.variables[fieldrain][:]
    crain=fc.variables[fieldrain][:]
    arain=np.squeeze(arain)
    brain=np.squeeze(brain)
    crain=np.squeeze(crain)
    ntimes,ny,nx=np.shape(arain)
    print(ntimes,ny,nx)


    asnow=fa.variables[fieldsnow][:]
    bsnow=fb.variables[fieldsnow][:]
    csnow=fc.variables[fieldsnow][:]
    asnow=np.squeeze(asnow)
    bsnow=np.squeeze(bsnow)
    csnow=np.squeeze(csnow)
   
    fa.close()
    fb.close()
    fc.close()

    
#average across the time dimension
    pi_arain_avg=np.mean(arain,axis=0)
    pi_brain_avg=np.mean(brain,axis=0)
    pi_crain_avg=np.mean(crain,axis=0)
   
    pi_asnow_avg=np.mean(asnow,axis=0)
    pi_bsnow_avg=np.mean(bsnow,axis=0)
    pi_csnow_avg=np.mean(csnow,axis=0)
    

    pi_seasrain=np.mean((pi_arain_avg,pi_brain_avg,pi_crain_avg),axis=0)
    pi_seasrain=pi_seasrain * 60. * 60. * 30. * 24.
  
    pi_seassnow=np.mean((pi_asnow_avg,pi_bsnow_avg,pi_csnow_avg),axis=0)
    pi_seassnow=pi_seassnow * 60. * 60. * 30. * 24.
    
    
    lonprecip=lon
    pi_seasrain,lon = shiftgrid(180.,pi_seasrain,lon,start=False)
    lon=lonprecip
    pi_seassnow,lon = shiftgrid(180.,pi_seassnow,lon,start=False)
    
    plotdata(pi_seasrain,0,lon,lat,'PI HadGEM2:rain '+seasname,0,50,5,0.0,'n','mm/month',mask_ind)
    
    plotdata(pi_seassnow,1,lon,lat,'PI HadGEM2:snow '+seasname,0,50,5,0.0,'n','mm/month',mask_ind)
    
     #==============
     # Pliocene

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+plio_expt+'/netcdf/pdfiles/'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+plio_expt+'_netcdf/'

    filestart=filestart+plio_expt+'a@pd'+extra+'*'

    print(filestart+m1+fileend)
    fa=MFDataset(filestart+m1+fileend)
    fb=MFDataset(filestart+m2+fileend)
    fc=MFDataset(filestart+m3+fileend)
    arain=fa.variables[fieldrain][:]
    brain=fb.variables[fieldrain][:]
    crain=fc.variables[fieldrain][:]
    arain=np.squeeze(arain)
    brain=np.squeeze(brain)
    crain=np.squeeze(crain)
    asnow=fa.variables[fieldsnow][:]
    bsnow=fb.variables[fieldsnow][:]
    csnow=fc.variables[fieldsnow][:]
    asnow=np.squeeze(asnow)
    bsnow=np.squeeze(bsnow)
    csnow=np.squeeze(csnow)

    # average across the time dimension    
    plio_arain_avg=np.mean(arain,axis=0)
    plio_brain_avg=np.mean(brain,axis=0)
    plio_crain_avg=np.mean(crain,axis=0)
   
    plio_asnow_avg=np.mean(asnow,axis=0)
    plio_bsnow_avg=np.mean(bsnow,axis=0)
    plio_csnow_avg=np.mean(csnow,axis=0)
    
    fa.close()
    fb.close()
    fc.close()



    plio_seasrain=np.mean((plio_arain_avg,plio_brain_avg,plio_crain_avg),axis=0)
    plio_seasrain=plio_seasrain * 60. * 60. * 30. * 24.
    plio_seassnow=np.mean((plio_asnow_avg,plio_bsnow_avg,plio_csnow_avg),axis=0)
    plio_seassnow=plio_seassnow * 60. * 60. * 30. * 24.

    lon=lonprecip
    plio_seasrain,lon = shiftgrid(180.,plio_seasrain,lon,start=False)
    lon=lonprecip
    plio_seassnow,lon = shiftgrid(180.,plio_seassnow,lon,start=False)
    
   
    plotdata(plio_seasrain,2,lon,lat,'Plio+2 HadGEM2: rain '+seasname,0,50,5,0.0,'n','mm/month',mask_ind)
 
    plotdata(plio_seassnow,3,lon,lat,'Plio+2 HadGEM2: snow '+seasname,0,50,5,0.0,'n','mm/month',mask_ind)
 

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/rainsnow_'+plio_expt+'_'+preind_expt+land_ocn_ext+seasname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()
    
    # Pliocene - preindustrial

    plio_rain_anom=plio_seasrain-pi_seasrain
    plio_snow_anom=plio_seassnow-pi_seassnow
   

    V=[-128,-64,-32,-16,-8,-4,0,4,8,16,32,64,128]
    V=0
   
    plotdata(plio_rain_anom,99,lon,lat,'Plio - PI Panom_HG2; rain '+seasname,-20,22,2,V,'a','mm/month',mask_ind)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/rainanom_'+plio_expt+'_'+preind_expt+land_ocn_ext+seasname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()
   
    plotdata(plio_snow_anom,99,lon,lat,'Plio - PI Panom_HG2; snow '+seasname,-20,22,2,V,'a','mm/month',mask_ind)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/snowanom_'+plio_expt+'_'+preind_expt+land_ocn_ext+seasname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()
    
    

#end def seasmean


def focusitcz(m1,m2,m3,figureno,seasname,preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind):

  
    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean
    #==============
    # preindustrial

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/precip_data/'
    fileend='_precip.nc'
    fieldreq='precip_1'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+preind_expt+'_netcdf/'
        fileend='.nc'
        fieldreq='precip'
    
    filestart=filestart+preind_expt+'a@pd'+extra+'*'

   
    # read in data from multiple files
    print(filestart+m1+fileend)
    fa=MFDataset(filestart+m1+fileend)
    fb=MFDataset(filestart+m2+fileend)
    fc=MFDataset(filestart+m3+fileend)
    lat = fa.variables['latitude'][:]
    lon = fa.variables['longitude'][:]
    aprecip=fa.variables[fieldreq][:]
    bprecip=fb.variables[fieldreq][:]
    cprecip=fc.variables[fieldreq][:]
    aprecip=np.squeeze(aprecip)
    bprecip=np.squeeze(bprecip)
    cprecip=np.squeeze(cprecip)
    ntimes,ny,nx=np.shape(aprecip)
    print(ntimes,ny,nx)

    fa.close()
    fb.close()
    fc.close()

    
#average across the time dimension
    pi_aprecip_avg=np.mean(aprecip,axis=0)
    pi_bprecip_avg=np.mean(bprecip,axis=0)
    pi_cprecip_avg=np.mean(cprecip,axis=0)
    

    pi_seasprecip=np.mean((pi_aprecip_avg,pi_bprecip_avg,pi_cprecip_avg),axis=0)
    pi_seasprecip=pi_seasprecip * 60. * 60. * 30. * 24.
    
    
    lonprecip=lon
    pi_seasprecip,lon = shiftgrid(180.,pi_seasprecip,lon,start=False)
    
    
    

     #==============
     # Pliocene+2

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/precip_data/'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+pliop2_expt+'_netcdf/'

    filestart=filestart+pliop2_expt+'a@pd'+extra+'*'

    print(filestart)
    fa=MFDataset(filestart+m1+fileend)
    fb=MFDataset(filestart+m2+fileend)
    fc=MFDataset(filestart+m3+fileend)
    aprecip=fa.variables[fieldreq][:]
    bprecip=fb.variables[fieldreq][:]
    cprecip=fc.variables[fieldreq][:]
    aprecip=np.squeeze(aprecip)
    bprecip=np.squeeze(bprecip)
    cprecip=np.squeeze(cprecip)

    fa.close()
    fb.close()
    fc.close()

    
    pliop2_aprecip_avg=np.mean(aprecip,axis=0)
    pliop2_bprecip_avg=np.mean(bprecip,axis=0)
    pliop2_cprecip_avg=np.mean(cprecip,axis=0)


    pliop2_seasprecip=np.mean((pliop2_aprecip_avg,pliop2_bprecip_avg,pliop2_cprecip_avg),axis=0)
    pliop2_seasprecip=pliop2_seasprecip * 60. * 60. * 30. * 24.

    lon=lonprecip
    pliop2_seasprecip,lon = shiftgrid(180.,pliop2_seasprecip,lon,start=False)
    

    
    # Pliocene+2 - preindustrial

    pliop2_anom=pliop2_seasprecip-pi_seasprecip
    #V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]
    V=[-128,-64,-32,-16,-8,-4,0,4,8,16,32,64,128]
   

  


    ###############################################
    # plot anomaly in mm/day for paper.
    # also calculate mean changes

    V=0
    titlename=seasname
    plotmap_itcz(pliop2_anom/30.,figureno,lon,lat,titlename,-2.0,2.2,0.2,V,'a','mm/day',mask_ind,10,12,figureno*2,0,9,2)
   

    # do line graph showing absolute value by latitude

   
    ax=plt.subplot2grid((10,12),(figureno*2,10),colspan=2,rowspan=2)
    ax.plot(np.mean(pliop2_seasprecip/30.,axis=1),lat,label='mPWP',linewidth=0.8)
    ax.plot(np.mean(pi_seasprecip/30.,axis=1),lat,label='PI',linewidth=0.8)
    # plot dotted lines at locations
    xmax=10
    xmin=0
    ax.plot([xmin,xmax],[-15,-15],color='black',linestyle='dotted',linewidth=0.8) 
    ax.plot([xmin,xmax],[0,0],color='black',linestyle='dotted',linewidth=0.8) 
    ax.plot([xmin,xmax],[15,15],color='black',linestyle='dotted',linewidth=0.8) 
   
    ax.tick_params(axis='y',which='both',labelleft='off')
    ax.set_ylim(-30,30)
    ax.set_xlim(xmin,xmax)
    if figureno !=3:
        ax.set_xticks([]) # disable xticks
    else:
        ax.set_xticks(list(range(xmin,xmax,2)))
        #ax.set_xlabel("mm/day",va='top')

    # add legend away from plot
    ax_leg=plt.subplot2grid((10,12),(9,10),colspan=2,rowspan=1)
    ax_leg.legend(*ax.get_legend_handles_labels(),loc='center')
    ax_leg.axis('off')
    plt.legend()
    


     


#end def focus itcz


def nonanomitcz(m1,figureno,monname,preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind):

  
    # m1 is the month 
    #==============
    # preindustrial

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/precip_data/'
    fileend='_precip.nc'
    fieldreq='precip_1'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+preind_expt+'_netcdf/'
        fileend='.nc'
        fieldreq='precip'
    
    filestart=filestart+preind_expt+'a@pd'+extra+'*'

   
    # read in data from multiple files
    fa=MFDataset(filestart+m1+fileend)
    lat = fa.variables['latitude'][:]
    lon = fa.variables['longitude'][:]
    aprecip=fa.variables[fieldreq][:]
    aprecip=np.squeeze(aprecip)
    ntimes,ny,nx=np.shape(aprecip)

    fa.close()

    
#average across the time dimension
    pi_monprecip=np.mean(aprecip,axis=0)
    pi_monprecip=pi_monprecip * 60. * 60. * 24. # mm/day

    lonprecip=lon
    pi_monprecip,lon = shiftgrid(180.,pi_monprecip,lon,start=False)
    
    
     #==============
     # Pliocene+2

    filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/precip_data/'
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/netcdf/'+pliop2_expt+'_netcdf/'

    filestart=filestart+pliop2_expt+'a@pd'+extra+'*'


    fa=MFDataset(filestart+m1+fileend)
    aprecip=fa.variables[fieldreq][:]
    aprecip=np.squeeze(aprecip)
    fa.close()
    
    pliop2_monprecip=np.mean(aprecip,axis=0)
    pliop2_monprecip=pliop2_monprecip * 60. * 60. * 24.

    lon=lonprecip
    pliop2_monprecip,lon = shiftgrid(180.,pliop2_monprecip,lon,start=False)
    


    ###############################################
    # plot precip in mm/day for paper.
    # also calculate mean changes

    V=0
    titlename=monname+' mPWP'
    plotmap_itcz(pliop2_monprecip,figureno,lon,lat,titlename,0,15,0.1,V,'n','mm/day',mask_ind,4,1,0,0,1,1)
   
    titlename=monname+' PI'
    plotmap_itcz(pi_monprecip,figureno,lon,lat,titlename,0,15,0.1,V,'n','mm/day',mask_ind,4,1,1,0,1,1)

    titlename=monname+' anomaly'
    plotmap_itcz(pliop2_monprecip-pi_monprecip,figureno,lon,lat,titlename,-2.0,2,0.1,V,'a','mm/day',mask_ind,4,1,2,0,1,1)
   
  
#end def nonanomitcz

################################
# main program

# annual mean
#figureno=0
preind_expt='xkvje'
plio_expt='xkvjf'
pliop2_expt='xkvjg'
extra='n'
HadCM3='n'

#preind_expt='xiboi'
#plio_expt='xibol'
#pliop2_expt='xibol'
#extra='y'
#HadCM3='y'

#plt.figure(figureno)
#annmean(figureno,preind_expt,plio_expt,pliop2_expt,extra,'b')
#figureno=figureno+1

#djf mean
#mask_ind='b'  # values are l=land, o=ocean b=both
#seasmean('dc','ja','fb',figureno,'djf',preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind)

#mam mean
#seasmean('mr','ar','my',figureno,'mam',preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind)

#jja mean
#seasmean('jn','jl','ag',figureno,'jja',preind_expt,plio_expt,pliop2_expt,extra,#HadCM3,mask_ind)


#son mean
#seasmean('sp','ot','nv',figureno,'son',preind_expt,plio_expt,pliop2_expt,extra,#HadCM3,mask_ind)

######################################################################
# try and focus on the tropical (ITCZ) region

#fig=plt.figure()   
#mask_ind='t'  # values are l=land, o=ocean b=both t=tropics
#figureno=0
#focusitcz('dc','ja','fb',figureno,'DJF',preind_expt,plio_expt,pliop2_expt,extra#,HadCM3,mask_ind)
#focusitcz('mr','ar','my',1,'MAM',preind_expt,plio_expt,pliop2_expt,extra,HadCM3#,mask_ind)
#focusitcz('jn','jl','ag',2,'JJA',preind_expt,plio_expt,pliop2_expt,extra,HadCM3#,mask_ind)
#focusitcz('sp','ot','nv',3,'SON',preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind)


#fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/ITCZ#_diag_'+pliop2_expt+'_allseasons.eps' 
#plt.savefig(fileout, bbox_inches='tight')  

#fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfprecip/ITCZ#_diag_'+pliop2_expt+'_allseasons.tiff' 
#plt.savefig(fileout, bbox_inches='tight')  
  

#plt.close()



######################################################################
# have a look at non anomaly plots for the tropical region.
# this will just check the shape of the ITCZ

#fig=plt.figure()   
#mask_ind='talt'  # values are l=land, o=ocean b=both t=tropics
#HadCM3='n'
#figureno=0
#nonanomitcz('ja',figureno,'january',preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind)   


#monthnames=['fb','mr','ar','my','jn','jl','ag','sp','t','nv','dc']

#for month in monthnames:
#    print(month)
#    fig=plt.figure()
#    nonanomitcz(month,0,month,preind_expt,plio_expt,pliop2_expt,extra,HadCM3,mask_ind)   


##############################################################
# compare snowfall vs rain (particularly in winter

#djf mean
mask_ind='l'  # values are l=land, o=ocean b=both
seas_snowrain('dc','ja','fb','djf',preind_expt,pliop2_expt,extra,HadCM3,mask_ind)

#mam mean
seas_snowrain('mr','ar','my','mam',preind_expt,pliop2_expt,extra,HadCM3,mask_ind)

#jja mean
seas_snowrain('jn','jl','ag','jja',preind_expt,pliop2_expt,extra,HadCM3,mask_ind)


#son mean
seas_snowrain('sp','ot','nv','son',preind_expt,pliop2_expt,extra,HadCM3,mask_ind)



#plt.show()
#sys.exit(0)

####

