#!/usr/bin/env python2.7
#NAME
#    PLOT_SURFPMINE
#PURPOSE
#    This program will plot the p-e (annual and seasonal) and
#    the p-e anomaly (annual and seasonal) for the final 30 years
#    of the HadGEM2 simulations
#
#    also by latitude
#
#
# search for 'main program' to find end of functions
# Julia 4/7/2017  (note annmean and seasmean not done, 
#  just done change with latitude



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
    if fileno != 99:
        plt.subplot(2,2,fileno+1)

   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='h')
   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='l')
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary
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
                print(np.shape(plotdata))
                cs = map.contourf(x,y,plotdata,V)
                cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
#end def plotdata


###########################################################
def plotmap_nh_precipevap(plotprecip,plotevap,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,mask_ind,ygrid,xgrid,yuse,xuse,yspan,xspan):

    # setup

    if fileno > 3:
        print('ERROR NOT ENOUGH SPACE ON PAGE ',fileno)
        sys.exit()
    lons, lats = np.meshgrid(lon,lat)
   
    if mask_ind =='n': # tropics mask
        northlat=90.0
        southlat=30.0
    elif mask_ind == 'nalt': # alternative tropical region:
        northlat=90.0
        southlat=0.0
    else:
        northlat=90.0
        southlat=-90.0

    if V == 0:
        V=np.arange(minval,maxval,valinc)


    # plot precip
    plt.subplot2grid((ygrid,xgrid),(yuse,xuse),colspan=yspan,rowspan=xspan)

    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=southlat,urcrnrlat=northlat,projection='cyl',resolution='c',fix_aspect=False)
    x, y = map(lons, lats)
    map.drawcoastlines()
        
    if uselog == 'n':
        cs = map.contourf(x,y,plotprecip,V,cmap='YlGnBu',extend='both')
    else:
        cs = map.contourf(x,y,plotprecip,V,cmap='RdBu',extend='both')
   
    parallels=(30,45,60,75)
    map.drawparallels(parallels,labels=[False,False,False,False]) # labels right
    meridians=np.arange(-180.,180.,90.)
    if fileno !=3:
        map.drawmeridians(meridians,labels=[False,False,False,False]) # nolabels
    else:
        map.drawmeridians(meridians,labels=[False,False,False,True]) # labels bottom
   
    fontsize=10
    plt.text(-180.0-6,northlat-fontsize-1,titlename,fontsize=fontsize,ha="right",bbox=dict(boxstyle="square,pad=0.1",color="white"))
 
    plt.text(-180.0-6,northlat-fontsize-1-15,'(precip)',fontsize=fontsize,ha="right",bbox=dict(boxstyle="square,pad=0.1",color="white"))

    # plot map boundary
    map.drawmapboundary


    # plot evap
    print(ygrid,xgrid,yuse,xuse,yspan,xspan)
    plt.subplot2grid((ygrid,xgrid),(yuse,xuse+yspan+1),colspan=yspan,rowspan=xspan)

    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=southlat,urcrnrlat=northlat,projection='cyl',resolution='c',fix_aspect=False)
    x, y = map(lons, lats)
    map.drawcoastlines()
        
    if uselog == 'n':
        cs = map.contourf(x,y,plotevap,V,cmap='YlGnBu',extend='both')
    else:
        cs = map.contourf(x,y,plotevap,V,cmap='RdBu',extend='both')
   
    parallels=(30,45,60,75)
    map.drawparallels(parallels,labels=[False,False,False,False]) # no labels
    meridians=np.arange(-180.,180.,90.)
    if fileno !=3:
        map.drawmeridians(meridians,labels=[False,False,False,False]) # nolabels
    else:
        map.drawmeridians(meridians,labels=[False,False,False,True]) # labels bottom
   
    fontsize=10
    plt.text(-180.0-6,northlat-fontsize-1,titlename,fontsize=fontsize,ha="right",bbox=dict(boxstyle="square,pad=0.1",color="white"))

    plt.text(-180.0-6,northlat-fontsize-1-15,'(evap)',fontsize=fontsize,ha="right",bbox=dict(boxstyle="square,pad=0.1",color="white"))

  
  # plot map boundary
    map.drawmapboundary


    # colorbar
    if fileno==0:
        plt.subplot2grid((10,11),(9,0),colspan=11,rowspan=1)
        plt.gca().set_visible(False)
        cbar = plt.colorbar(cs,orientation="horizontal",fraction=1.0)         
        cbar.set_label(cbarname)
        cbar.ax.tick_params(labelsize=10)
   
    
  

    
    
   


#end def plotmap_nh_precipevap


def annmean(figureno,HadCM3):
    #==============
    # preindustrial


    # read in data from multiple files
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
        aprecip=f.variables['precip'][:]
        titlepi='PI-Ann_HadCM3'
        titleplio='Plio-Ann_HadCM3'
        titlediff='Plio - preind  Ann_HadCM3'
    else:
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/xkvjea@pdn[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
        aprecip=f.variables['precip_1'][:]
        titlepi='PI-Ann_HadGEM2'
        titleplio='Plio-Ann_HadGEM2'
        titlediff='Plio - preind Ann_HadGEM2'


    aevap=np.squeeze(aevap)    
    aprecip=np.squeeze(aprecip)
    ntimes,ny,nx=np.shape(aevap)
    print(ntimes,ny,nx)
    f.close()
    print('read in preindustrial')
  
    apme=aprecip-aevap
#average across the time dimension
    pi_pme_ann=np.mean(apme,axis=0)
    print('new shape',np.shape(pi_pme_ann))
    
    pi_pme_ann=pi_pme_ann * 60. * 60. * 24. * 30.
    
    plt.figure(0)
    lonevap=lon
    pi_pme_ann,lon = shiftgrid(180.,pi_pme_ann,lon,start=False)
    
    plotdata(pi_pme_ann,0,lon,lat,titlepi,0,275,25.0,0.0,'n','mm/month')
    print('first plot done')


     #==============
     # Pliocene

    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aevap=np.squeeze(aevap)
        aprecip=f.variables['precip'][:]
        aprecip=np.squeeze(aprecip)
    else:
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjf/netcdf/xkvjfa@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aevap=np.squeeze(aevap)
        aprecip=f.variables['precip_1'][:]
        aprecip=np.squeeze(aprecip)
    ntimes,ny,nx=np.shape(aevap)
    print(ntimes,ny,nx)

    apme=aprecip-aevap
    plio_pme_ann=np.mean(apme,axis=0)
    plio_pme_ann=plio_pme_ann * 60. * 60. * 24. * 30.

    lon=lonevap
    plio_pme_ann,lon = shiftgrid(180.,plio_pme_ann,lon,start=False)

    plotdata(plio_pme_ann,1,lon,lat,titleplio,0,275,25,0.0,'n','mm/month')
    f.close()
    print('read in plio')


     #==============
     # Pliocene+2

    if HadCM3 != 'y':
        # read in data from multiple files
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/xkvjga@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aevap=np.squeeze(aevap)
        aprecip=f.variables['precip_1'][:]
        aprecip=np.squeeze(aprecip)
        ntimes,ny,nx=np.shape(aevap)
        print(ntimes,ny,nx)
        print(np.shape(aprecip))
        print(np.shape(aevap))
        
    #average across the time dimension
        apme=aprecip-aevap
        plio_pmep2_ann=np.mean(apme,axis=0)
        plio_pmep2_ann=plio_pmep2_ann * 60. * 60. * 24. * 30.
        lon=lonevap
        plio_pmep2_ann,lon = shiftgrid(180.,plio_pmep2_ann,lon,start=False)


    f.close()
    print('read in plio+2')


    # Pliocene - preindustrial

    plio_anom=plio_pme_ann-pi_pme_ann

    V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]
    plotdata(plio_anom,2,lon,lat,titlediff,0,275,25,V,'n','mm/month')
    
    # Pliocene+2 - preindustrial

    if HadCM3 != 'y':
        pliop2_anom=plio_pmep2_ann-pi_pme_ann
        V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]

        plotdata(pliop2_anom,3,lon,lat,'PlioP2 - PI Panom_HG2',0,275,25,V,'n','mm/month')


    if HadCM3 == 'y':
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_HadCM3.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_HadCM3.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  
    else:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  

    plt.close()



    # Pliocene - preindustrial percentage change for paper

    #print(pi_pme_ann)
    #plio_peranom=((plio_pme_ann-pi_pme_ann)/pi_pme_ann)*100.
    #plotdata(plio_peranom,0,lon,lat,titlediff,-50,60,5,0,'a','%')

    # quickly check that globally averaged p-e is 0
    # create weighting array
    weightarr=np.zeros(np.shape(plio_pmep2_ann))
    for i in range(0,len(lon)):
        weightarr[:,i]=np.cos(np.deg2rad(lat))

    print('mean plio pme',np.average(plio_pmep2_ann,weights=weightarr))
    print('mean pi pme',np.average(pi_pme_ann,weights=weightarr))
    print('mean plio-pi pme',np.average(plio_pmep2_ann-pi_pme_ann,weights=weightarr))


    if HadCM3 != 'y':
        # without mask
        plio_peranomp2=((plio_pmep2_ann-pi_pme_ann)/pi_pme_ann)*100.
        plotdata(plio_peranomp2,99,lon,lat,'Plio - PI Panom_HG2+2 %',-75,78,5,0,'a','%')

        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  
        plt.close()

        # with mask for if pi -pme is close to zero
        plio_peranomp2_mask=plio_peranomp2

        for i in range(0,len(lon)):
            for j in range(0,len(lat)):
                if pi_pme_ann[j,i] < 2 and plio_pmep2_ann[j,i]<2:
                    plio_peranomp2_mask[j,i]=0.

        
        plotdata(plio_peranomp2_mask,99,lon,lat,'Plio - PI Panom_HG2+2 %',-75,78,5,0,'a','%')
     
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent_mask.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent_mask.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  



        plt.close()

        # in mm/month

        plio_anomp2=(plio_pmep2_ann-pi_pme_ann)
        
        plotdata(plio_anomp2,99,lon,lat,'Plio - PI Panom_HG2+2',-20,25,5,0,'a','mm/month')
     
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom2.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom2.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  



        plt.close()



#end def annmean


def seasmean(m1,m2,m3,figureno,seasname,HadCM3):
    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean
    #==============
    # preindustrial

   
    # read in data from multiple files
    if HadCM3 == 'y':
        fa=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*'+m1+'.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*'+m2+'.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*'+m3+'.nc')
        lat = fa.variables['latitude'][:]
        lon = fa.variables['longitude'][:]
        aevap=fa.variables['field184'][:]
        bevap=fb.variables['field184'][:]
        cevap=fc.variables['field184'][:]
        pititle='PI HadCM3: '+seasname
        pliotitle='Plio HadCM3: '+seasname
        difftitle='Plio-PI HadCM3: '+seasname
        
    else:
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvja/netcdf/xkvjaa@pdm[7-9]*'+m1+'_evap.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvja/netcdf/xkvjaa@pdm[7-9]*'+m2+'_evap.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvja/netcdf/xkvjaa@pdm[7-9]*'+m3+'_evap.nc')
        lat = fa.variables['latitude'][:]
        lon = fa.variables['longitude'][:]
        aevap=fa.variables['field184'][:]
        bevap=fb.variables['field184'][:]
        cevap=fc.variables['field184'][:]
        pititle='PI HadGEM2: '+seasname
        pliotitle='Plio HadGEM2: '+seasname
        difftitle='Plio-PI HadGEM2: '+seasname




    aevap=np.squeeze(aevap)
    bevap=np.squeeze(bevap)
    cevap=np.squeeze(cevap)
    ntimes,ny,nx=np.shape(aevap)
    print(ntimes,ny,nx)
    
#average across the time dimension
    pi_aevap_avg=np.mean(aevap,axis=0)
    pi_bevap_avg=np.mean(bevap,axis=0)
    pi_cevap_avg=np.mean(cevap,axis=0)
    
    pi_seasevap=np.mean((pi_aevap_avg,pi_bevap_avg,pi_cevap_avg),axis=0)
    pi_seasevap=pi_seasevap * 60. * 60. * 30. * 24.
    
    
    lonevap=lon
    pi_seasevap,lon = shiftgrid(180.,pi_seasevap,lon,start=False)
    
    plotdata(pi_seasevap,0,lon,lat,pititle,0,275,25,0.0,'n','mm/month')
    
     #==============
     # Pliocene

    if HadCM3 == 'y':
        fa=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[7-9]*'+m1+'.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[7-9]*'+m2+'.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[7-9]*'+m3+'.nc')
        aevap=fa.variables['field184'][:]
        bevap=fb.variables['field184'][:]
        cevap=fc.variables['field184'][:]
    else:
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjb/netcdf/xkvjba@pdm[7-9]*'+m1+'_evap.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjb/netcdf/xkvjba@pdm[7-9]*'+m2+'_evap.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjb/netcdf/xkvjba@pdm[7-9]*'+m3+'_evap.nc')
        aevap=fa.variables['field184'][:]
        bevap=fb.variables['field184'][:]
        cevap=fc.variables['field184'][:]

    aevap=np.squeeze(aevap)
    bevap=np.squeeze(bevap)
    cevap=np.squeeze(cevap)

    # average across the time dimension    
    plio_aevap_avg=np.mean(aevap,axis=0)
    plio_bevap_avg=np.mean(bevap,axis=0)
    plio_cevap_avg=np.mean(cevap,axis=0)
    
    plio_seasevap=np.mean((plio_aevap_avg,plio_bevap_avg,plio_cevap_avg),axis=0)
    plio_seasevap=plio_seasevap * 60. * 60. * 30. * 24.

    lon=lonevap
    plio_seasevap,lon = shiftgrid(180.,plio_seasevap,lon,start=False)
    
    
    plotdata(plio_seasevap,1,lon,lat,pliotitle,0,275,25,0.0,'n','mm/month')



     #==============
     # Pliocene+2

    if HadCM3 != 'y':
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjc/netcdf/xkvjca@pdm[7-9]*'+m1+'_evap.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjc/netcdf/xkvjca@pdm[7-9]*'+m2+'_evap.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjc/netcdf/xkvjca@pdm[7-9]*'+m3+'_evap.nc')
        aevap=fa.variables['field184'][:]
        bevap=fb.variables['field184'][:]
        cevap=fc.variables['field184'][:]
        aevap=np.squeeze(aevap)
        bevap=np.squeeze(bevap)
        cevap=np.squeeze(cevap)
        
        pliop2_aevap_avg=np.mean(aevap,axis=0)
        pliop2_bevap_avg=np.mean(bevap,axis=0)
        pliop2_cevap_avg=np.mean(cevap,axis=0)
        
        pliop2_seasevap=np.mean((pliop2_aevap_avg,pliop2_bevap_avg,pliop2_cevap_avg),axis=0)
        pliop2_seasevap=pliop2_seasevap * 60. * 60. * 30. * 24.
        
        lon=lonevap
        pliop2_seasevap,lon = shiftgrid(180.,pliop2_seasevap,lon,start=False)
    
 

    # Pliocene - preindustrial

    plio_anom=plio_seasevap-pi_seasevap
    V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]
    plotdata(plio_anom,2,lon,lat,difftitle,0,275,25,V,'la','mm/month')
    
    # Pliocene+2 - preindustrial

    if HadCM3 != 'y':
        pliop2_anom=pliop2_seasevap-pi_seasevap
        V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]

        plotdata(pliop2_anom,3,lon,lat,'PlioP2 - PI Panom_HG2',0,275,25,V,'la','mm/month')


    if HadCM3 == 'y':
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_'+seasname+'anom_HadCM3.eps' 
    else:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_'+seasname+'anom.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

    # Pliocene - preindustrial percentage change

    plio_peranom=((plio_seasevap-pi_seasevap)/pi_seasevap)*100.
    plotdata(plio_peranom,0,lon,lat,difftitle,-50,50,5,0,'a','%')


    if HadCM3 !='y':
        plio_peranomp2=((pliop2_seasevap-pi_seasevap)/pi_seasevap)*100.
        titlename='Plio+2-PI evap %'+seasname
        plotdata(plio_peranomp2,1,lon,lat,titlename,-50,50,5,0,'a','%')

    if HadCM3 =='y':
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent_'+seasname+'HadCM3.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent_'+seasname+'HadCM3.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  
    else:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_anom_percent_'+seasname+'.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  


    plt.close()




#end def seasmean


def seasmean_get(m1,m2,m3,figureno,seasname,HadCM3,moses2,field):
    # this program is a bit like seasmean but it will just get the data
    # it won't analyse it or print anything out or plot anything

    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean
    #==============
    # preindustrial

   
    # read in data from multiple files
    if HadCM3 == 'y':
        if moses2 == 'y':
            fa=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*'+m1+'.nc')
            fb=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*'+m2+'.nc')
            fc=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*'+m3+'.nc')
       
        else: # fergus simulation xgrad
            fa=MFDataset('/nfs/hera1/earjcti/um/netcdf/xgrad_netcdf/pdfiles/xgrada@pdq[7-9]*'+m1+'.nc')
            fb=MFDataset('/nfs/hera1/earjcti/um/netcdf/xgrad_netcdf/pdfiles/xgrada@pdq[7-9]*'+m2+'.nc')
            fc=MFDataset('/nfs/hera1/earjcti/um/netcdf/xgrad_netcdf/pdfiles/xgrada@pdq[7-9]*'+m3+'.nc')
      
        lat = fa.variables['latitude'][:]
        lon = fa.variables['longitude'][:]
        
    else:
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/pdfiles/xkvjea@pdn[7-9]*'+m1+'.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/pdfiles/xkvjea@pdn[7-9]*'+m2+'.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/pdfiles/xkvjea@pdn[7-9]*'+m3+'.nc')
        lat = fa.variables['latitude'][:]
        lon = fa.variables['longitude'][:]
        if field=='precip':
            field='precip_1'
    aprecip=fa.variables[field][:]
    bprecip=fb.variables[field][:]
    cprecip=fc.variables[field][:]




    aprecip=np.squeeze(aprecip)
    bprecip=np.squeeze(bprecip)
    cprecip=np.squeeze(cprecip)
    ntimes,ny,nx=np.shape(aprecip)
    
#average across the time dimension
    pi_aprecip_avg=np.mean(aprecip,axis=0)
    pi_bprecip_avg=np.mean(bprecip,axis=0)
    pi_cprecip_avg=np.mean(cprecip,axis=0)
    
    pi_seasprecip=np.mean((pi_aprecip_avg,pi_bprecip_avg,pi_cprecip_avg),axis=0)
    pi_seasprecip=pi_seasprecip * 60. * 60. * 24. # mm/day
    pi_seasprecip,lon = shiftgrid(180.,pi_seasprecip,lon,start=False)
   
    
     #==============
     # Pliocene

    if HadCM3 == 'y':
        if moses2 =='y':
            fa=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[7-9]*'+m1+'.nc')
            fb=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[7-9]*'+m2+'.nc')
            fc=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[7-9]*'+m3+'.nc')
        else:
            fa=MFDataset('/nfs/hera1/earjcti/um/netcdf/xgrac_netcdf/pdfiles/xgraca@pdt[7-9]*'+m1+'.nc')
            fb=MFDataset('/nfs/hera1/earjcti/um/netcdf/xgrac_netcdf/pdfiles/xgraca@pdt[7-9]*'+m2+'.nc')
            fc=MFDataset('/nfs/hera1/earjcti/um/netcdf/xgrac_netcdf/pdfiles/xgraca@pdt[7-9]*'+m3+'.nc')
       
    else:
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/pdfiles/xkvjga@pdn[7-9]*'+m1+'.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/pdfiles/xkvjga@pdn[7-9]*'+m2+'.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/pdfiles/xkvjga@pdn[7-9]*'+m3+'.nc')
    aprecip=fa.variables[field][:]
    bprecip=fb.variables[field][:]
    cprecip=fc.variables[field][:]
    lon = fa.variables['longitude'][:]
       

    aprecip=np.squeeze(aprecip)
    bprecip=np.squeeze(bprecip)
    cprecip=np.squeeze(cprecip)

    # average across the time dimension    
    plio_aprecip_avg=np.mean(aprecip,axis=0)
    plio_bprecip_avg=np.mean(bprecip,axis=0)
    plio_cprecip_avg=np.mean(cprecip,axis=0)
    
    plio_seasprecip=np.mean((plio_aprecip_avg,plio_bprecip_avg,plio_cprecip_avg),axis=0)
    plio_seasprecip=plio_seasprecip * 60. * 60.  * 24.
   
    plio_seasprecip,lon = shiftgrid(180.,plio_seasprecip,lon,start=False)
   
    retdata=[lon,lat,pi_seasprecip,plio_seasprecip]
    return retdata
    
#end def seasmean_get

#####################################
def nh_allseasons_precipevap(HadMC3,moses2):
# plot precip and evap amount for all seasons
    if HadCM3=='n':
        plio_expt='xkvjg'
        precip_field='precip_1'
    else:
        plio_expt='HadCM3'+moses2
        precip_field='precip'

    # getdata for all seasons in nh we are going to put it on a nice figure
    retdata=seasmean_get('dc','ja','fb',figureno,'djf',HadCM3,moses2,precip_field)
    lon=retdata[0]
    lat=retdata[1]
    pi_djf_precip=retdata[2]
    plio_djf_precip=retdata[3]
    retdata=seasmean_get('dc','ja','fb',figureno,'djf',HadCM3,moses2,'field184')
    pi_djf_evap=retdata[2]
    plio_djf_evap=retdata[3]
    

    mask_ind='n' # northernhemisphere
    plotmap_nh_precipevap(plio_djf_precip-pi_djf_precip,plio_djf_evap-pi_djf_evap,0,lon,lat,'DJF',-1.0,1.1,0.1,0,'a','mm/day',mask_ind,10,11,0,0,5,2)
   
    retdata=seasmean_get('mr','ar','my',figureno,'mam',HadCM3,moses2,precip_field)
    pi_mam_precip=retdata[2]
    plio_mam_precip=retdata[3]
    retdata=seasmean_get('mr','ar','my',figureno,'mam',HadCM3,moses2,'field184')
    pi_mam_evap=retdata[2]
    plio_mam_evap=retdata[3]
    plotmap_nh_precipevap(plio_mam_precip-pi_mam_precip,plio_mam_evap-pi_mam_evap,1,lon,lat,'MAM',-1.0,1.1,0.1,0,'a','mm/day',mask_ind,10,11,2,0,5,2)
      
 
    retdata=seasmean_get('jn','jl','ag',figureno,'jja',HadCM3,moses2,precip_field)
    pi_jja_precip=retdata[2]
    plio_jja_precip=retdata[3] 
    retdata=seasmean_get('jn','jl','ag',figureno,'jja',HadCM3,moses2,'field184')
    pi_jja_evap=retdata[2]
    plio_jja_evap=retdata[3]
    plotmap_nh_precipevap(plio_jja_precip-pi_jja_precip,plio_jja_evap-pi_jja_evap,2,lon,lat,'JJA',-1.0,1.1,0.1,0,'a','mm/day',mask_ind,10,11,4,0,5,2)
       

    retdata=seasmean_get('sp','ot','nv',figureno,'son',HadCM3,moses2,precip_field)
    pi_son_precip=retdata[2]
    plio_son_precip=retdata[3]
    retdata=seasmean_get('sp','ot','nv',figureno,'son',HadCM3,moses2,'field184')
    pi_son_evap=retdata[2]
    plio_son_evap=retdata[3]
    plotmap_nh_precipevap(plio_son_precip-pi_son_precip,plio_son_evap-pi_son_evap,3,lon,lat,'SON',-1.0,1.1,0.1,0,'a','mm/day',mask_ind,10,11,6,0,5,2)
   
  
   
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/NH_precipevap_'+plio_expt+'_allseasons.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/NH_precipevap'+plio_expt+'_allseasons.tiff' 
    plt.savefig(fileout, bbox_inches='tight')  



#enddef nh_allseasons_precipevap



def pmine_chg_by_latitude(HadCM3,land_ocn_ind,abs_pcent):


    if land_ocn_ind == 'l' or land_ocn_ind == 'o':  # land or ocean
        # get land mask
        if HadCM3 == 'y' :
            fm=Dataset('/nfs/hera2/apps/metadata/ancil/preind2/qrparm.mask.nc')
        else:
            fm=Dataset('/nfs/hera1/earjcti/um/HadGEM_ancils/qrparm.mask.nc')
        lsmlon=fm.variables['longitude'][:]
        lsmlat=fm.variables['latitude'][:]
        lsm=fm.variables['lsm'][:]
        lsm=np.squeeze(lsm)
        fm.close()

    #==============
    # preindustrial

    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
        aprecip=f.variables['precip'][:]
    else:
    # read in data from multiple files
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/xkvjea@pdn[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
        aprecip=f.variables['precip_1'][:]
    
    f.close()
    aevap=np.squeeze(aevap) 
    aprecip=np.squeeze(aprecip)
    pmine=aprecip-aevap

    ntimes,ny,nx=np.shape(aevap)
    print(ntimes,ny,nx)

#average across the time dimension and the latitude dimension
    pi_pmine_ann=np.mean(pmine,axis=0)

    if land_ocn_ind == 'l': # mask out all non land points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            pi_pmine_ann=pi_pmine_ann / (lsm)
            pi_pmine_ann[pi_pmine_ann == float('Inf')] = float('NaN')
            pi_pmine_ann[pi_pmine_ann == float('-Inf')] = float('NaN')
        else:
            print('error lon/lat of land sea mask dont match')

    if land_ocn_ind == 'o': # mask out all non ocean points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            pi_pmine_ann=pi_pmine_ann / (np.abs(lsm-1.0))
            pi_pmine_ann[pi_pmine_ann == float('Inf')] = float('NaN')
            pi_pmine_ann[pi_pmine_ann == float('-Inf')] = float('NaN')
        else:
            print('error lon/lat of land sea mask dont match')

  
    

    pi_pmine_lat_ann=np.nanmean(pi_pmine_ann,axis=1)
    print('new shape pi',np.shape(pi_pmine_lat_ann))
    
     #==============
     # Pliocene

    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
        aprecip=f.variables['precip'][:]
    else:
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjf/netcdf/xkvjfa@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aprecip=f.variables['precip_1'][:]

    f.close()
    aevap=np.squeeze(aevap)
    aprecip=np.squeeze(aprecip)
    pmine=aprecip-aevap

    ntimes,ny,nx=np.shape(aevap)

    plio_pmine_ann=np.mean(pmine,axis=0)


    if land_ocn_ind == 'l': # mask out all non land points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            plio_pmine_ann=plio_pmine_ann / (lsm)
            plio_pmine_ann[plio_pmine_ann == float('Inf')] = float('NaN')
            plio_pmine_ann[plio_pmine_ann == float('-Inf')] = float('NaN')
        else:
            print('error lon/lat of land sea mask dont match')

    if land_ocn_ind == 'o': # mask out all non ocean points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            plio_pmine_ann=plio_pmine_ann / (np.abs(lsm-1.0))
            plio_pmine_ann[plio_pmine_ann == float('Inf')] = float('NaN')
            plio_pmine_ann[plio_pmine_ann == float('-Inf')] = float('NaN')
        else:
            print('error lon/lat of land sea mask dont match')


    plio_pmine_lat_ann=np.nanmean(plio_pmine_ann,axis=1)


     #==============
     # Pliocene+2

    if HadCM3 != 'y':
        # read in data from multiple files
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/xkvjga@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aprecip=f.variables['precip_1'][:]
        aevap=np.squeeze(aevap)
        aprecip=np.squeeze(aprecip)
        pmine=aprecip-aevap

        ntimes,ny,nx=np.shape(aevap)
        print(ntimes,ny,nx)
        
    #average across the time dimension
        plio_evapp2_ann=np.mean(pmine,axis=0)
        if land_ocn_ind == 'l': # mask out all non land points.
            if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
                plio_pmep2_ann=plio_pmep2_ann / (lsm)
                plio_pmep2_ann[plio_pmep2_ann == float('Inf')] = float('NaN')
                plio_pmep2_ann[plio_pmep2_ann == float('-Inf')] = float('NaN')
            
            else:
                print('error lon/lat of land sea mask dont match')

        if land_ocn_ind == 'o': # mask out all non ocean points.
            if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
                plio_pmep2_ann=plio_pmep2_ann / (np.abs(lsm-1.0))
                plio_pmep2_ann[plio_pmep2_ann == float('Inf')] = float('NaN')
                plio_pmep2_ann[plio_pmep2_ann == float('-Inf')] = float('NaN')

            else:
                print('error lon/lat of land sea mask dont match')


        plio_pmep2_lat_ann=np.nanmean(plio_pmep2_ann,axis=1)
        f.close()


    # Pliocene - preindustrial

    if abs_pcent == 'a':
        plio_anom=(plio_pmine_lat_ann-pi_pmine_lat_ann)* 60. *  60. *24.
    if abs_pcent == 'p':
        plio_anom=(plio_pmine_lat_ann-pi_pmine_lat_ann)/ pi_pmine_lat_ann
        plio_anom=plio_anom * 100.

    # Pliocene+2 - preindustrial

    if HadCM3 != 'y':
        if abs_pcent == 'a':
          pliop2_anom=(plio_pmep2_lat_ann-pi_pmine_lat_ann)* 60.*  60. *24.
        if abs_pcent == 'p':
          pliop2_anom=(plio_pmep2_lat_ann-pi_pmine_lat_ann)/ pi_pmine_lat_ann
          pliop2_anom=pliop2_anom * 100.



    if HadCM3 == 'y':
        allanom=[lat,plio_anom]
        return allanom
    else:
        allanom=[lat,plio_anom,pliop2_anom]
        return allanom
    
    



#end def evap_chg_by_latitude



################################
# main program

# annual mean
figureno=0

HadCM3='y'
moses2='y'
#plt.figure(figureno)
#annmean(figureno,HadCM3)
#figureno=figureno+1

#djf mean
#plt.figure(figureno)
#seasmean('dc','ja','fb',figureno,'djf',HadCM3)
#figureno=figureno+1

#jja mean
#plt.figure(figureno)
#seasmean('jn','jl','ag',figureno,'jja',HadCM3)
#figureno=figureno+1

############################################################
# just plot poleward of 30N for all seasons
# we would like precipitation and evaporation on the same plot
# ie all seasons on one page
#nh_allseasons(HadCM3,moses2) # get precipitation
nh_allseasons_precipevap(HadCM3,moses2) # get precipevap

##########################################################
# to see what the polar amplification is
#HadCM3='y'
#land_ocn_ind='l'    # valid values are l - land , o- ocean b=both
#abs_pcent='a'       # valid values are p - percentage chg, a=absolute change
                    # note percent change is not meaningful in p-e because 
                     # p-e can be negative for the preindustrial
#evap_ret=pmine_chg_by_latitude(HadCM3,land_ocn_ind,abs_pcent)
#lats_HadCM3=evap_ret[0]
#pmine_anom_HadCM3=evap_ret[1]
#plt.plot(pmine_anom_HadCM3,lats_HadCM3)

#HadCM3='n'
#evap_ret=pmine_chg_by_latitude(HadCM3,land_ocn_ind,abs_pcent)
#lats_HadGEM=evap_ret[0]
#pmine_anom_HadGEM=evap_ret[1]
#pmine_anom_HadGEMp2=evap_ret[2]


#mp.rcParams.update({'font.size':15})
#plt.plot(pmine_anom_HadGEMp2,lats_HadGEM,'g',label='HadGEM2')
#plt.plot(pmine_anom_HadCM3,lats_HadCM3,'r',label='HadCM3')
#plt.plot(pmine_anom_HadGEM,lats_HadGEM,'b')

#if abs_pcent == 'a':
#    plt.xlabel('mm/day',fontsize=20)
#if abs_pcent == 'p': 
#    plt.xlabel('%',fontsize=20)
#plt.ylabel('latitude',fontsize=20)

#if land_ocn_ind =='l':
#    plt.title('f) mPWP - preind land p-e anom',loc='left',fontsize=25)
#if land_ocn_ind =='o':
#    plt.title('g) mPWP - preind ocean p-e anom',loc='left',fontsize=25)
#if land_ocn_ind =='b':
#    plt.title('e) mPWP - preind p-e anom',loc='left',fontsize=25)
#axes=plt.gca()
#if abs_pcent == 'p':
#    axes.set_xlim(xmin=-60.0,xmax=100.0)
#axes.set_ylim([-80,80])
#legend=plt.legend()
#if land_ocn_ind =='l':
#    if abs_pcent =='a':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_polar_amp_land.eps' 
#if land_ocn_ind =='o':
#    if abs_pcent =='a':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmine/P-E_polar_amp_ocean.eps' 


#if land_ocn_ind =='b':
#    if abs_pcent == 'a':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfpmi#ne/P-E_polar_amp.eps' 
#plt.savefig(fileout, bbox_inches='tight')  

#if abs_pcent =='p':
#    plt.show()



sys.exit(0)

####

