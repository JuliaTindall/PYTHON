#!/usr/bin/env python2.7
#NAME
#    PLOT_SURFEVAP
#PURPOSE
#    This program will plot the evapitation (annual and seasonal) and
#    the evapitation anomaly (annual and seasonal) for the final 30 years
#    of the HadGEM2 simulations
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
   # this is good for the eurasia
   # map=Basemap(llcrnrlon=0.0,urcrnrlon=180.0,llcrnrlat=0.0,urcrnrlat=90.0,projection='cyl',resolution='l')
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary
    x, y = map(lons, lats)
    map.drawcoastlines()
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.),extend='both')
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
                print(np.shape(plotdata))
                cs = map.contourf(x,y,plotdata,V,cmap='spectral',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")

    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        cbar.set_label(cbarname,labelpad=-60,size=15)
        cbar.ax.tick_params(labelsize=15)
        plt.title(titlename,loc='left',fontsize=18)
  

#end def plotdata

def annmean(figureno,HadCM3,land_ocn_ind):

    print('j0')



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
    if land_ocn_ind == 'l':
        valmax=100.
        valmaxanom=12.
    else:
        valmax=250.
        valmaxanom=30.
        
    valdiff=valmax/10.
    valdiffanom=valmaxanom/6.

    #==============
    # preindustrial


    # read in data from multiple files
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
        titlepi='PI-Ann_HadCM3'
        titleplio='Plio-Ann_HadCM3'
        titlediff='Plio - preind  Ann_HadCM3'
    else:
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/pdfiles/xkvjea@pdn[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
        titlepi='PI-Ann_HadGEM2'
        titleplio='Plio-Ann_HadGEM2'
        titlediff='Plio - preind Ann_HadGEM2'


    aevap=np.squeeze(aevap)
    ntimes,ny,nx=np.shape(aevap)
    print(ntimes,ny,nx)
    
#average across the time dimension
    pi_evap_ann=np.mean(aevap,axis=0)

    if land_ocn_ind == 'l': # mask out all non land points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            pi_evap_ann=pi_evap_ann / (lsm)
            pi_evap_ann[pi_evap_ann == float('Inf')] = float('NaN')
            pi_evap_ann[pi_evap_ann == float('-Inf')] = float('NaN')
 
        else:
            print('error lon/lat of land sea mask dont match')

    if land_ocn_ind == 'o': # mask out all non ocean points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            pi_evap_ann=pi_evap_ann / (np.abs(lsm-1.0))
            pi_evap_ann[pi_evap_ann == float('Inf')] = float('NaN')
            pi_evap_ann[pi_evap_ann == float('-Inf')] = float('NaN')
 
        else:
            print('error lon/lat of land sea mask dont match')




    print('new shape',np.shape(pi_evap_ann))
    
    pi_evap_ann=pi_evap_ann * 60. * 60. * 24. * 30.
    
    plt.figure(0)
    lonevap=lon
    pi_evap_ann,lon = shiftgrid(180.,pi_evap_ann,lon,start=False)
    
    plotdata(pi_evap_ann,0,lon,lat,titlepi,0,valmax,valdiff,0.0,'n','mm/month')
    f.close()
    
     #==============
     # Pliocene

    print('j1')


    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aevap=np.squeeze(aevap)
    else:
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjf/netcdf/pdfiles/xkvjfa@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aevap=np.squeeze(aevap)
    ntimes,ny,nx=np.shape(aevap)
    print(ntimes,ny,nx)

    plio_evap_ann=np.mean(aevap,axis=0)
    plio_evap_ann=plio_evap_ann * 60. * 60. * 24. * 30.
    lon=lonevap

    if land_ocn_ind == 'l': # mask out all non land points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            plio_evap_ann=plio_evap_ann / (lsm)
            plio_evap_ann[plio_evap_ann == float('Inf')] = float('NaN')
            plio_evap_ann[plio_evap_ann == float('-Inf')] = float('NaN')
 
        else:
            print('error lon/lat of land sea mask dont match')

    if land_ocn_ind == 'o': # mask out all non ocean points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            plio_evap_ann=plio_evap_ann / (np.abs(lsm-1.0))
            plio_evap_ann[plio_evap_ann == float('Inf')] = float('NaN')
            plio_evap_ann[plio_evap_ann == float('-Inf')] = float('NaN')
 
        else:
            print('error lon/lat of land sea mask dont match')



    plio_evap_ann,lon = shiftgrid(180.,plio_evap_ann,lon,start=False)

    plotdata(plio_evap_ann,1,lon,lat,titleplio,0,valmax,valdiff,0.0,'n','mm/month')
    f.close

     #==============
     # Pliocene+2

    print('j2')

    if HadCM3 != 'y':
        # read in data from multiple files
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/pdfiles/xkvjga@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aevap=np.squeeze(aevap)
        ntimes,ny,nx=np.shape(aevap)
        print(ntimes,ny,nx)
        
    #average across the time dimension
        plio_evapp2_ann=np.mean(aevap,axis=0)
        plio_evapp2_ann=plio_evapp2_ann * 60. * 60. * 24. * 30.
        lon=lonevap
                

        if land_ocn_ind == 'l': # mask out all non land points.
            if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
                plio_evapp2_ann=plio_evapp2_ann / (lsm)
                plio_evapp2_ann[plio_evapp2_ann == float('Inf')] = float('NaN')
                plio_evapp2_ann[plio_evapp2_ann == float('-Inf')] = float('NaN')
                
            else:
                print('error lon/lat of land sea mask dont match')
                
        if land_ocn_ind == 'o': # mask out all non ocean points.
            if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
                plio_evapp2_ann=plio_evapp2_ann / (np.abs(lsm-1.0))
                plio_evapp2_ann[plio_evapp2_ann == float('Inf')] = float('NaN')
                plio_evapp2_ann[plio_evapp2_ann == float('-Inf')] = float('NaN')
 
            else:
                print('error lon/lat of land sea mask dont match')


        plio_evapp2_ann,lon = shiftgrid(180.,plio_evapp2_ann,lon,start=False)
        f.close()



    # Pliocene - preindustrial

    print('j3')

    plio_anom=plio_evap_ann-pi_evap_ann

    V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]
    V=0
    print(valmaxanom,valdiffanom)
    plotdata(plio_anom,2,lon,lat,titlediff,-1.0*valmaxanom,valmaxanom,valdiffanom,V,'n','mm/month')
    
    # Pliocene+2 - preindustrial

    print('j4')

    if HadCM3 != 'y':
        pliop2_anom=plio_evapp2_ann-pi_evap_ann
        V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]
        V=0

        plotdata(pliop2_anom,3,lon,lat,'PlioP2 - PI Eanom_HG2',-1.0*valmaxanom,valmaxanom,valdiffanom,V,'n','mm/month')


    if HadCM3 == 'y':
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_HadCM3.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_HadCM3.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  
    else:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

    if HadCM3 == 'y':
        plotdata(plio_anom,99,lon,lat,'b) HadCM3 evaporation anomaly',-1.0*valmaxanom,valmaxanom+valdiffanom,valdiffanom,0,'n','mm/month')
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_HadCM3.eps' 
    else:
        plotdata(pliop2_anom,99,lon,lat,'a) HadGEM2 evaporation anomaly',-1.0*valmaxanom,valmaxanom+valdiffanom,valdiffanom,0,'n','mm/month')
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_HadGEM2.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()


    # Pliocene - preindustrial percentage change

    print('j5')

    print(pi_evap_ann)
    plio_peranom=((plio_evap_ann-pi_evap_ann)/pi_evap_ann)*100.
    plotdata(plio_peranom,0,lon,lat,titlediff,-70,78,5,0,'a','%')

    if HadCM3 != 'y':
        plio_peranomp2=((plio_evapp2_ann-pi_evap_ann)/pi_evap_ann)*100.
        plotdata(plio_peranomp2,99,lon,lat,'c) Plio - PI Eanom_HG2+2 %',-70,78,5,0,'a','%')
    

    if HadCM3 == 'y':
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent_HadCM3.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent_HadCM3.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  
    else:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  



    plt.close()


    # big plot for paper
    print('j5a')

    if HadCM3 != 'y':
        V=0
        plotdata(plio_peranomp2,99,lon,lat,'c) HadGEM2: mPWP-PI Evap anomaly',-70.,75.,5.,V,'a','%')

        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_HadGEM2_pcent_single.eps' 
        plt.savefig(fileout, bbox_inches='tight')  

        plt.close()

    else:

        V=0
        plotdata(plio_peranom,99,lon,lat,'d) HadCM3: mPWP-PI Evap anomaly',-70.,75.,5.,V,'a','%')

        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_HadCM3_pcent_single.eps' 
        plt.savefig(fileout, bbox_inches='tight')  

        plt.close()



#end def annmean


def seasmean(m1,m2,m3,figureno,seasname,HadCM3):
    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean


 # get land mask and put on correct grid

    if HadCM3 != 'y':
        fm=Dataset('/nfs/hera1/earjcti/um/HadGEM_ancils/qrparm.mask.nc')
    else:
        fm=Dataset('/nfs/hera2/apps/metadata/ancil/preind2/qrparm.mask.nc')
    lsmlon=fm.variables['longitude'][:]
    lsmlat=fm.variables['latitude'][:]
    lsm=fm.variables['lsm'][:]
    lsm=np.squeeze(lsm)
    lsm,lsmlon = shiftgrid(180.,lsm,lsmlon,start=False)
    fm.close()

   
 
 

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
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/pdfiles/xkvjea@pdn[7-9]*'+m1+'.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/pdfiles/xkvjea@pdn[7-9]*'+m2+'.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/pdfiles/xkvjea@pdn[7-9]*'+m3+'.nc')
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
 
    if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
        pi_land=pi_seasevap/lsm
        pi_sea=pi_seasevap / (np.abs(lsm-1.0))                              
    else:
        print('error lon/lat of land sea mask dont match')
        anom_land=plio_anom * lsm
        plotdata(anom_land,99,lon,lat,'a) mPWP temperature anomaly',0,10,1.0,V,'i',degC)
        plt.show()
        sys.exit()

   
    plotdata(pi_land,0,lon,lat,pititle,0,20,1,0.0,'n','mm/month')
    
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
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjf/netcdf/pdfiles/xkvjfa@pdn[7-9]*'+m1+'.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjf/netcdf/pdfiles/xkvjfa@pdn[7-9]*'+m2+'.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjf/netcdf/pdfiles/xkvjfa@pdn[7-9]*'+m3+'.nc')
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
    
    plio_land=plio_seasevap/lsm
    plio_sea=plio_seasevap / (np.abs(lsm-1.0))                              
    
    plotdata(plio_land,1,lon,lat,pliotitle,0,20,1,0.0,'n','mm/month')



     #==============
     # Pliocene+2

    if HadCM3 != 'y':
        fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/pdfiles/xkvjga@pdn[7-9]*'+m1+'.nc')
        fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/pdfiles/xkvjga@pdn[7-9]*'+m2+'.nc')
        fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/pdfiles/xkvjga@pdn[7-9]*'+m3+'.nc')
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
    
 

    # Pliocene+2 - preindustrial

    if HadCM3 != 'y':
        pliop2_anom=pliop2_seasevap-pi_seasevap
        V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]

        plotdata((pliop2_anom/pi_seasevap)*100.,2,lon,lat,'PlioP2 - PI Eanom_HG2',-100,100,1,0,'a','%')
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_'+seasname+'anom.eps' 
        plt.savefig(fileout, bbox_inches='tight')  

        plt.close()
        plotdata(pliop2_anom,99,lon,lat,'PlioP2 - PI Eanom_HG2',-20,22,2,0,'a','mm/month')
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_'+seasname+'anom_mm.eps' 
        plt.savefig(fileout, bbox_inches='tight')  

        plt.close()
      

    else:
        pliop_anom=pliop_seasevap-pi_seasevap
        V=[-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64]

        plotdata(pliop_anom,3,lon,lat,'Plio - PI Eanom_HG2',-10,10,1,0,'a','mm/month')
        plotdata((pliop_anom/pi_seasevap)*100.,2,lon,lat,'Plio - PI Eanom_HG2',-100,100,1,0,'a','%')
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_'+seasname+'anom_HadCM3.eps' 
        plt.savefig(fileout,bbox_inches='tight')
        plt.close()

       
  
    # Pliocene - preindustrial percentage change

    plio_peranom=((plio_seasevap-pi_seasevap)/pi_seasevap)*100.
    plotdata(plio_peranom,0,lon,lat,difftitle,-50,50,5,0,'a','%')


    if HadCM3 !='y':
        plio_peranomp2=((pliop2_seasevap-pi_seasevap)/pi_seasevap)*100.
        titlename='Plio+2-PI evap %'+seasname
        plotdata(plio_peranomp2,1,lon,lat,titlename,-50,50,5,0,'a','%')

    if HadCM3 =='y':
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent_'+seasname+'HadCM3.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent_'+seasname+'HadCM3.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  
    else:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_anom_percent_'+seasname+'.tiff' 
        plt.savefig(fileout, bbox_inches='tight')  


    plt.close()




#end def seasmean


def evap_chg_by_latitude(HadCM3,land_ocn_ind,abs_pcent):


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
    else:
    # read in data from multiple files
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvje/netcdf/xkvjea@pdn[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]
    
    f.close()
    aevap=np.squeeze(aevap)
    ntimes,ny,nx=np.shape(aevap)
    print(ntimes,ny,nx)

#average across the time dimension and the latitude dimension
    pi_evap_ann=np.mean(aevap,axis=0)
    pi_evap_ann_orig=pi_evap_ann

    if land_ocn_ind == 'l': # mask out all non land points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            pi_evap_ann=pi_evap_ann / (lsm)
            pi_evap_ann[pi_evap_ann == float('Inf')] = float('NaN')
            pi_evap_ann[pi_evap_ann == float('-Inf')] = float('NaN')
 
        else:
            print('error lon/lat of land sea mask dont match')

    if land_ocn_ind == 'o': # mask out all non ocean points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            pi_evap_ann=pi_evap_ann / (np.abs(lsm-1.0))
            pi_evap_ann[pi_evap_ann == float('Inf')] = float('NaN')
            pi_evap_ann[pi_evap_ann == float('-Inf')] = float('NaN')
 
        else:
            print('error lon/lat of land sea mask dont match')

    pi_evap_lat_ann=np.nanmean(pi_evap_ann,axis=1)
  
    
    
     #==============
     # Pliocene

    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy[5-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aevap=f.variables['field184'][:]

    else:
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjf/netcdf/xkvjfa@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]

    f.close()
    aevap=np.squeeze(aevap)
    ntimes,ny,nx=np.shape(aevap)

    plio_evap_ann=np.mean(aevap,axis=0)
    plio_evap_ann_orig=plio_evap_ann


    if land_ocn_ind == 'l': # mask out all non land points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            plio_evap_ann=plio_evap_ann / (lsm)
            plio_evap_ann[plio_evap_ann == float('Inf')] = float('NaN')
            plio_evap_ann[plio_evap_ann == float('-Inf')] = float('NaN')

        else:
            print('error lon/lat of land sea mask dont match')

    if land_ocn_ind == 'o': # mask out all non ocean points.
        if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
            plio_evap_ann=plio_evap_ann / (np.abs(lsm-1.0))
            plio_evap_ann[plio_evap_ann == float('Inf')] = float('NaN')
            plio_evap_ann[plio_evap_ann == float('-Inf')] = float('NaN')

        else:
            print('error lon/lat of land sea mask dont match')




    # find average
    weightarr=np.zeros(np.shape(plio_evap_ann))
    for i in range (0,len(lon)):
        weightarr[:,i]=np.cos(np.deg2rad(lsmlat))

    

    if land_ocn_ind == 'l':
        weightarr=weightarr * lsm
    if land_ocn_ind == 'o':
        weightarr=weightarr * np.abs(lsm-1.0)

    # set values over poles to be 0
    weightarr[0,:]=0.0
    weightarr[len(lat)-1,:]=0.0
    pi_evap_ann[len(lat)-1,:]=0.0
    plio_evap_ann[len(lat)-1,:]=0.0

    

    print('HadCM3 is',HadCM3,'land ocean ind is ',land_ocn_ind,'mean Plioanom ',np.average(plio_evap_ann_orig - pi_evap_ann_orig,weights=weightarr)*60.*60.*24)


    plio_evap_lat_ann=np.nanmean(plio_evap_ann,axis=1)


     #==============
     # Pliocene+2

    if HadCM3 != 'y':
        # read in data from multiple files
        f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/xkvjg/netcdf/xkvjga@pdn[5-9]*.nc')
        aevap=f.variables['field184'][:]
        aevap=np.squeeze(aevap)
        ntimes,ny,nx=np.shape(aevap)
        print(ntimes,ny,nx)
        
    #average across the time dimension
        plio_evapp2_ann=np.mean(aevap,axis=0)
        plio_evapp2_ann_orig=plio_evapp2_ann
        if land_ocn_ind == 'l': # mask out all non land points.
            if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
                plio_evapp2_ann=plio_evapp2_ann / (lsm)
                plio_evapp2_ann[plio_evapp2_ann == float('Inf')] = float('NaN')
            
            else:
                print('error lon/lat of land sea mask dont match')

        if land_ocn_ind == 'o': # mask out all non ocean points.
            if (np.array_equal(lsmlon,lon)) and (np.array_equal(lsmlat,lat)):
                plio_evapp2_ann=plio_evapp2_ann / (np.abs(lsm-1.0))
                plio_evapp2_ann[plio_evapp2_ann == float('Inf')] = float('NaN')

            else:
                print('error lon/lat of land sea mask dont match')


        # find average using weighting array found earlier
        print('HadCM3 is',HadCM3,'land ocean ind is ',land_ocn_ind,'mean Plioanomp2 ',np.average(plio_evapp2_ann_orig-pi_evap_ann_orig,weights=weightarr)*60.*60.*24)

        plio_evapp2_lat_ann=np.nanmean(plio_evapp2_ann,axis=1)
        f.close()


    # Pliocene - preindustrial

    if abs_pcent == 'a':
        plio_anom=(plio_evap_lat_ann-pi_evap_lat_ann)* 60. *  60. *24.

    if abs_pcent == 'p':
        plio_anom=(plio_evap_lat_ann-pi_evap_lat_ann)/ pi_evap_lat_ann
        plio_anom=plio_anom * 100.

    # Pliocene+2 - preindustrial

    if HadCM3 != 'y':
        if abs_pcent == 'a':
          pliop2_anom=(plio_evapp2_lat_ann-pi_evap_lat_ann)* 60.*  60. *24.
         
        if abs_pcent == 'p':
          pliop2_anom=(plio_evapp2_lat_ann-pi_evap_lat_ann)/ pi_evap_lat_ann
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

HadCM3='n'
land_ocn_ind='b'
#plt.figure(figureno)
annmean(figureno,HadCM3,land_ocn_ind)
#HadCM3='n'
#annmean(figureno,HadCM3,land_ocn_ind)
#figureno=figureno+1

#djf mean
#plt.figure(figureno)
#seasmean('dc','ja','fb',figureno,'djf',HadCM3)
#figureno=figureno+1

#jja mean
#plt.figure(figureno)
#seasmean('jn','jl','ag',figureno,'jja',HadCM3)
#figureno=figureno+1

#seasmean('mr','ar','my',figureno,'mam',HadCM3)
#seasmean('sp','ot','nv',figureno,'son',HadCM3)

# to see what the polar amplification is
#HadCM3='y'
#land_ocn_ind='o'    # valid values are l - land , o- ocean b=both
#abs_pcent='a'       # valid values are p - percentage chg, a=absolute change
#evap_ret=evap_chg_by_latitude(HadCM3,land_ocn_ind,abs_pcent)
#lats_HadCM3=evap_ret[0]
#evap_anom_HadCM3=evap_ret[1]
#plt.plot(evap_anom_HadCM3,lats_HadCM3)

#HadCM3='n'
#evap_ret=evap_chg_by_latitude(HadCM3,land_ocn_ind,abs_pcent)
#lats_HadGEM=evap_ret[0]
#evap_anom_HadGEM=evap_ret[1]
#evap_anom_HadGEMp2=evap_ret[2]

#mp.rcParams.update({'font.size':15})
#plt.plot(evap_anom_HadGEMp2,lats_HadGEM,'g',label='HadGEM2')
#plt.plot(evap_anom_HadCM3,lats_HadCM3,'r',label='HadCM3')
#plt.plot(evap_anom_HadGEM,lats_HadGEM,'b')

#if abs_pcent == 'a':
#    plt.xlabel('mm/day',fontsize=20)
#if abs_pcent == 'p':
#    plt.xlabel('%',fontsize=20)
#plt.ylabel('latitude',fontsize=20)

#if land_ocn_ind =='l':
#    plt.title('f) mPWP - preind land evap anom',loc='left',fontsize=25)
#if land_ocn_ind =='o':
#    plt.title('g) mPWP - preind ocean evap anom',loc='left',fontsize=25)
#if land_ocn_ind =='b':
#    plt.title('e) mPWP - preind evap anom',loc='left',fontsize=25)
#axes=plt.gca()
#if abs_pcent == 'p':
#    axes.set_xlim(xmin=-20.0,xmax=100.0)
#axes.set_ylim([-80,80])
#legend=plt.legend(loc='center right')
#if land_ocn_ind =='l':
#    if abs_pcent =='a':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_polar_amp_land.eps' 
#    if abs_pcent =='p':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_polar_amp_land_pcent.eps' 
#if land_ocn_ind =='o':
#    if abs_pcent =='a':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_polar_amp_ocean.eps' 
#    if abs_pcent =='p':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfeva#p/evap_polar_amp_ocean_pcent.eps' 

#if land_ocn_ind =='b':
#    if abs_pcent == 'a':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_polar_amp.eps' 
#    if abs_pcent == 'p':
#        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_surfevap/evap_polar_amp_pcent.eps' 
#plt.savefig(fileout, bbox_inches='tight')  



sys.exit(0)

####

