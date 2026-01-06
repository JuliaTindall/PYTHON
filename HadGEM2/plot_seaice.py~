#!/usr/bin/env python2.7
#NAME
#    PLOT_SEAICE
#PURPOSE
#    This program will plot the sea ice change mPWP-PI in the simulations
#
# search for 'main program' to find end of functions
# Julia 28/11/2018


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
def plotdata(region,plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,oplot_pi,pi_seaice):

    if region == 'NP':
        proj='npstere'
        latbb=45.
    if region == 'SP':
        proj='spstere'
        latbb=-45.
    #proj='stere'

    lons, lats = np.meshgrid(lon,lat)
    if fileno != 99:
        plt.subplot(2,2,fileno+1)

   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe

    map=Basemap(projection=proj,resolution='c',lon_0=0,boundinglat=latbb,round=True,lat_0=90)
    #map=Basemap(projection='stere', lat_0=90., lon_0=0.,
    #                       llcrnrlon=-180, llcrnrlat=60.,
    #                       urcrnrlon=180., urcrnrlat=60.)


    #map.drawmapboundary(fill_color='green')
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
                    print(V)
                    cs = map.contourf(x,y,plotdata,V,norm=mp.colors.LogNorm(vmin=0,vmax=32),cmap='Reds')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    mycmap=mp.cm.get_cmap('Reds',len(V+2))
                    newcolors=mycmap(np.linspace(0,1,len(V+2)))
                    white=([1,1,1,1])
                    newcolors[0:2,:]=white
                    mycmap=ListedColormap(newcolors)
                    cs = map.contourf(x,y,plotdata,V,cmap=mycmap)
                   # map.drawparallels(np.arange(-80.,81.,20.))
                   # map.drawmeridians(np.arange(-180.,181.,20.))
                    map.fillcontinents()
                  
                    cbar = plt.colorbar(cs,orientation="horizontal")
                    # overplot in contours where preindustrial sea ice was over 50% (or 0.5)
                    if oplot_pi == 'y':
                        map.contour(x,y,pi_seaice,[0,0.5],colors='lime',linestyles=[':','--'],linewidth=10)
  


    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        cbar.set_label(cbarname,labelpad=-70,size=20)
        cbar.ax.tick_params(labelsize=20)
        plt.title(titlename,loc='left',fontsize=20)
   


#end def plotdata

def seasmean(m1,m2,m3,figureno,seasname,preind_expt,plio_expt,pliop2_expt,extra):
    # m1 m2 m3 are the month neames needed to reproduce the seasonal mean
    #==============
    # preindustrial

    # read in temperature from a single file in order to get land mask
    f=Dataset('/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/netcdf/pffiles/'+preind_expt+'o@pf'+extra+'76ja.nc')
    temp=f.variables['temp'][:]
    mask=temp/temp # ie temp is 1 everywhere except where it is masked
    mask=np.squeeze(mask)
    f.close()
   

   
    # read in data from multiple files
    fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/netcdf/pffiles/'+preind_expt+'o@pf'+extra+'*'+m1+'.nc')
    fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/netcdf/pffiles/'+preind_expt+'o@pf'+extra+'*'+m2+'.nc')
    fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+preind_expt+'/netcdf/pffiles/'+preind_expt+'o@pf'+extra+'*'+m3+'.nc')
    lat = fa.variables['latitude'][:]
    lon = fa.variables['longitude'][:]
    atemp=fa.variables['iceconc'][:]
    btemp=fb.variables['iceconc'][:]
    ctemp=fc.variables['iceconc'][:]
    atemp=np.squeeze(atemp)
    btemp=np.squeeze(btemp)
    ctemp=np.squeeze(ctemp)
    ntimes,ny,nx=np.shape(atemp)
    
    #average across the time dimension
    pi_atemp_avg=np.mean(atemp,axis=0)
    pi_btemp_avg=np.mean(btemp,axis=0)
    pi_ctemp_avg=np.mean(ctemp,axis=0)

    #stdev across the time dimension
    pi_atemp_stdev=np.std(atemp,axis=0)
    pi_btemp_stdev=np.std(btemp,axis=0)
    pi_ctemp_stdev=np.std(ctemp,axis=0)
    
    pi_seaice=np.mean((pi_atemp_avg,pi_btemp_avg,pi_ctemp_avg),axis=0)
   
    
    pi_seaice=pi_seaice * mask
    
  
    plotdata('NP',pi_seaice,0,lon,lat,'PI HadGEM2: '+seasname,0,1.1,0.1,0,'n','fraction','y',pi_seaice)
    plotdata('SP',pi_seaice,1,lon,lat,'PI HadGEM2: '+seasname,0,1.1,0.1,0,'n','fraction','y',pi_seaice)
    

    fa.close()
    fb.close()
    fc.close()
   
     #==============
     # Pliocene+2


    fa=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/netcdf/pffiles/'+pliop2_expt+'o@pf'+extra+'*'+m1+'.nc')
    fb=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/netcdf/pffiles/'+pliop2_expt+'o@pf'+extra+'*'+m2+'.nc')
    fc=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+pliop2_expt+'/netcdf/pffiles/'+pliop2_expt+'o@pf'+extra+'*'+m3+'.nc')
    atemp=fa.variables['iceconc'][:]
    btemp=fb.variables['iceconc'][:]
    ctemp=fc.variables['iceconc'][:]
    atemp=np.squeeze(atemp)
    btemp=np.squeeze(btemp)
    ctemp=np.squeeze(ctemp)
    
    pliop2_atemp_avg=np.mean(atemp,axis=0)
    pliop2_btemp_avg=np.mean(btemp,axis=0)
    pliop2_ctemp_avg=np.mean(ctemp,axis=0)
    
    pliop2_seaice=np.mean((pliop2_atemp_avg,pliop2_btemp_avg,pliop2_ctemp_avg),axis=0)

    
    fa.close()
    fb.close()
    fc.close()
   


    # plot data
    
    #pliop2_seaice=pliop2_seaice * mask

 
    plotdata('NP',pliop2_seaice,2,lon,lat,'Plio HadGEM2: '+seasname,0,1.1,0.1,0,'n','fraction','n',pi_seaice)
    plotdata('SP',pliop2_seaice,3,lon,lat,'Plio HadGEM2: '+seasname,0,1.1,0.1,0,'n','fraction','n',pi_seaice)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/'+seasname+'_'+pliop2_expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()
    


    # Pliocene - preindustrial

    perc_lost=((pi_seaice-pliop2_seaice)/pi_seaice)*100.
    plotdata('NP',perc_lost,99,lon,lat,'Sea ice loss:'+seasname,0,105.,5.0,0,'n','%','y',pi_seaice)
    
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/'+seasname+'_'+pliop2_expt+'_NP_iceanom.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()
  

    plotdata('SP',perc_lost,99,lon,lat,'Sea ice loss:'+seasname,0,105.,5.0,0,'n','%','y',pi_seaice)
    
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/'+seasname+'_'+pliop2_expt+'_SP_iceanom.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()
    
    plt.show()
    


######################################################   
def seascyc(preind_expt,pliop2_expt,extra,HadCM3):
# get seasonal cycle of sea ice

    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']

    if HadCM3=='y':
        filestart='/nfs/hera1/earjcti/um/'
        filemid='/netcdf/'
    else:
        filestart='/nfs/hera1/earjcti/um/HadGEM_data/'
        filemid='/netcdf/pffiles/'

    # read in temperature from a single file in order to get land mask
    f=Dataset(filestart+preind_expt+filemid+preind_expt+'o@pf'+extra+'76ja.nc')
    temp=f.variables['temp'][:]
    mask=temp/temp # ie temp is 1 everywhere except where it is masked
   
    if HadCM3 == 'y':
      mask=np.ma.array(temp,mask=temp > 1E10)
      mask=mask/mask

    mask=np.squeeze(mask)
    f.close()
   
   
    # read in data from pi and pliocene files
    for mon in range(0,len(monthnames)):
        print(mon)
        fa=MFDataset(filestart+preind_expt+filemid+preind_expt+'o@pf'+extra+'[7-9]*'+monthnames[mon]+'.nc')
        lat = fa.variables['latitude'][:]
        lon = fa.variables['longitude'][:]
        atemp=fa.variables['iceconc'][:]
        atemp=np.squeeze(atemp)
        htemp=fa.variables['icedepth'][:]
        htemp=np.squeeze(htemp)
        fa.close()
        print('array size pi',np.shape(atemp))
        pi_atemp_avg=np.mean(atemp,axis=0)
        pi_htemp_avg=np.mean(htemp,axis=0)
      

        fa=MFDataset(filestart+pliop2_expt+filemid+pliop2_expt+'o@pf'+extra+'[7-9]*'+monthnames[mon]+'.nc')
        atemp=fa.variables['iceconc'][:]
        atemp=np.squeeze(atemp)
        htemp=fa.variables['icedepth'][:]
        htemp=np.squeeze(htemp)
       
        print('array size plio',np.shape(atemp))
        pliop2_atemp_avg=np.mean(atemp,axis=0)
        pliop2_htemp_avg=np.mean(htemp,axis=0)
       
        fa.close()
        print('gotdata',mon)
  

          
        ny,nx=np.shape(pi_atemp_avg)
        if mon == 0:
            pi_seaice=np.zeros((len(monthnames),ny,nx))
            plio_seaice=np.zeros((len(monthnames),ny,nx))
            pi_seaice_d=np.zeros((len(monthnames),ny,nx))
            plio_seaice_d=np.zeros((len(monthnames),ny,nx))
   
        pi_seaice[mon,:,:]=pi_atemp_avg
        plio_seaice[mon,:,:]=pliop2_atemp_avg
        pi_seaice_d[mon,:,:]=pi_htemp_avg
        plio_seaice_d[mon,:,:]=pliop2_htemp_avg
        
   
    #======================================================
    # get average area of seaice

    xres=lon[1]-lon[0]
    yres=lat[1]-lat[0]
    a=40075. # circumference of earth in km
    onedeg=a/360.
    gridbox_nonweight=xres * yres * onedeg * onedeg
  
    pi_avg_ice_nh=np.zeros(len(monthnames))
    plio_avg_ice_nh=np.zeros(len(monthnames))
    pi_avg_ice_sh=np.zeros(len(monthnames))
    plio_avg_ice_sh=np.zeros(len(monthnames))

    pi_vol_ice_nh=np.zeros(len(monthnames))
    plio_vol_ice_nh=np.zeros(len(monthnames))
    pi_vol_ice_sh=np.zeros(len(monthnames))
    plio_vol_ice_sh=np.zeros(len(monthnames))

    for mon in range(0,len(monthnames)):
        print('averaging for',mon)
        for j in range(0,ny):
            coslat=np.cos(np.radians(lat[j]))
            if lat[j] > 0 :
                for i in range(0,nx):
                    pi_avg_ice_nh[mon]=(pi_avg_ice_nh[mon] + 
                        (pi_seaice[mon,j,i] * coslat * gridbox_nonweight))
                    plio_avg_ice_nh[mon]=(plio_avg_ice_nh[mon] + 
                        (plio_seaice[mon,j,i] * coslat * gridbox_nonweight))

                    pi_vol_ice_nh[mon]=(pi_vol_ice_nh[mon] + 
                       (pi_seaice_d[mon,j,i] * coslat * gridbox_nonweight))
                    plio_vol_ice_nh[mon]=(plio_vol_ice_nh[mon] + 
                        (plio_seaice_d[mon,j,i] * coslat * gridbox_nonweight))

            if lat[j] < 0 :
                for i in range(0,nx):
                    pi_avg_ice_sh[mon]=(pi_avg_ice_sh[mon] + 
                        (pi_seaice[mon,j,i] * coslat * gridbox_nonweight))
                    plio_avg_ice_sh[mon]=(plio_avg_ice_sh[mon] + 
                        (plio_seaice[mon,j,i] * coslat * gridbox_nonweight))

                    pi_vol_ice_sh[mon]=(pi_vol_ice_sh[mon] + 
                       (pi_seaice_d[mon,j,i] * coslat * gridbox_nonweight))
                    plio_vol_ice_sh[mon]=(plio_vol_ice_sh[mon] + 
                        (plio_seaice_d[mon,j,i] * coslat * gridbox_nonweight))

       
    
    plt.plot(np.arange(1,13),pi_avg_ice_nh,label='PI')
    plt.plot(np.arange(1,13),plio_avg_ice_nh,label='mPWP')
    plt.title('Arctic Sea Ice - Areal extent',fontsize=20)
    plt.ylabel('km^2',fontsize=20)
    plt.xlabel('month',fontsize=20)
    plt.legend(fontsize=20)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/'+pliop2_expt+'_NH_areal_extent.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()
 

    plt.plot(np.arange(1,13),pi_avg_ice_sh,label='PI')
    plt.plot(np.arange(1,13),plio_avg_ice_sh,label='mPWP')
    plt.title('Antarctic Sea Ice - Areal extent')
    plt.ylabel('km^2')
    plt.xlabel('month')
    plt.legend()

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/'+pliop2_expt+'_SH_areal_extent.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()
  
   
    retdata=[pi_avg_ice_nh,plio_avg_ice_nh,pi_avg_ice_sh,plio_avg_ice_sh,pi_vol_ice_nh,plio_vol_ice_nh,pi_vol_ice_sh,plio_vol_ice_sh]
    return(retdata)

# end def seascyc








################################
# main program

# annual mean
figureno=0

#preind_expt='xkvje'
#plio_expt='xkvjg'
#pliop2_expt='xkvjg'
#extra='n'
#HadCM3='n'

preind_expt='xiboi'
plio_expt='xibol'
pliop2_expt='xibol'
extra='y'
HadCM3='y'



#djf mean
#plt.figure(figureno)
#seasmean('dc','ja','fb',figureno,'djf',preind_expt,plio_expt,pliop2_expt,extra)

#mam mean
#plt.figure(figureno)
#seasmean('mr','ar','my',figureno,'mam',preind_expt,plio_expt,pliop2_expt,extra)
#figureno=figureno+1

#jja mean
#plt.figure(figureno)
#seasmean('jn','jl','ag',figureno,'jja',preind_expt,plio_expt,pliop2_expt,extra)
#figureno=figureno+1

#son mean
#plt.figure(figureno)
#seasmean('sp','ot','nv',figureno,'son',preind_expt,plio_expt,pliop2_expt,extra)
#figureno=figureno+1

#####################################
# plot annual cycle of sea ice loss

preind_expt='xkvje'
plio_expt='xkvjg'
pliop2_expt='xkvjg'
extra='n'
HadCM3='n'
retdata=seascyc(preind_expt,pliop2_expt,extra,HadCM3)
pi_ice_nh_HadGEM=retdata[0]
plio_ice_nh_HadGEM=retdata[1]
pi_ice_sh_HadGEM=retdata[2]
plio_ice_sh_HadGEM=retdata[3]
pi_vol_nh_HadGEM=retdata[4]
plio_vol_nh_HadGEM=retdata[5]
pi_vol_sh_HadGEM=retdata[6]
plio_vol_sh_HadGEM=retdata[7]
print('got HadGEM data') 

preind_expt='xiboi'
plio_expt='xibol'
pliop2_expt='xibol'
extra='y'
HadCM3='y'
retdata=seascyc(preind_expt,pliop2_expt,extra,HadCM3)
pi_ice_nh_HadCM3=retdata[0]
plio_ice_nh_HadCM3=retdata[1]
pi_ice_sh_HadCM3=retdata[2]
plio_ice_sh_HadCM3=retdata[3]
pi_vol_nh_HadCM3=retdata[4]
plio_vol_nh_HadCM3=retdata[5]
pi_vol_sh_HadCM3=retdata[6]
plio_vol_sh_HadCM3=retdata[7]
print('got HadCM3 data') 

# HadCM3 does not agree with fergus' paper where it is 10million
# here we have 7.6million.

print('Arctic sea ice extent in HadGEM',np.mean(plio_ice_nh_HadGEM))
print('Arctic sea ice extent in HadCM3',np.mean(plio_ice_nh_HadCM3))
 
# sea ice areal extent SH

fig=plt.figure()
ax=plt.subplot(111)
ax.plot(np.arange(1,13),pi_ice_sh_HadGEM/1000000.,'b',label='PI HadGEM')
ax.plot(np.arange(1,13),plio_ice_sh_HadGEM/1000000.,'r',label='mPWP HadGEM')
ax.plot(np.arange(1,13),pi_ice_sh_HadCM3/1000000,'b--',label='PI HadCM3')
ax.plot(np.arange(1,13),plio_ice_sh_HadCM3/1000000,'r--',label='mPWP HadCM3')
plt.title('f) Antarctic Sea Ice - Areal Extent',fontsize=20,loc='left')
plt.ylabel('million km^2',fontsize=15)
plt.xlabel('month',fontsize=15)
plt.legend(fontsize=15)
plt.tick_params(axis='both',labelsize=15)
box=ax.get_position()
ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
ax.legend(loc='lower center',bbox_to_anchor=(0.5,-0.3),ncol=4)

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_SH_areal_extent.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_SH_areal_extent.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()


# sea ice vol extent SH

plt.plot(np.arange(1,13),pi_vol_sh_HadGEM,'b',label='PI HadGEM')
plt.plot(np.arange(1,13),plio_vol_sh_HadGEM,'r',label='mPWP HadGEM')
plt.plot(np.arange(1,13),pi_vol_sh_HadCM3,'b--',label='PI HadCM3')
plt.plot(np.arange(1,13),plio_vol_sh_HadCM3,'r--',label='mPWP HadCM3')
plt.title('Antarctic Sea Ice - volume')
plt.ylabel('m^3')
plt.xlabel('month')
plt.legend()

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_SH_vol.eps' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

#sea ice loss SH

plt.plot(np.arange(1,13),pi_ice_sh_HadGEM-plio_ice_sh_HadGEM,label='HadGEM')
plt.plot(np.arange(1,13),pi_ice_sh_HadCM3-plio_ice_sh_HadCM3,label='HadCM3')
plt.title('Antarctic Sea Ice - Loss')
plt.ylabel('km^2')
plt.xlabel('month')
plt.legend()

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_SH_loss.eps' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# NH sea ice extent
fig=plt.figure()
ax=plt.subplot(111)
ax.plot(np.arange(1,13),pi_ice_nh_HadGEM/1000000.,color='b',label='PI HadGEM')
ax.plot(np.arange(1,13),plio_ice_nh_HadGEM/1000000.,color='r',label='mPWP HadGEM')
ax.plot(np.arange(1,13),pi_ice_nh_HadCM3/1000000.,'b--',label='PI HadCM3')
ax.plot(np.arange(1,13),plio_ice_nh_HadCM3/1000000.,'r--',label='mPWP HadCM3')
plt.title('e) Arctic Sea Ice - Areal Extent',fontsize=20,loc='left')
plt.ylabel('million km^2',fontsize=15)
plt.xlabel('month',fontsize=15)
plt.tick_params(axis='both',labelsize=15)
box=ax.get_position()
ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
ax.legend(loc='lower center',bbox_to_anchor=(0.5,-0.3),ncol=4)

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_NH_areal_extent.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_NH_areal_extent.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# sea ice vol extent NH

plt.plot(np.arange(1,13),pi_vol_nh_HadGEM,'b',label='PI HadGEM')
plt.plot(np.arange(1,13),plio_vol_nh_HadGEM,'r',label='mPWP HadGEM')
plt.plot(np.arange(1,13),pi_vol_nh_HadCM3,'b--',label='PI HadCM3')
plt.plot(np.arange(1,13),plio_vol_nh_HadCM3,'r--',label='mPWP HadCM3')
plt.title('e) Arctic Sea Ice - volume')
plt.ylabel('m^3')
plt.xlabel('month')
plt.legend()

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_NH_vol.eps' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

  
#sea ice loss NH

plt.plot(np.arange(1,13),pi_ice_nh_HadGEM-plio_ice_nh_HadGEM,label='HadGEM')
plt.plot(np.arange(1,13),pi_ice_nh_HadCM3-plio_ice_nh_HadCM3,label='HadCM3')
plt.title('Arctic Sea Ice - Loss')
plt.ylabel('km^2')
plt.xlabel('month')
plt.legend()

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_seaice/HadCM3_and_HadGEM_NH_loss.eps' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()


sys.exit(0)

####

