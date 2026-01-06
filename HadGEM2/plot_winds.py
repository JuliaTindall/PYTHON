#!/usr/bin/env python2.7
#NAME
#    PLOT_RADIATION
#PURPOSE
#    This program will plot the radiation budget for the pliocene simulations
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

sys.path.append('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/COMMON')
from jumaps import contourglobe


#functions are:
#  def plotquiver
#  def annmean
#  def seasmean

# functions start here
def plotquiver(udata,vdata,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)

   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary
    x, y = map(lons, lats)
    map.drawcoastlines()
    # quiver plot every nth arrow
    n=5
    qv = map.quiver(x[::n,::n],y[::n,::n],udata[::n,::n],vdata[::n,::n],pivot='mid')
    plt.title(titlename)

#end def plotquiver

########################################################
def annmean_surf(switch,HadCM3,expt,extra):
    # switch is a dummy variable to allow the program to be called

    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt+'/netcdf/'
        filename=expt+'a@pd'+extra+'[7-9]*.nc'
        os.system('ls '+dirname+filename)
        os.system('ls '+dirname+filename+' | wc -l')
        nval=os.system('ls '+dirname+filename+' | wc -l')
        allfiles=subprocess.check_output('ls '+dirname+filename+' | wc -l',shell=True)
        print('number of files=',allfiles)
        f=MFDataset(dirname+filename)
        lat = f.variables['latitude_1'][:]
        lon = f.variables['longitude_1'][:]
        
        # we will plot 10 m winds (called u and v)
        
        u=f.variables['u'][:]
        v=f.variables['v'][:]


    u=np.squeeze(u)
    v=np.squeeze(v)
    ntimes,ny,nx=np.shape(u)
    print(ntimes,ny,nx)
    
#average across the time dimension
    u_ann=np.mean(u,axis=0)
    v_ann=np.mean(v,axis=0)
    
    
    plt.figure(0)
    ms='m/s'
    lontemp=lon
    titlename=expt+' winds'
    u_ann,lon = shiftgrid(180.,u_ann,lon,start=False)    
    lon=lontemp
    v_ann,lon = shiftgrid(180.,v_ann,lon,start=False)    
    plotquiver(u_ann,v_ann,lon,lat,titlename,0,400,40.0,0.0,'n',ms)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/annmean_'+expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


#end def annmean



########################################################
def seasmean_surf(m1,m2,m3,seasname,HadCM3,expt,extra):
 
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt+'/winds_data/'
        filea=expt+'a@pd*'+m1+'_winds.nc'
        fileb=expt+'a@pd*'+m2+'_winds.nc'
        filec=expt+'a@pd*'+m3+'_winds.nc'
        os.system('ls '+dirname+filea+' | wc -l')
        os.system('ls '+dirname+fileb+' | wc -l')
        os.system('ls '+dirname+filec+' | wc -l')
        print(dirname+filea)
        fa=MFDataset(dirname+filea)
        fb=MFDataset(dirname+fileb)
        fc=MFDataset(dirname+fileb)

        lat = fa.variables['latitude_1'][:]
        lon = fa.variables['longitude_1'][:]
        au=fa.variables['u'][:]
        bu=fb.variables['u'][:]
        cu=fc.variables['u'][:]
        au=np.squeeze(au)
        bu=np.squeeze(bu)
        cu=np.squeeze(cu)


        av=fa.variables['v'][:]
        bv=fb.variables['v'][:]
        cv=fc.variables['v'][:]
        av=np.squeeze(av)
        bv=np.squeeze(bv)
        cv=np.squeeze(cv)
        ntimes,ny,nx=np.shape(au)
        print(ntimes,ny,nx)
    
        print(np.shape(au),np.shape(bu),np.shape(cu))
        u=au+bu+cu
        v=av+bv+cv

    u=np.squeeze(u)
    v=np.squeeze(v)
    ntimes,ny,nx=np.shape(u)
    print(ntimes,ny,nx)
    
#average across the time dimension
    u_seas=np.mean(u,axis=0)
    v_seas=np.mean(v,axis=0)
    
    
    plt.figure(0)
    ms='m/s'
    lontemp=lon
    titlename=expt+' winds:'+seasname
    u_seas,lon = shiftgrid(180.,u_seas,lon,start=False)    
    lon=lontemp
    v_seas,lon = shiftgrid(180.,v_seas,lon,start=False)    
    plotquiver(u_seas,v_seas,lon,lat,titlename,0,400,40.0,0.0,'n',ms)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_'+expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()

    retdata=[u_seas, v_seas,lon,lat]

    return retdata

#end def seasmean


########################################################
def seasmean_height(m1,m2,m3,seasname,HadCM3,expt,extra,pressreq):
 
    import numpy as np
    print(pressreq)
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/xiboia@pdy[7-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt+'/netcdf/pcfiles/'
        filea=expt+'a@pc*'+m1+'.nc'
        fileb=expt+'a@pc*'+m2+'.nc'
        filec=expt+'a@pc*'+m3+'.nc'
        os.system('ls '+dirname+filea+' | wc -l')
        os.system('ls '+dirname+fileb+' | wc -l')
        os.system('ls '+dirname+filec+' | wc -l')
        print(dirname+filea)
        fa=MFDataset(dirname+filea)
        fb=MFDataset(dirname+fileb)
        fc=MFDataset(dirname+fileb)

        lat = fa.variables['latitude_1'][:]
        lon = fa.variables['longitude_1'][:]
        pressure=fa.variables['p'][:]
        au=fa.variables['u'][:]
        bu=fb.variables['u'][:]
        cu=fc.variables['u'][:]
        au=np.squeeze(au)
        bu=np.squeeze(bu)
        cu=np.squeeze(cu)


        av=fa.variables['v'][:]
        bv=fb.variables['v'][:]
        cv=fc.variables['v'][:]
        av=np.squeeze(av)
        bv=np.squeeze(bv)
        cv=np.squeeze(cv)
        ntimes,npress,ny,nx=np.shape(au)
      
        for press in range(0,len(pressure)):
            if pressure[press] == pressreq:
                u=au[:,press,:,:]+bu[:,press,:,:]+cu[:,press,:,:]
                v=av[:,press,:,:]+bv[:,press,:,:]+cv[:,press,:,:]

    u=np.squeeze(u)
    v=np.squeeze(v)
    
#average across the time dimension
    u_seas=np.mean(u,axis=0)
    v_seas=np.mean(v,axis=0)
    
    print('j1')
 
    plt.figure(0)
    ms='m/s'
    lontemp=lon
    titlename=expt+' winds:'+seasname
    u_seas,lon = shiftgrid(180.,u_seas,lon,start=False)    
    lon=lontemp
    v_seas,lon = shiftgrid(180.,v_seas,lon,start=False)    
    plotquiver(u_seas,v_seas,lon,lat,titlename,0,400,40.0,0.0,'n',ms)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_'+expt+'_'+np.str(np.int(pressreq))+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()

    retdata=[u_seas, v_seas,lon,lat]

    return retdata

#end def seasmean

##################################
def allseason_surf(exptname,extra,controlname,controlextra,HadCM3):
    #djf mean
    seasname='djf'
    retdata=seasmean('dc','ja','fb',seasname,HadCM3,exptname,extra)
    u_seas=retdata[0]
    v_seas=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # get djfmean from control
    retdata=seasmean('dc','ja','fb',seasname,HadCM3,controlname,controlextra)
    u_ctl=retdata[0]
    v_ctl=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # plot djf anomaly
    titlename=exptname+'-'+controlname+':'+seasname+': wind anomaly'
    plotquiver(u_seas-u_ctl,v_seas-v_ctl,lon,lat,titlename,0,400,40.0,0.0,'n','m/s)')

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_'+exptname+'-'+controlname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


    # do some contour plots of u and v for djf

    titlename='u: '+exptname
    contourglobe(u_seas,3,2,0,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+exptname
    contourglobe(v_seas,3,2,1,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='u: '+controlname
    contourglobe(u_ctl,3,2,2,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+controlname
    contourglobe(v_ctl,3,2,3,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='uanom'
    contourglobe(u_seas-u_ctl,3,2,4,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')
    titlename='vanom'
    contourglobe(v_seas-v_ctl,3,2,5,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_and_'+exptname+'-'+controlname+'u_and_v.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()



    #jja mean
    seasname='jja'
    retdata=seasmean('jn','jl','ag',seasname,HadCM3,exptname,extra)
    u_seas=retdata[0]
    v_seas=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # get jjamean from control
    retdata=seasmean('jn','jl','ag',seasname,HadCM3,controlname,controlextra)
    u_ctl=retdata[0]
    v_ctl=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # plot jja anomaly
    titlename=exptname+'-'+controlname+':'+seasname+': wind anomaly'
    plotquiver(u_seas-u_ctl,v_seas-v_ctl,lon,lat,titlename,0,400,40.0,0.0,'n','m/s)')

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_'+exptname+'-'+controlname+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


    # do some contour plots of u and v for jja

    titlename='u: '+exptname
    contourglobe(u_seas,3,2,0,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+exptname
    contourglobe(v_seas,3,2,1,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='u: '+controlname
    contourglobe(u_ctl,3,2,2,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+controlname
    contourglobe(v_ctl,3,2,3,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='uanom'
    contourglobe(u_seas-u_ctl,3,2,4,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')
    titlename='vanom'
    contourglobe(v_seas-v_ctl,3,2,5,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_and_'+exptname+'-'+controlname+'u_and_v.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()

#


##################################
def allseason_height(exptname,extra,controlname,controlextra,HadCM3,pressreq):
    #djf mean
    seasname='djf'
    retdata=seasmean_height('dc','ja','fb',seasname,HadCM3,exptname,extra,pressreq)
    u_seas=retdata[0]
    v_seas=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # get djfmean from control
    retdata=seasmean_height('dc','ja','fb',seasname,HadCM3,controlname,controlextra,pressreq)
    u_ctl=retdata[0]
    v_ctl=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # plot djf anomaly
    titlename=exptname+'-'+controlname+':'+seasname+': wind anomaly'
    plotquiver(u_seas-u_ctl,v_seas-v_ctl,lon,lat,titlename,0,400,40.0,0.0,'n','m/s)')

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_'+exptname+'-'+controlname+'_'+np.str(np.int(pressreq))+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


    # do some contour plots of u and v for djf

    titlename='u: '+exptname
    contourglobe(u_seas,3,2,0,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+exptname
    contourglobe(v_seas,3,2,1,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='u: '+controlname
    contourglobe(u_ctl,3,2,2,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+controlname
    contourglobe(v_ctl,3,2,3,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='uanom'
    contourglobe(u_seas-u_ctl,3,2,4,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')
    titlename='vanom'
    contourglobe(v_seas-v_ctl,3,2,5,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_and_'+exptname+'-'+controlname+'_'+np.str(np.int(pressreq))+'_u_and_v.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()



    #jja mean
    seasname='jja'
    retdata=seasmean_height('jn','jl','ag',seasname,HadCM3,exptname,extra,pressreq)
    u_seas=retdata[0]
    v_seas=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # get jjamean from control
    retdata=seasmean_height('jn','jl','ag',seasname,HadCM3,controlname,controlextra,pressreq)
    u_ctl=retdata[0]
    v_ctl=retdata[1]
    lon=retdata[2]
    lat=retdata[3]

    # plot jja anomaly
    titlename=exptname+'-'+controlname+':'+seasname+': wind anomaly'
    plotquiver(u_seas-u_ctl,v_seas-v_ctl,lon,lat,titlename,0,400,40.0,0.0,'n','m/s)')

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_'+exptname+'-'+controlname+'_'+np.str(np.int(pressreq))+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


    # do some contour plots of u and v for jja

    titlename='u: '+exptname
    contourglobe(u_seas,3,2,0,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+exptname
    contourglobe(v_seas,3,2,1,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='u: '+controlname
    contourglobe(u_ctl,3,2,2,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='v: '+controlname
    contourglobe(v_ctl,3,2,3,lon,lat,titlename,-20,22,2,0,'n','m/s','c')
    titlename='uanom'
    contourglobe(u_seas-u_ctl,3,2,4,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')
    titlename='vanom'
    contourglobe(v_seas-v_ctl,3,2,5,lon,lat,titlename,-3,3.3,0.3,0,'n','m/s','c')


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/'+seasname+'_and_'+exptname+'-'+controlname+'+'+pressreq+'u_and_v.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()


def zonalmean_ann_surf(exptname,controlname,extra,controlextra,HadCM3):
    # get experiment data
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+'a@pdy[7-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/winds_data/'
        filename=exptname+'a@pd'+extra+'[7-9]*winds.nc'
        os.system('ls '+dirname+filename+' | wc -l')
        f=MFDataset(dirname+filename)
        lat = f.variables['latitude_1'][:]
        lon = f.variables['longitude_1'][:]
        # we will plot zonal wind 10 m winds (called)
        u=f.variables['u'][:]
     

    u=np.squeeze(u)
    ntimes,ny,nx=np.shape(u)
    print(ntimes,ny,nx)
    
    #average across the time dimension and longitude direction
    u_ann=np.mean(u,axis=0)
    u_zm_pi=np.mean(u_ann,axis=1)
   

  # get control data
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/'+controlname+'_netcdf/'+controlname+'a@pdy[7-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+controlname+'/winds_data/'
        filename=controlname+'a@pd'+extra+'[7-9]*winds.nc'
        os.system('ls '+dirname+filename+' | wc -l')
        f=MFDataset(dirname+filename)
        lat = f.variables['latitude_1'][:]
        lon = f.variables['longitude_1'][:]
        # we will plot zonal wind 10 m winds (called)
        u=f.variables['u'][:]
     

    u=np.squeeze(u)
    ntimes,ny,nx=np.shape(u)
    print(ntimes,ny,nx)
    
    #average across the time dimension and longitude direction
    u_ann=np.mean(u,axis=0)
    u_zm_plio=np.mean(u_ann,axis=1)
   
    # plot data

    plt.subplot(1,2,1)
    plt.plot(u_zm_pi,lat,label='preindustrial')
    plt.plot(u_zm_plio,lat,label='mPWP')
    plt.plot([0,0],[-90.,90.])
    plt.title('zonal mean u')
    plt.xlabel('m/s')
    plt.ylabel('latitude')
    plt.legend()
   
    plt.subplot(1,2,2)
    plt.plot(u_zm_plio-u_zm_pi,lat,label='mPWP-PI')
    plt.plot([0,0],[-90.,90.])
    plt.title('zm u_anom')
    plt.xlabel('m/s')
    plt.ylabel('latitude')
    plt.legend()

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/zonal_meanu_'+exptname+'-'+controlname+'_ann.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()
  
    #check for each month individually
    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    for i in range(0,len(monthnames)):
        if HadCM3 == 'y':
            f=MFDataset('/nfs/hera1/earjcti/um/netcdf/'+controlname+'_netcdf/'+controlname+'a@pdy[7-9]*.nc')
            lat = f.variables['latitude'][:]
            lon = f.variables['longitude'][:]
        else:
            #setup filename
            direxpt='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/winds_data/'
            dircntl='/nfs/hera1/earjcti/um/HadGEM_data/'+controlname+'/winds_data/'
            fileexpt=exptname+'a@pd'+extra+'[7-9]*'+monthnames[i]+'_winds.nc'
            filecntl=controlname+'a@pd'+extra+'[7-9]*'+monthnames[i]+'_winds.nc'

            # read data from experiment and control
            f=MFDataset(direxpt+fileexpt)
            lat = f.variables['latitude_1'][:]
            lon = f.variables['longitude_1'][:]
            u_ex=f.variables['u'][:]
            f.close()
            
            f=MFDataset(dircntl+filecntl)
            u_ctl=f.variables['u'][:]
            f.close()

            # average it
            u_ex=np.squeeze(u_ex)
            u_ann=np.mean(u_ex,axis=0)
            u_zm_ex=np.mean(u_ann,axis=1)
            u_ctl=np.squeeze(u_ctl)
            u_ann=np.mean(u_ctl,axis=0)
            u_zm_ct=np.mean(u_ann,axis=1)
  
            # plot it
     
            plt.subplot(1,2,1)
            plt.plot(u_zm_ct,lat,label='preindustrial')
            plt.plot(u_zm_ex,lat,label='mPWP')
            plt.plot([0,0],[-90.,90.])
            plt.title('zonal mean u '+monthnames[i])
            plt.xlabel('m/s')
            plt.ylabel('latitude')
            plt.legend()
   
            plt.subplot(1,2,2)
            plt.plot(u_zm_ex-u_zm_ct,lat,label='mPWP-PI')
            plt.plot([0,0],[-90.,90.])
            plt.title('zm u_anom')
            plt.xlabel('m/s')
            plt.ylabel('latitude')
            plt.legend()

          
            fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/zonal_meanu_'+exptname+'-'+controlname+'_'+monthnames[i]+'.eps' 
            plt.savefig(fileout, bbox_inches='tight')  
    
            plt.close()

    
   


#enddef zonalmean_ann


def zonalmean_ann_height(exptname,controlname,extra,controlextra,HadCM3,pressreq):
    # get experiment data
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/'+exptname+'_netcdf/'+exptname+'a@pcy[7-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/netcdf/pcfiles/'
        filename=exptname+'a@pc'+extra+'[7-9]*.nc'
        os.system('ls '+dirname+filename+' | wc -l')
        f=MFDataset(dirname+filename)
        lat = f.variables['latitude_1'][:]
        lon = f.variables['longitude_1'][:]
        pressure=f.variables['p'][:]
        # we will plot zonal wind 10 m winds (called)
        u=f.variables['u'][:]
    
    ntimes,npress,ny,nx=np.shape(u)
    print(ntimes,npress,ny,nx)
  
    for press in range(0,len(pressure)):
        if pressure[press] == pressreq:
            u_p=u[:,press,:,:]

    u=0
    
    #average across the time dimension and longitude direction
    u_ann=np.mean(u_p,axis=0)
    u_zm_pi=np.mean(u_ann,axis=1)
   

  # get control data
    if HadCM3 == 'y':
        f=MFDataset('/nfs/hera1/earjcti/um/netcdf/'+controlname+'_netcdf/'+controlname+'a@pdy[7-9]*.nc')
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
    else:
    # read in data from multiple files
        dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+controlname+'/netcdf/pcfiles/'
        filename=controlname+'a@pc'+extra+'[7-9]*.nc'
        os.system('ls '+dirname+filename+' | wc -l')
        f=MFDataset(dirname+filename)
        lat = f.variables['latitude_1'][:]
        lon = f.variables['longitude_1'][:]
        # we will plot zonal wind 10 m winds (called)
        u=f.variables['u'][:]
     

    for press in range(0,len(pressure)):
        if pressure[press] == pressreq:
            u_p=u[:,press,:,:]

    u_ann=np.mean(u_p,axis=0)
    u_zm_plio=np.mean(u_ann,axis=1)
   
    # plot data

    plt.subplot(1,2,1)
    plt.plot(u_zm_pi,lat,label='preindustrial')
    plt.plot(u_zm_plio,lat,label='mPWP')
    plt.plot([0,0],[-90.,90.])
    plt.title('zonal mean u')
    plt.xlabel('m/s')
    plt.ylabel('latitude')
    plt.legend()
   
    plt.subplot(1,2,2)
    plt.plot(u_zm_plio-u_zm_pi,lat,label='mPWP-PI')
    plt.plot([0,0],[-90.,90.])
    plt.title('zm u_anom')
    plt.xlabel('m/s')
    plt.ylabel('latitude')
    plt.legend()

   
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/zonal_meanu_'+exptname+'-'+controlname+'_'+np.str(pressreq)+'_ann.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    
    plt.close()
  
    #check for each month individually
    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    for i in range(0,len(monthnames)):
        if HadCM3 == 'y':
            f=MFDataset('/nfs/hera1/earjcti/um/netcdf/'+controlname+'_netcdf/'+controlname+'a@pdy[7-9]*.nc')
            lat = f.variables['latitude'][:]
            lon = f.variables['longitude'][:]
        else:
            #setup filename
            direxpt='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/netcdf/pcfiles/'
            dircntl='/nfs/hera1/earjcti/um/HadGEM_data/'+controlname+'/netcdf/pcfiles/'
            fileexpt=exptname+'a@pc'+extra+'[7-9]*'+monthnames[i]+'.nc'
            filecntl=controlname+'a@pc'+extra+'[7-9]*'+monthnames[i]+'.nc'

            # read data from experiment and control
            f=MFDataset(direxpt+fileexpt)
            lat = f.variables['latitude_1'][:]
            lon = f.variables['longitude_1'][:]
            pressure=f.variables['p'][:]
            u_ex=f.variables['u'][:]
            f.close()
            
            f=MFDataset(dircntl+filecntl)
            u_ctl=f.variables['u'][:]
            f.close()

            for press in range(0,len(pressure)):
                if pressure[press] == pressreq:
                    u_p_ctl=u_ctl[:,press,:,:]
                    u_p_ex=u_ex[:,press,:,:]

            # average it
            u_ex=np.squeeze(u_p_ex)
            u_ann=np.mean(u_ex,axis=0)
            u_zm_ex=np.mean(u_ann,axis=1)
            u_ctl=np.squeeze(u_p_ctl)
            u_ann=np.mean(u_ctl,axis=0)
            u_zm_ct=np.mean(u_ann,axis=1)
  
            # plot it
     
            plt.subplot(1,2,1)
            plt.plot(u_zm_ct,lat,label='preindustrial')
            plt.plot(u_zm_ex,lat,label='mPWP')
            plt.plot([0,0],[-90.,90.])
            plt.title('zonal mean u '+monthnames[i])
            plt.xlabel('m/s')
            plt.ylabel('latitude')
            plt.legend()
   
            plt.subplot(1,2,2)
            plt.plot(u_zm_ex-u_zm_ct,lat,label='mPWP-PI')
            plt.plot([0,0],[-90.,90.])
            plt.title('zm u_anom')
            plt.xlabel('m/s')
            plt.ylabel('latitude')
            plt.legend()

          
            fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_winds/zonal_meanu_'+exptname+'-'+controlname+'_'+monthnames[i]+'_'+np.str(np.int(pressreq))+'.eps' 
            plt.savefig(fileout, bbox_inches='tight')  
    
            plt.close()

    
   


#enddef zonalmean_ann_height


################################
# main program

# annual mean
figureno=0

HadCM3='n'
#exptname='ximuq'
#extra='l'

exptname='xkvjg'  #xkvje xkvjf xkvjg
extra='n'
controlname='xkvje'
controlextra='n'


#allseason_surf(exptname,extra,controlname,controlextra,HadCM3)

#zonalmean_ann_surf(exptname,controlname,extra,controlextra,HadCM3)

pressreq=200.
#zonalmean_ann_height(exptname,controlname,extra,controlextra,HadCM3,pressreq)
allseason_height(exptname,extra,controlname,controlextra,HadCM3,pressreq)






sys.exit(0)

####

