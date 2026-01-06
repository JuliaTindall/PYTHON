#!/usr/bin/env python2.7
#|N|A|ME
#    PLOT_d18o_timeslice
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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import iris
from iris.cube import CubeList
import iris.plot as iplt
import iris.quickplot as qplt
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid



# functions start here
def plotdata(cube,fileno,minval,maxval,valinc,V,uselog,cbarname,
             latmin,latmax,lonmin,lonmax,titlename):
   
    if fileno != 99:
        plt.subplot(2,2,fileno+1)

    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = iplt.contourf(cube,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal")
    else:
        if uselog =='la':
            cs = iplt.contourf(cube,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu_r')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='max')

        else:
            if uselog =='a':
                print(V)
                mycmap=mp.cm.get_cmap('bwr',len(V+2))
                newcolors=mycmap(np.linspace(0,1,len(V+1)))
                print(newcolors)
                white=([1,1,1,1])
                print(len(V),np.int((len(V)+1)/2))
                newcolors[(np.int((len(V)+1)/2)-2)
                          :np.int(((len(V)+1)/2)+1),:]=white
                mycmap=ListedColormap(newcolors)
                print(mycmap)
                
                cs = iplt.contourf(cube,V,cmap=mycmap,extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='i': #increasing
                    print(V)
                    cs = iplt.contourf(cube,V,norm=mp.colors.LogNorm(vmin=0,vmax=32),cmap='Reds')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    cs = iplt.contourf(cube,V,extend='both',cmap='rainbow')
                    cbar = plt.colorbar(cs,orientation="horizontal")


    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        #cbar.set_label(cbarname,labelpad=-70,size=15)
        cbar.set_label(cbarname,size=15)
        cbar.ax.tick_params(labelsize=15)
        plt.title(titlename,loc='left',fontsize=15)
    plt.gca().coastlines()
   


#end def plotdata

def d18o_precip(expt1,extra1,startyear1,endyear1,timeslice1,expt2,extra2,startyear2,endyear2,timeslice2,icevolcorr1,icevolcorr2):

    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']

    # read in data from expt1 files
    filestart='/nfs/hera1/earjcti/um/'+expt1+'/pcpd/'+expt1+'a#pd'+extra1

    cubes = CubeList([])
    for year in range(startyear1,endyear1):
        for month in months:
            filename = filestart + str(year) + month + '+.nc'
            cube = iris.load_cube(filename,'QCL')
            cubes.append(cube)
           
    iris.util.equalise_attributes(cubes)
    allcubes = cubes.concatenate_cube()
    meancube1 = allcubes.collapsed('t',iris.analysis.MEAN)
    data1=meancube1.data
    print(meancube1)
   

    # read in data from expt2 files
    filestart='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'a#pd'+extra2

    cubes = CubeList([])
    for year in range(startyear2,endyear2):
        for month in months:
            filename = filestart + str(year) + month + '+.nc'
            cube = iris.load_cube(filename,'QCL')
            cubes.append(cube)
           
    iris.util.equalise_attributes(cubes)
    allcubes = cubes.concatenate_cube()
    meancube2 = allcubes.collapsed('t',iris.analysis.MEAN)
    data2=meancube2.data
    print(meancube2)

    #######################################
    # calculate d18o from cube1 and cube2 abd okit
    tot_18o1=data1[1,:,:]+data1[4,:,:]+data1[7,:,:]+data1[10,:,:]
    tot_16o1=data1[0,:,:]+data1[3,:,:]+data1[6,:,:]+data1[9,:,:]
    d18o1 = ((tot_18o1 / tot_16o1) - 2005.2E-6) / 2005.2E-9
    d18o1_cube = meancube1[0,:,:].copy(data=d18o1)
    d18o1_cube.long_name='d18o'
    print(d18o1_cube)
    d18o1_cube.units = None

    # plot expt 1
    titlename = 'd18op: ' + expt1
    plotdata(d18o1_cube,99,-50,10,5,0.0,'n','permille',-90,90,-180,180,
             titlename)
    fileout='/nfs/hera1/earjcti/um/' + expt1 + '/isotopes/d18o.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/hera1/earjcti/um/' + expt1 + '/isotopes/d18o.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
    iris.save(d18o1_cube,
              '/nfs/hera1/earjcti/um/' + expt1 + '/isotopes/d18o.nc')
    
    
    tot_18o2=data2[1,:,:]+data2[4,:,:]+data2[7,:,:]+data2[10,:,:]
    tot_16o2=data2[0,:,:]+data2[3,:,:]+data2[6,:,:]+data2[9,:,:]
    d18o2 = ((tot_18o2 / tot_16o2) - 2005.2E-6) / 2005.2E-9
    d18o2_cube = meancube2[0,:,:].copy(data=d18o2)
    d18o2_cube.long_name='d18o'
    d18o2_cube.units = None


    # plot expt 2
    titlename = 'd18op: ' + expt2
    plotdata(d18o2_cube,99,-50,10,5,0.0,'n','permille',-90,90,-180,180,
             titlename)
    fileout='/nfs/hera1/earjcti/um/' + expt2 + '/isotopes/d18o.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/hera1/earjcti/um/' + expt2 + '/isotopes/d18o.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
    iris.save(d18o2_cube,
              '/nfs/hera1/earjcti/um/' + expt2 + '/isotopes/d18o.nc')
   

    # do anomaly
    d18o_anom=d18o2_cube - d18o1_cube
    titlename = 'd18op: ' + expt2 + '-' + expt1
    plotdata(d18o_anom,99,-10,11,1,0.0,'a','permille',-90,90,-180,180,
             titlename)
    fileout=('/nfs/hera1/earjcti/um/' + expt2 + '/isotopes/' + 
             expt2 + '-' + expt1 + 'd18o.eps')
    plt.savefig(fileout,bbox_inches='tight')
    fileout=('/nfs/hera1/earjcti/um/' + expt2 + '/isotopes/' +       
             expt2 + '-' + expt1 + 'd18o.png')
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
    


    filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'a@pd'+extra2+'['+startdec2+'-'+enddec2+']*.nc'
    f=MFDataset(filenames)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['QCL'][:] # get d18o
    f.close()
    atemp_avg=np.mean(atemp,axis=0)  # average by time
    tot_18o=atemp_avg[1,:,:]+atemp_avg[4,:,:]+atemp_avg[7,:,:]+atemp_avg[10,:,:]
    tot_16o=atemp_avg[0,:,:]+atemp_avg[3,:,:]+atemp_avg[6,:,:]+atemp_avg[9,:,:]
    
    d18o_expt2=((tot_18o/tot_16o)-2005.2E-6)/2005.2E-9
    d18o_expt2,lon = shiftgrid(180.,d18o_expt2,lon,start=False)

   
    titlename='d18O_p ('+timeslice2+'-'+timeslice1+')'
    plotdata(d18o_expt2-d18o_expt1,99,lon1,lat,titlename,-3,3.5,0.5,0.0,'a','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt2+'-'+expt1+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18op_'+expt2+'-'+expt1+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
    
 

#end def d18o_precip



def d18o_runoff(expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2):

    # read in data from expt1 files
    filenames='/nfs/hera1/earjcti/um/'+expt1+'/netcdf/'+expt1+'a@pd'+extra1+'['+startdec1+'-'+enddec1+']*.nc'
    print(filenames)
    f=MFDataset(filenames)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    #sruoff_18o=f.variables['slowrunoff_1'][:] 
    #fruoff_18o=f.variables['fastrunoff_1'][:]
    sruoff_16o=f.variables['slowrunoff'][:] 
    fruoff_16o=f.variables['fastrunoff'][:]
    f.close()

    #runoff_18o=np.mean(sruoff_18o,axis=0)+np.mean(fruoff_18o,axis=0)
    runoff_16o=np.mean(sruoff_16o,axis=0)+np.mean(fruoff_16o,axis=0)

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

    filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'a@pd'+extra2+'['+startdec2+'-'+enddec2+']*.nc'
    f=MFDataset(filenames)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    #sruoff_18o=f.variables['slowrunoff_1'][:] 
    #fruoff_18o=f.variables['fastrunoff_1'][:]
    sruoff_16o=f.variables['slowrunoff'][:] 
    fruoff_16o=f.variables['fastrunoff'][:]
    f.close()

    #runoff_18o=np.mean(sruoff_18o,axis=0)+np.mean(fruoff_18o,axis=0)
    runoff_16o=np.mean(sruoff_16o,axis=0)+np.mean(fruoff_16o,axis=0)
    
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
    plotdata(runoff_16o_expt2-runoff_16o_expt1,99,lon1,lat,titlename,-3,3.1,0.1,0.0,'a','mm/day',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/runoff_'+expt2+'-'+expt1+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/runoff_'+expt2+'-'+expt1+'.png'
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



def d18o_ocean(expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2):

    latmin=0.
    latmax=40.
    lonmin=50.
    lonmax=130.

    icevolcorr_diff=icevolcorr2-icevolcorr1
   
    # read in data from expt1 files
    filenames='/nfs/hera1/earjcti/um/'+expt1+'/netcdf/'+expt1+'o@pg'+extra1+'['+startdec1+'-'+enddec1+']*.nc'
    print(filenames)
    f=MFDataset(filenames)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['otracer1'][:] # get d18o
    f.close()
    atemp_avg=np.mean(atemp,axis=0)  # average by time
    atemp_top=atemp_avg[0,:,:]
   
    d18o_expt1=((atemp_top - 2005.2E-6)/2005.2E-9)+icevolcorr1
    d18o_expt1,lon1 = shiftgrid(180.,d18o_expt1,lon,start=False)

    
  
    # read in data from expt2 files

    filenames='/nfs/hera1/earjcti/um/'+expt2+'/netcdf/'+expt2+'o@pg'+extra2+'['+startdec2+'-'+enddec2+']*.nc'
    print(filenames)
    f=MFDataset(filenames)
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['otracer1'][:] # get d18o
    f.close()
    atemp_avg=np.mean(atemp,axis=0)  # average by time
    atemp_top=atemp_avg[0,:,:]
   
    d18o_expt2=((atemp_top - 2005.2E-6)/2005.2E-9)+icevolcorr2
    d18o_expt2,lon = shiftgrid(180.,d18o_expt2,lon,start=False)

    latmin=0.
    latmax=40.
    lonmin=50.
    lonmax=130.

    titlename='d18O_sw ('+timeslice2+'-'+timeslice1+')'
    plotdata(d18o_expt2-d18o_expt1,99,lon1,lat,titlename,-0.5+icevolcorr_diff,0.55+icevolcorr_diff,0.05,0.0,'a','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'-'+expt1+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'-'+expt1+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
  
    titlename='d18O_sw ('+timeslice2+')'
    plotdata(d18o_expt2,0,lon1,lat,titlename,-1.0,1.1,0.1,0.0,'n','permille',latmin,latmax,lonmin,lonmax)
    titlename='d18O_sw ('+timeslice1+')'
    plotdata(d18o_expt1,1,lon1,lat,titlename,-1.0,1.1,0.1,0.0,'n','permille',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'_and_'+expt1+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/d18o_ocn_'+expt2+'_and_'+expt1+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
    
 

#end def d18o_ocean


def salinity_ocean(expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2):

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
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'-'+expt1+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'-'+expt1+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
  
    titlename='sal_sw ('+timeslice2+')'
    plotdata(sal_expt2,0,lon1,lat,titlename,30.0,40.0,1.0,0.0,'n','psu',latmin,latmax,lonmin,lonmax)
    titlename='sal_sw ('+timeslice1+')'
    plotdata(sal_expt1,1,lon1,lat,titlename,30.0,40.0,1.0,0.0,'n','psu',latmin,latmax,lonmin,lonmax)
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'_and_'+expt1+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/d18o_timeslice/sal_ocn_'+expt2+'_and_'+expt1+'.png'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()
   
    
 

#end def salinity_ocean

################################
# main program

figureno=0

# timeslices are xiboi=preindustrial, xibol=3205 - km5c', xoekc=3205-km5c, xoekd=3060 (K1), xoeke=2950 (G17), xoekf=3155 (KM3)
#expt1='xiboi'
#timeslice1='preindustrial'
#extra1='y'
#startdec1='5'
#enddec1='9'
#icevolcorr1=0.0


#expt2='xibol'
#extra2='y'
#startdec2='5'
#enddec2='9'
#timeslice2='KM5C'
#icevolcorr2=-0.3 # expt2-expt1, so km5c-pi is negative because pliocene
                # oceans were reduced in d18o.



##########################################################################

# timeslices are xiboi=preindustrial, xibol=3205 - km5c', xoekc=3205-km5c, xoekd=3060 (K1), xoeke=2950 (G17), xoekf=3155 (KM3)

#expt1='xjplc'
#extra1='6'
#startdec1='7'
#enddec1='9'
#timeslice1='KM5C'
#icevolcorr1=-0.3

#extra2='6'
#startdec2='7'
#enddec2='9'
#icevolcorr2=-0.3


#expt2='xjpld'
#timeslice2='K1'

#expt2='xjple'
#timeslice2='G17'

#expt2='xjplf'
#timeslice2='KM3'

expt1='xqbwa'
extra1='0000039'
startyear1=20
endyear1=50
timeslice1='pi'
icevolcorr1=0.0

expt2='xqetc'
extra2='0000059'
startyear2=70
endyear2=99
timeslice2='18ka'
icevolcorr2=0.0

 

################################################################
d18o_precip(expt1,extra1,startyear1,endyear1,timeslice1,expt2,extra2,startyear2,endyear2,timeslice2,icevolcorr1,icevolcorr2)

d18o_ocean(expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2)

salinity_ocean(expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2)

d18o_runoff(expt1,extra1,startdec1,enddec1,timeslice1,expt2,extra2,startdec2,enddec2,timeslice2,icevolcorr1,icevolcorr2)



sys.exit(0)

####

