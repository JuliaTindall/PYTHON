#!/usr/bin/env python2.7
#NAME
#    PLOT_regional_precip_indyear
#PURPOSE
#    This program will plot the precipitation (annual and seasonal) and
#    the precipitation anomaly (annual and seasonal) for each individual 
#    year for a specified region of the globe
#
# search for 'main program' to find end of functions
# Julia 9/2/2017



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
#  def get_annmean_precip
#  def get_seasmean_precip

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,lonmin,lonmax,latmin,latmax,plotbox):
    lons, lats = np.meshgrid(lon,lat)
    if fileno != 99:
        plt.subplot(3,3,fileno+1)
    else:
        plt.subplot(2,1,1)

   # this is good for the globe
    map=Basemap(llcrnrlon=lonmin-10,urcrnrlon=lonmax+10,llcrnrlat=latmin-10,urcrnrlat=latmax+10,projection='cyl',resolution='c')
    map.drawmapboundary
    x, y = map(lons, lats)
    map.drawcoastlines()
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
#        cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu')
#            cbar = plt.colorbar(cs,orientation="horizontal",extend='both')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
            else:
                print(np.shape(plotdata))
                cs = map.contourf(x,y,plotdata,V,extend='max')
#                cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)

    if plotbox =='y':  # overplot a box showing all the lats and lons
        plt.plot([lonmin,lonmin],[latmin,latmax],'white')
        plt.plot([lonmax,lonmax],[latmin,latmax],'white')
        plt.plot([lonmin,lonmax],[latmin,latmin],'white')
        plt.plot([lonmin,lonmax],[latmax,latmax],'white')
        

    if fileno == 99:  # single plot so add colorbar
        cbar = plt.colorbar(cs,orientation="horizontal")
        cbar.set_label(cbarname,labelpad=-40)
#end def plotdata

##############################################
def get_annmean_precip(exptname,yearstart,yearend,lonmin,lonmax,latmin,latmax,regionname):


    # setup initial values before loop
    plotno=0
    fileform='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/precip_data/'+exptname

    # loop over each year and plot
    for year in range(yearstart,yearend):

        # get fileyear and extra value
        century=np.floor(year/100.)

        choices = {10: 'a', 11: 'b', 12: 'c', 13: 'd', 14: 'e', 
                   15: 'f', 16: 'g', 17: 'h', 18: 'i', 19: 'j', 
                   20: 'k', 21: 'l', 22: 'm', 23: 'n', 24: 'o', 
                   25: 'p', 26: 'q', 27: 'r', 28: 's', 29: 't', 
                   30: 'u', 31: 'v', 32: 'w', 33: 'x', 34: 'y', 
                   35: 'z'}

        extra=choices.get(century,99) # the second value is the default v
                                         # value for if it is not found in
                                         # the choices list
      
        yearuse=np.int(year-(century * 100))
        yearuse=str("%02d"%yearuse)
        
        fname=fileform+'a@pd'+extra+np.str(yearuse)+'*.nc'

        print(yearuse)
        print(fname)
        f=MFDataset(fname)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        aprecip=f.variables['precip_1'][:]
        aprecip=np.squeeze(aprecip)
        ntimes,ny,nx=np.shape(aprecip)
        
        if ntimes !=12:
            print('you have not got 12 months in year ',year)
            print('you have',ntimes)
            sys.exit()
        
    
       #average across the time dimension
        pi_precip_ann=np.mean(aprecip,axis=0)
    
        pi_precip_ann=pi_precip_ann * 60. * 60. * 24. * 30.
       
      

        # shift grid
        lonprecip=lon
        pi_precip_ann,lon = shiftgrid(180.,pi_precip_ann,lon,start=False)
    

        # set up array with all data for averaging
        if year == yearstart:
            pi_precip_ann_allyears=np.zeros((yearend-yearstart,ny,nx))
            areamean_precip=np.zeros(yearend-yearstart)

            
        pi_precip_ann_allyears[year-yearstart,:,:]=pi_precip_ann


        # mask out latitudes required and get average over region

        ix1=(lon >=lonmin) & (lon <=lonmax)
        ix2=(lat >=latmin) & (lat <=latmax)
        lats_reg=lat[ix2]
        lons_reg=lon[ix1]
    
        mask_precip=pi_precip_ann[ix2]
        mask_precip=mask_precip[:,ix1]
        
        areamean_precip[year-yearstart]=np.mean(mask_precip)

        # plot the data and highlight the region of the average
       
        titlename=np.str(year)
        plotbox='y'
        plotdata(pi_precip_ann,plotno,lon,lat,titlename,0,300,1.0,0.0,'y','mm/month',lonmin,lonmax,latmin,latmax,plotbox)
        
        plotno=(plotno+1)%9

        if plotno ==0 or year==yearend:
            plt.subplots_adjust(bottom=0.2)
            cax=plt.axes([0.1,0.1,0.8,0.045])
            plt.colorbar(cax=cax,orientation='horizontal')
            axistitle='precip for '+exptname+ 'mm/month'
            plt.title(axistitle)
            
            fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_regional_precip_indyear/'+exptname+'_'+np.str(year)+'_'+regionname+'.eps'
            plt.savefig(fileout, bbox_inches='tight')  

            plt.close()

    # now plot average
    titlename=exptname+'_allavg_'+regionname
    plotdata(np.mean(pi_precip_ann_allyears,axis=0),99,lon,lat,titlename,0,300,1.0,0.0,'y','mm/month',lonmin,lonmax,latmin,latmax,plotbox)
        
    # plot average over region
    plt.subplot(2,1,2)
    plt.plot(areamean_precip)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_regional_precip_indyear/'+exptname+'_annaverage_'+regionname+'.eps'
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()




#end def annmean


def seasmean(m1,m2,m3,seasname,exptname,yearstart,yearend,lonmin,lonmax,latmin,latmax,regionname):

    # set up stuff for plotting depending on region
    uselog='n'
    plotmax=500
    if regionname=='Sahara':
        uselog='y'
        plotmax=300

    # setup initial values before loop
    plotno=0
    fileform='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/precip_data/'+exptname

    # loop over each year and plot
    for year in range(yearstart,yearend):

        # get fileyear and extra value
        century=np.floor(year/100.)

        choices = {10: 'a', 11: 'b', 12: 'c', 13: 'd', 14: 'e', 
                   15: 'f', 16: 'g', 17: 'h', 18: 'i', 19: 'j', 
                   20: 'k', 21: 'l', 22: 'm', 23: 'n', 24: 'o', 
                   25: 'p', 26: 'q', 27: 'r', 28: 's', 29: 't', 
                   30: 'u', 31: 'v', 32: 'w', 33: 'x', 34: 'y', 
                   35: 'z'}

        extra=choices.get(century,99) # the second value is the default v
                                         # value for if it is not found in
                                         # the choices list
      
        yearuse=np.int(year-(century * 100))
        yearuse=str("%02d"%yearuse)
        
        fname=fileform+'a@pd'
        y1=np.str(yearuse)
        y2=np.str(yearuse)
        y3=np.str(yearuse)
        e1=extra
        e2=extra
        e3=extra


        # if 1 is november or december then we need to use previous year
        if m1 == 'nv' or m1 == 'dc':
            print('reformatting month1')
            century=np.floor((year-1)/100.)
            e1=choices.get(century,99) 
            yearuse=np.int((year-1)-(century * 100))
            y1=str("%02d"%yearuse)

        if m2 == 'dc':
            print('reformatting month2')
            century=np.floor(year-1/100.)
            e1=choices.get(century,99) 
            yearuse=np.int(year-1-(century * 100))
            y2=str("%02d"%yearuse)

        

        
        print('yearuse=',yearuse)
        print(fname+e1+np.str(y1)+m1+'_precip.nc')
        print(fname+e2+np.str(y2)+m2+'_precip.nc')
        print(fname+e3+np.str(y3)+m3+'_precip.nc')
        print(' ')
        fa=Dataset(fname+e1+np.str(y1)+m1+'_precip.nc')
        fb=Dataset(fname+e2+np.str(y2)+m2+'_precip.nc')
        fc=Dataset(fname+e3+np.str(y3)+m3+'_precip.nc')

        lat = fa.variables['latitude'][:]
        lon = fa.variables['longitude'][:]
        aprecip=fa.variables['precip_1'][:]
        bprecip=fb.variables['precip_1'][:]
        cprecip=fc.variables['precip_1'][:]
        aprecip=np.squeeze(aprecip)
        bprecip=np.squeeze(bprecip)
        cprecip=np.squeeze(cprecip)
        ny,nx=np.shape(aprecip)


        # average across the time dimension
        precip_seas=np.mean([aprecip,bprecip,cprecip],axis=0)
        precip_seas=precip_seas * 60. * 60. * 24. * 30.
       
      

        # shift grid
        lonprecip=lon
        precip_seas,lon = shiftgrid(180.,precip_seas,lon,start=False)
    

        # set up array with all data for averaging
        if year == yearstart:
            precip_seas_allyears=np.zeros((yearend-yearstart,ny,nx))
            areamean_precip=np.zeros(yearend-yearstart)

            
        precip_seas_allyears[year-yearstart,:,:]=precip_seas


        # mask out latitudes required and get average over region

        ix1=(lon >=lonmin) & (lon <=lonmax)
        ix2=(lat >=latmin) & (lat <=latmax)
        lats_reg=lat[ix2]
        lons_reg=lon[ix1]
    
        mask_precip=precip_seas[ix2]
        mask_precip=mask_precip[:,ix1]
        
        areamean_precip[year-yearstart]=np.mean(mask_precip)

        # plot the data and highlight the region of the average
       
        titlename=np.str(year)
        plotbox='y'
        plotdata(precip_seas,plotno,lon,lat,titlename,0,plotmax,1.0,0.0,uselog,'mm/month',lonmin,lonmax,latmin,latmax,plotbox)
        
        plotno=(plotno+1)%9

        if plotno ==0 or year==yearend:
            plt.subplots_adjust(bottom=0.2)
            cax=plt.axes([0.1,0.1,0.8,0.045])
            plt.colorbar(cax=cax,orientation='horizontal')
            axistitle='precip for '+exptname+ 'mm/month'
            plt.title(axistitle)
            
            fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_regional_precip_indyear/'+exptname+'_'+np.str(year)+'_'+seasname+'_'+regionname+'.eps'
            plt.savefig(fileout, bbox_inches='tight')  

            plt.close()

    # now plot average
    titlename=exptname+'_allavg_'+regionname
    plotdata(np.mean(precip_seas_allyears,axis=0),99,lon,lat,titlename,0,plotmax,1.0,0.0,uselog,'mm/month',lonmin,lonmax,latmin,latmax,plotbox)
        
    # plot average over region
    plt.subplot(2,1,2)
    plt.plot(areamean_precip)

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_regional_precip_indyear/'+exptname+'_'+seasname+'_'+regionname+'.eps'
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()





#end def seasmean

################################
# main program

# annual mean
figureno=0

#plt.figure(figureno)
#get_annmean_precip('xkvjf',2301,2400,-10.0,30.0,15.0,30.0,'Sahara')
#figureno=figureno+1

#djf mean
#plt.figure(figureno)
#seasmean('dc','ja','fb','djf','xkvje',2301,2400,15.0,30.0,-30.0,-5.0,'SAfrica')
seasmean('jl','ag','sp','jas','xkvje',2301,2400,75.0,85.0,10.0,25.0,'India')
seasmean('jl','ag','sp','jas','xkvjg',2301,2400,75.0,85.0,10.0,25.0,'India')
#figureno=figureno+1



#jja mean
#plt.figure(figureno)
#seasmean('jn','jl','ag',figureno,'jja')
#figureno=figureno+1


sys.exit(0)

####

