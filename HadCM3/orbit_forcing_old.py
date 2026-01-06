#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
#NAME
#    orbit_forcing
#PURPOSE
#    This program will assess the orbital forcing between two different 
# experiments.  We are actually trying to replicate the difference in incoming
# shortwave radiation
#
#
# search for 'main program' to find end of functions
# Julia January 2019



import os
import numpy as np
import scipy as sp
#import cf
import iris
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import netCDF4
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
#import basemap
#from mpl_toolkits.basemap import Basemap, shiftgrid




def get_incom_sw(filestart1,filestart2,expt1,expt2,daily):

    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    
  

    for i in range(0,len(monthnames)):
        # read in data from expt1 files
        filename=filestart1+monthnames[i]+'.nc'
        f=Dataset(filename)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        days=f.variables['t'][:]
        atemp=f.variables['field200'][:] # get incoming sw radiation
        atemp=np.squeeze(atemp)
        f.close()

        if i==0:
            data1=np.zeros((len(lat),len(monthnames)*len(days)))
            data2=np.zeros((len(lat),len(monthnames)*len(days)))
            timesteps=len(monthnames)*len(days)

        daystart=i*len(days)
        dayend=(i+1)*len(days)

    
        # read in data from expt2 files
        f=Dataset(filestart2+monthnames[i]+'.nc')
        atemp2=f.variables['field200'][:] # get incoming sw radiation
        atemp2=np.squeeze(atemp2)
        f.close()

        if daily == 'y':
            data1[:,daystart:dayend]=np.swapaxes(atemp[0:len(days),:,0],0,1)
            data2[:,daystart:dayend]=np.swapaxes(atemp2[0:len(days),:,0],0,1)
        else:
            data1[:,i]=atemp[:,0]
            data2[:,i]=atemp2[:,0]

            
    
    # get maximum from first experiment
    max_data1_np=np.argmax(data1[0,:]) # find max at northpole
    max_data1_sp=np.argmax(data1[len(lat)-1,:]) # find max at southpole
    print('max expt1')
    print('np',max_data1_np,data1[0,max_data1_np])
    print('sp',max_data1_sp,data1[0,max_data1_sp])


    # get maximum from second experiment
    max_data2_np=np.argmax(data2[0,:]) # find max at northpole
    max_data2_sp=np.argmax(data2[len(lat)-1,:]) # find max at southpole
    print('max expt2')
    print('np',max_data2_np,data2[0,max_data2_np])
    print('sp',max_data2_sp,data2[0,max_data2_sp])

    nh_shift=max_data2_np-max_data1_np
    sh_shift=max_data2_sp-max_data1_sp

    # if the shift calculated from the NH and the SH are within two days 
    # use the mean

    if np.abs(nh_shift-sh_shift) > 300:
        if sh_shift < 0:
            sh_shift=sh_shift+360.
        if nh_shift < 0:
            nh_shift=nh_shift+360.

    if np.abs(nh_shift-sh_shift) <= 4:
        shiftreq=np.int(np.fix((nh_shift+sh_shift)/2.0))
    else:
        print('we are shifting differently in NH and SH')
        print('nh shift',max_data2_np-max_data1_np)
        print('sh shift',max_data2_sp-max_data1_sp)
        sys.exit(0)


    print('nh shift',nh_shift)
    print('sh shift',sh_shift)
    print('shift',shiftreq)

    data2_raw=data2
    data2_shifted=np.zeros(np.shape(data2))

    # shift data2 so that the maximum value in the northern hemisphere is in 
    # the same place as it is in data1

    shiftreq=shiftreq*(-1)
    
    for j in range(0,len(lat)):
        # shift values using np.roll
        data2_shifted[j,:]=np.roll(data2_raw[j,:],shiftreq,axis=0) 
        
    

    anom=data2_shifted-data1

    V=np.arange(0,600,50)
    
    plt.subplot(2,2,1)
    print('j1',np.shape(len(monthnames)))
    print('j2',np.shape(lat))
    print('j3',np.shape(data1))

    plt.contourf(np.arange(0,timesteps),lat,data1,V,extend='max')
    plt.title('PI')
   
    plt.colorbar()

    plt.subplot(2,2,2)
    plt.contourf(np.arange(0,timesteps),lat,data2_raw,V,extend='max')
    plt.title('new orbit '+expt2+' - raw')
    plt.colorbar()

    plt.subplot(2,2,3)
    plt.contourf(np.arange(0,timesteps),lat,data2_shifted,V,extend='max')
    plt.title('new orb (month+'+np.str(shiftreq)+')')

    plt.colorbar()

    plt.subplot(2,2,4)
    V=np.arange(-100,105,5)
    plt.contourf(np.arange(0,timesteps),lat,anom,V,cmap='RdBu_r',extend='both')
    plt.title('new orbit - PI')
    if daily == 'y':
        plt.plot([151,151],[-90,90])
        plt.plot([270,270],[-90,90])
    else:
        plt.plot([5,5],[-90,90])
        plt.plot([9,9],[-90,90])

    plt.colorbar()


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/orbit_forcing/'+expt2+'-'+expt1+daily+'.eps'
    plt.savefig(fileout,bbox_inches='tight')
    plt.close()

    
    return(shiftreq)
 

#end def get_incom_sw

#######################################
def shiftexpt(shiftreq,expt2,filestart):


# read in 30 years of data for each month
    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    startyear=41
    endyear=89
    for mon in range(0,len(monthnames)):
        print('setting up month',mon)
        daystart=mon*30
        dayend=((mon+1)*30)
        # testing with iris
        # note there is no obvious way of automatically averaging becaues
        # our time dimension is 30
        
        for year in range(startyear,endyear):
            if mon == 0: print('setting up year ',year)
            # set up for loop
            if year < 100:
                extra='7'
                yearuse=year
            else:
                extra='8'
                yearuse=year-100
            if yearuse < 10:
                yearchar='0'+np.str(yearuse)
            else:
                yearchar=np.str(yearuse)
             
            # get data into iris cubes
            # atm
            filename=(filestart+expt2+'/netcdf/'+
                expt2+'a@pa'+extra+yearchar+monthnames[mon]+'.nc')
            
            atmcube=iris.load(filename)
            ncubes=len(atmcube)
            
            #ocean we have two fields with the same name get the first one
            filename=(filestart+expt2+'/netcdf/'+
                     expt2+'o@pf'+extra+yearchar+monthnames[mon]+'.nc')
            
            ocnallcube=iris.load(filename,'OCN TOP-LEVEL TEMPERATURE          K')
            ocncube=ocnallcube[0]
           
           
           
           
            if (year == startyear) & (mon==0):
                # set up empty cube from ocean
                days,odepths,olats,olons=np.shape(ocncube)
                ndays=days*12
                dataocn=np.empty((ndays,odepths,olats,olons))
                dataocn
                # loop over all cubes and set up empty cubes from atmosphere
                for c in range(0,ncubes):
                    # get dimensions
                    days,ndepths,nlats,nlons=np.shape(atmcube[c])
                   
                    # get lats in case we want to plot
                    #for coord in atmcube[0].coords():
                    #    #get latitudes
                    #    if coord.name()=='latitude':
                    #        lats=coord.points
                
                    # this is a horrible way of programming
                    # doing it like this to save time
                    if c == 0:
                        dataarr0=np.zeros((ndays,ndepths,nlats,nlons)) 
                        # get latitude from array 1 for testing
                    if c == 1:
                        dataarr1=np.zeros((ndays,ndepths,nlats,nlons))
                    if c == 2:
                        dataarr2=np.zeros((ndays,ndepths,nlats,nlons))
                    if c == 3:
                        dataarr3=np.zeros((ndays,ndepths,nlats,nlons))
                    if c == 4:
                        dataarr4=np.zeros((ndays,ndepths,nlats,nlons))
                    if c == 5:
                        dataarr5=np.zeros((ndays,ndepths,nlats,nlons))
                    if c == 6:
                        dataarr6=np.zeros((ndays,ndepths,nlats,nlons))
                    if c == 7:
                        print('error you need to add more files to store cube data')
                        sys.exit(0)
            else:
                 # check names in new cube match that in previous cube
                 if ocncube.long_name !=cubeprevocn.long_name:
                     print('fields in ocean are not in consistent order')
                     print(ocncube.long_name,'does not equal',cubeprevocn[c].long_name)
                     sys.exit(0)
               
                 for c in range(0,ncubes):
                     if atmcube[c].long_name != cubeprev[c].long_name:
                         print('fields in file are not in consistent order')
                         print(atmcube[c].long_name,'does not equal',cubeprev[c].long_name)
                         sys.exit(0)
                        
            
            # put data into cube
            # add ocean and atm data
            if year ==startyear:
                  dataocn[daystart:dayend,:,:,:]=ocncube.data
                  dataarr0[daystart:dayend,:,:,:]=atmcube[0].data
                  dataarr1[daystart:dayend,:,:,:]=atmcube[1].data
                  dataarr2[daystart:dayend,:,:,:]=atmcube[2].data
                  dataarr3[daystart:dayend,:,:,:]=atmcube[3].data
                  dataarr4[daystart:dayend,:,:,:]=atmcube[4].data
                  dataarr5[daystart:dayend,:,:,:]=atmcube[5].data
                  count=1
                  
            else:
                 dataocn[daystart:dayend,:,:,:]=(
                            (dataocn[daystart:dayend,:,:,:]+ocncube.data)) 
            
                 dataarr0[daystart:dayend,:,:,:]=(
                            (dataarr0[daystart:dayend,:,:,:]+atmcube[0].data))
                 dataarr1[daystart:dayend,:,:,:]=(
                            (dataarr1[daystart:dayend,:,:,:]+atmcube[1].data))
                 dataarr2[daystart:dayend,:,:,:]=(
                            (dataarr2[daystart:dayend,:,:,:]+atmcube[2].data))
                 dataarr3[daystart:dayend,:,:,:]=(
                            (dataarr3[daystart:dayend,:,:,:]+atmcube[3].data))
                 dataarr4[daystart:dayend,:,:,:]=(
                            (dataarr4[daystart:dayend,:,:,:]+atmcube[4].data))
                 dataarr5[daystart:dayend,:,:,:]=(
                            (dataarr5[daystart:dayend,:,:,:]+atmcube[5].data))
                 count=count+1
                
            
            cubeprev=atmcube
            cubeprevocn=ocncube
            
    # data should all be loaded in now
    # shift it using np.roll
    dataocn=np.where(dataocn > 10000.,dataocn,dataocn/count)
    dataarr0=dataarr0/count 
    dataarr1=dataarr1/count 
    dataarr2=dataarr2/count 
    dataarr3=dataarr3/count 
    dataarr4=dataarr4/count  
    dataarr5=dataarr5/count 
    
    dataocn_sh=np.roll(dataocn,shiftreq,axis=0)
    for c in range(0,ncubes):
        if c ==0:
            dataarr0_sh=np.roll(dataarr0,shiftreq,axis=0)
        if c ==1:
            dataarr1_sh=np.roll(dataarr1,shiftreq,axis=0)
        if c ==2:
            dataarr2_sh=np.roll(dataarr2,shiftreq,axis=0)
        if c ==3:
            dataarr3_sh=np.roll(dataarr3,shiftreq,axis=0)
        if c ==4:
            dataarr4_sh=np.roll(dataarr4,shiftreq,axis=0)
        if c ==5:
            dataarr5_sh=np.roll(dataarr5,shiftreq,axis=0)
        if c ==6:
            dataarr6_sh=np.roll(dataarr6,shiftreq,axis=0)
        
        
        
        
    # looks like we have put data in correct place so try and split into
    # monthly chunks
    
    np_slice = atmcube[4].extract(iris.Constraint(latitude=50,longitude=0))
   

    for mon in range(0,len(monthnames)):
        daystart=mon*30
        dayend=((mon+1)*30)
        fileout=(filestart+expt2+'/netcdf/'+
                expt2+'a@pa_avg'+monthnames[mon]+'.nc')
        atmcube[0].data=dataarr0_sh[daystart:dayend,:,:,:]
        atmcube[1].data=dataarr1_sh[daystart:dayend,:,:,:]
        atmcube[2].data=dataarr2_sh[daystart:dayend,:,:,:]
        atmcube[3].data=dataarr3_sh[daystart:dayend,:,:,:]
        atmcube[4].data=dataarr4_sh[daystart:dayend,:,:,:]
        atmcube[5].data=dataarr5_sh[daystart:dayend,:,:,:]
        iris.save(atmcube, fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E20)
        ocncube.data=dataocn_sh[daystart:dayend,:,:,:]
      
        fileout=(filestart+expt2+'/netcdf/'+
                expt2+'o@pf_avg'+monthnames[mon]+'.nc')
        iris.save(ocncube,fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=2.0E+20)
        
    #print('j2',atmcube[4].data[:,0,0,0]))
    #plt.plot(np_slice.data)
    #np_slice_new = atmcube[4].extract(iris.Constraint(latitude=50,longitude=0))
    #plt.plot(np_slice_new.data)
    #print('j2',np_slice_new.data)
    #plt.show()
    #print('may have sucessfully overwritten cube')
    #sys.exit(0)
    #plt.subplot(2,2,1)
    #plt.plot(dataarr4[:,0,0,0])
    #plt.plot(dataarr4[:,0,72,0])
    #plt.plot(dataarr4[:,0,35,0])
    #plt.title(atmcube[4].long_name)
    
    #plt.subplot(2,2,2)
    #plt.plot(dataarr4_sh[:,0,0,0])
    #plt.plot(dataarr4_sh[:,0,72,0])
    #plt.plot(dataarr4_sh[:,0,35,0])
    
   
    #plt.subplot(2,2,3)
    #plt.contourf(np.arange(360),lats,np.ndarray.transpose(dataarr2[:,0,:,0]))
    #plt.title(atmcube[1].long_name)
    #plt.colorbar()
    
    #plt.subplot(2,2,4)
    #plt.contourf(np.arange(360),lats,np.ndarray.transpose(dataarr2_sh[:,0,:,0]))
    #plt.title(atmcube[1].long_name)
    #plt.colorbar()
    
    
    
    #plt.show()
    
    
    # shift data by required amounts
    
    
    
        

    

#enddef shiftexpt
    



################################
# main program


# timeslices are xiboi=preindustrial, xibol=3205 - km5c', xoekc=3205-km5c, xoekd=3060 (K1), xoeke=2950 (G17), xoekf=3155 (KM3)
#
# others are xogzb, xogzc, xogzd, xogze, xogzf


#############################
# from monthly data         #
#############################
#expt1='xogzs'
#filestart1='/nfs/hera1/earjcti/um/xogzs/netcdf/xogzsa@pdz39'

#expt2='xogzw'
# doesnt matter what year we use
#filestart2='/nfs/hera1/earjcti/um/xogzw/netcdf/xogzwa@pd729'
#filestart2='/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy99'
#filestart2='/nfs/hera1/earjcti/Xiaofang/xhgfk_netcdf/pdfiles/xhgfka@pdy99'

#shiftreq=get_incom_sw(filestart1,filestart2,expt1,expt2,'n')
#print('new shift required',shiftreq)



######################################################
# from daily data                                    #
######################################################
##############################################
# get the shift required from the daily data
expt1='xogzs'  # preindustrial = xogzs
filestart1='/nfs/hera1/earjcti/um/xogzs/netcdf/xogzsa@paz99'

expt2='xogzt'
filestart2='/nfs/hera1/earjcti/um/xogzt/netcdf/xogzta@pa760'
#filestart2='\\Users\\julia\\OneDrive\\DATA\\HadCM3_DATA\\'

shiftreq=get_incom_sw(filestart1,filestart2,expt1,expt2,'y')
print('old shift required',shiftreq)


########################################
#shift by appropriate amount
filestart2='/nfs/hera1/earjcti/um/'
shiftexpt(shiftreq,expt2,filestart2)

# check shift required is now zero with new files
filestart2='/nfs/hera1/earjcti/um/xogzt/netcdf/xogzta@pa_avg'
shiftreq=get_incom_sw(filestart1,filestart2,expt1,expt2,'y')
print('new shift required',shiftreq)
################################################################



#sys.exit(0)

####

