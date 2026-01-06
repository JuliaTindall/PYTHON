#!/usr/bin/env python2.7
#NAME
#    PLOT_NINOindex
#PURPOSE
#    This program will plot
#
# search for 'main program' to find end of functions
# Julia 22/11/2016



import os
import numpy as np
import scipy as sp
import scipy.signal as sig
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid
import glob

#functions are:
#  def plotdata   plot data (on a x-y grid)
#  def indexplot  plot nino34 index on a line graph
#  def fullprint  print full array for debugging
#  get_NINO34_temperatures gets the NINO34 temperatures and writes to a file
#  plot_NINO34_temperatures   plots NINO34 temperatures from above file

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(3,2,fileno+1)

   # this is good for a ELNINO region
    map=Basemap(llcrnrlon=100.0,urcrnrlon=300.0,llcrnrlat=-20.0,urcrnrlat=20.0,projection='cyl',resolution='c')
   # this is good for the globe
   # map=Basemap(llcrnrlon=0.0,urcrnrlon=360.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
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
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='ar':
                    cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    cs = map.contourf(x,y,plotdata,V,extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
#end def plotdata

def indexplot(toplot,fileno,nyears,nmonths,data_sm,elninoarr,
              laninaarr,xmin,xmax,expt,control,meanval):
    plt.subplot(4,1,fileno+1)

    plt.xlim([xmin,xmax])
    plt.ylim([-2.0,2.0])
    xdata=np.arange(nyears*nmonths)
    xdata=xdata / float(nmonths)
    # plot data
    plt.plot(xdata,toplot)
    if fileno==0:
        titlename='ONI index '+expt+' based on '+control+' mean is '+np.str(meanval)
        plt.title(titlename)
    # overplot smoothed data
    plt.plot(xdata,data_sm,'-')
    # overplot zero line and +-0.5deg line
    plt.plot(xdata,np.zeros(nyears*nmonths))
    plt.plot(xdata,np.zeros(nyears*nmonths)+0.5)
    plt.plot(xdata,np.zeros(nyears*nmonths)-0.5)
    bar_width=1.0/12.0
    plt.bar(xdata,elninoarr,bar_width,color='red',edgecolor="none")
    plt.bar(xdata,laninaarr,bar_width,color='blue',edgecolor="none")
   
# 


   
# end def indexplot

def fullprint(printarr):
  from pprint import pprint
  opt = np.get_printoptions()
  np.set_printoptions(threshold='nan')
  pprint(printarr)
  np.set_printoptions(**opt)
#end def fullprint


def get_NINO34_temperatures(exptname):
# this subroutine will get the NINO3.4 temperature and write them to a file

    print('get nino34_temperatures')

    dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/temp_data/'
    os.chdir(dirname)
    #full_file_paths =  glob.glob(dirname+"*o@pf*")
    full_file_paths =  glob.glob("*o@pf*")
    nfiles=len(full_file_paths)
    nyears=np.ceil(nfiles/12.)
    nmonths=12
    print(nfiles,nyears,nmonths)


    # prepare output files

    fileout='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/NINODATA/'+exptname+'_NINO3.4_temp.txt'
    f1=open(fileout,'w+')
    f1.write('nyears='+str(nyears)+' nmonths='+str(nmonths)+'\n')
    f1.write('filename              year   month   NINO3.4 temperature\n')
        

    count=0
    for fname in full_file_paths:
        # extract year and month from file

        # read in data from file
        f=Dataset(fname,mode='r')
        if count == 0:
            lat = f.variables['latitude'][:]
            latsize=np.shape(lat)
            lon = f.variables['longitude'][:]
            lonsize=np.shape(lon)
        
        count==count+1
        temperature=f.variables['temp'][:] 
       
        temperature=np.squeeze(temperature)
# mask out all data that is not in nino3.4 region (120W-170W, 5N-5S) 

        ix1=(lon >=190) & (lon <=240)
        ix2=(lat <= 5) & (lat >=-5) 
        lats_reg = lat[ix2]
        lons_reg = lon[ix1]
        
        ninoint=temperature[ix2]
     
        ninotemp=ninoint[:,ix1]

        # to check that correct region is found
        #plotdata(ninotemp,0,lons_reg,lats_reg,'test',20,30,1,0,'n','tempdeg')
        #plt.show()


        # calculate average temperature in region
        weights=np.cos(np.radians(lats_reg))
        nino34_zonT=np.average(ninotemp,axis=0,weights=weights)
        nino34_T=np.average(nino34_zonT)
        #print(count,fname,nino34_T)

        # from filename obtain year and month
        year=int(float(fname[10:12]))
        month=fname[12:14]
        extra=fname[9:10]
        
        choices = {'a': 10, 'b': 11, 'c': 12, 'd': 13, 'e': 14, 
                   'f': 15, 'g': 16, 'h': 17, 'i': 18, 'j': 19, 
                   'k': 20, 'l': 21, 'm': 22, 'n': 23, 'o': 24, 
                   'p': 25, 'q': 26, 'r': 27, 's': 28, 't': 29, 
                   'u': 30, 'v': 31, 'w': 32, 'x': 33, 'y': 34, 
                   'z': 35}

        century=choices.get(extra,extra) # the second extra is the default v
                                         # value for if it is not found in
                                         # the choices list
       

        choices = {'ja': 0, 'fb': 1, 'mr': 2, 'ar': 3, 'my': 4, 
                   'jn': 5, 'jl': 6, 'ag': 7, 'sp': 8, 'ot': 9, 
                   'nv': 10, 'dc': 11}

        monthno=choices.get(month,-99) # the second extra is the default v
                                         # value for if it is not found in
                                         # the choices list
       
        year=(century * 100) + year
        print(fname,year,monthno,nyears,nmonths)
        f1.write(fname+';'+str(year)+';'+str(monthno)+';'+str(nino34_T)+'\n')

    
        
    f1.close()



#===NEW PROG=====================================
def plot_NINO34_temperatures(exptname,controlname,tempanomaly):
# this program will plot the temperatures obtained in get_NINO34_temperatures
# it will look like an el nino index
# if will also see if the index is different if we are using the reference 
# temperatures from a control simulation

    fig=plt.figure()

    # get data from experiment
    filein='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/NINODATA/'+exptname+'_NINO3.4_temp.txt'
    f=open(filein,'r')
    # read file to get number of years and months of data
    firstline=f.readline()
    b=firstline.split() # split into nyears and nmonths title

    c=b[0] # get nyears title
    b1=c.split('=')  # split by removing text before equals sign
    nyears=int(float(b1[1])) # this is nyears

    c=b[1] # get nmonths title
    b1=c.split('=')  # split by removing text before equals sign
    nmonths=int(float(b1[1])) # this is nmonths
    print('nyears',nyears)
    print('nmonths',nmonths)


    # discard second titleline
    secondline=f.readline()


    # read over all the rest of the data in the file, find the startyear and
    # the endyear

    yearmin=1E10
    yearmax=-1E10
    for line in f:
        # extract year month and nino index
        linesplit=line.split(';')
        year=int(float(linesplit[1]))
        if year < yearmin:
            yearmin=year
        if year > yearmax:
            yearmax=year
        
    f.close()

    # make 2 arrays to store data
    
    nyears=yearmax-yearmin+1
    allyears=np.empty(nyears)
    allyears[:]=np.NAN
    nino34index=np.empty((nyears,nmonths))
    nino34index[:,:]=np.NAN
    filelist=np.empty((nyears,nmonths),dtype=object)



    # open the file again, get the NINO3.4 index and add to the arrays
    print('yearminmax',yearmax,yearmin)

    f=open(filein,'r')
    # discard two title lines
    firstline=f.readline()
    firstline=f.readline()
    for line in f:
        # extract year month and nino index
        linesplit=line.split(';')
        year=int(float(linesplit[1]))
        month=int(float(linesplit[2]))
        ninoval=float(linesplit[3])

        allyears[year-yearmin]=year
        nino34index[year-yearmin,month]=ninoval
        filelist[year-yearmin,month]=linesplit[0]

    f.close()


    #####################################################
    # get data from control if appropriate

    if controlname != exptname:
        filein='/nfs/hera1/earjcti/um/HadGEM_data/'+controlname+'/NINODATA/'+controlname+'_NINO3.4_temp.txt'
        f=open(filein,'r')
        # read file to get number of years and months of data
        firstline=f.readline()
        b=firstline.split() # split into nyears and nmonths title

        c=b[0] # get nyears title
        b1=c.split('=')  # split by removing text before equals sign
        nyears_ct=int(float(b1[1])) # this is nyears

        c=b[1] # get nmonths title
        b1=c.split('=')  # split by removing text before equals sign
        nmonths_ct=int(float(b1[1])) # this is nmonths
        print('nyears_ct',nyears_ct)
        print('nmonths__ct',nmonths_ct)


        # discard second titleline
        secondline=f.readline()


        # read over all the rest of the data in the file, find the startyear and
        # the endyear

        yearmin_ct=1E10
        yearmax_ct=-1E10
        for line in f:
            # extract year month and nino index
            linesplit=line.split(';')
            year_ct=int(float(linesplit[1]))
            if year_ct < yearmin_ct:
                yearmin_ct=year_ct
            if year_ct > yearmax_ct:
                yearmax_ct=year_ct
           
                    
        f.close()

        # make 2 arrays to store data

        nyears_ct=yearmax_ct-yearmin_ct+1
        allyears_ct=np.empty(nyears_ct)
        allyears_ct[:]=np.NAN
        nino34index_ct=np.empty((nyears,nmonths))
        nino34index_ct[:,:]=np.NAN
        filelist_ct=np.empty((nyears,nmonths),dtype=object)
        

        
        # open the file again, get the NINO3.4 index and add to the arrays
        
        f=open(filein,'r')
        # discard two title lines
        firstline=f.readline()
        firstline=f.readline()
        for line in f:
            # extract year month and nino index
            linesplit=line.split(';')
            year_ct=int(float(linesplit[1]))
            month_ct=int(float(linesplit[2]))
            ninoval=float(linesplit[3])
            
            allyears_ct[year_ct-yearmin_ct]=year_ct
            nino34index_ct[year_ct-yearmin_ct,month_ct]=ninoval
            filelist_ct[year_ct-yearmin_ct,month_ct]=linesplit[0]
            
        f.close()


    

    # get the average annual cycle and remove
    if exptname == controlname:
        meananncyc=np.nanmean(nino34index,axis=0)
    else:
        # reduce temperature anomaly to account for the 
        # fact that the Pliocene is a warmer climate
        print(tempanomaly)
        meananncyc=np.nanmean(nino34index_ct,axis=0)-tempanomaly

    nino34index=nino34index-meananncyc

   
    toplot=np.reshape(nino34index,nyears*nmonths)
    allfiles=np.reshape(filelist,nyears*nmonths)
    smoothednino34index=[np.mean(toplot[i-1:i+1]) for i in range (0,nyears*nmonths)]


    
    # setup an el nino array and a la nina array
    # condition must be met for 5 months

    elninoarr=np.zeros(nyears*nmonths)
    laninaarr=np.zeros(nyears*nmonths)
    for i in range(1,nyears*nmonths-4):
        if elninoarr[i-1] != 0:    # if previous month is el nino
            if smoothednino34index[i] >= 0.5:
                elninoarr[i]=smoothednino34index[i]
            else:
                elninoarr[i]=0
        else:
            flag='y'
            for i2 in range(i,i+5):
                if smoothednino34index[i] <= 0.5:
                    flag='n'
            if flag == 'y':
                elninoarr[i]=smoothednino34index[i]
            else:
                elninoarr[i]=0
        
        if laninaarr[i-1] != 0:    # if previous month is la nino
            if smoothednino34index[i] <= -0.5:
                laninaarr[i]=smoothednino34index[i]
            else:
                laninaarr[i]=0
        else:
            flag='y'
            for i2 in range(i,i+5):
                if smoothednino34index[i] >= -0.5:
                    flag='n'
            if flag == 'y':
                laninaarr[i]=smoothednino34index[i]
            else:
                laninaarr[i]=0
                
            
    meanval=np.nanmean(smoothednino34index)
    indexplot(toplot,0,nyears,nmonths,smoothednino34index,elninoarr,
              laninaarr,0,25,exptname,controlname,meanval)
    indexplot(toplot,1,nyears,nmonths,smoothednino34index,elninoarr,
              laninaarr,25,50,exptname,controlname,meanval)
    indexplot(toplot,2,nyears,nmonths,smoothednino34index,elninoarr,
              laninaarr,50,75,exptname,controlname,meanval)
    indexplot(toplot,3,nyears,nmonths,smoothednino34index,elninoarr,
              laninaarr,75,100,exptname,controlname,meanval)

                

    if exptname != controlname:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_NINOindex/'+exptname+'_'+controlname+'_NINO3.4_ONI.eps'
    else:
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_NINOindex/'+exptname+'_NINO3.4_ONI.eps'

        fileout='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/NINODATA/'+exptname+'_NINO3.4_index.txt'
        f1=open(fileout,'w+')
        f1.write('filename  year month  smoothed NINO3.4 temperature    el nino ind   la nina ind\n')

        print(allfiles)

        for i in range(0,nyears*nmonths):
            if allfiles[i]:
                fname=allfiles[i]
                print(fname)
                print(fname[10:12])
                year=int(float(fname[10:12]))
                month=fname[12:14]
                extra=fname[9:10]
                
                choices = {'a': 10, 'b': 11, 'c': 12, 'd': 13, 'e': 14, 
                       'f': 15, 'g': 16, 'h': 17, 'i': 18, 'j': 19, 
                       'k': 20, 'l': 21, 'm': 22, 'n': 23, 'o': 24, 
                       'p': 25, 'q': 26, 'r': 27, 's': 28, 't': 29, 
                       'u': 30, 'v': 31, 'w': 32, 'x': 33, 'y': 34, 
                       'z': 35}
            
                century=choices.get(extra,extra) # the second extra is the default v
                # value for if it is not found in
                # the choices list
            
            
                choices = {'ja': 0, 'fb': 1, 'mr': 2, 'ar': 3, 'my': 4, 
                       'jn': 5, 'jl': 6, 'ag': 7, 'sp': 8, 'ot': 9, 
                       'nv': 10, 'dc': 11}
            
                monthno=choices.get(month,-99) # the second extra is the default v
                # value for if it is not found in
                # the choices list
                
                print(century,year)
                year=(century * 100) + year
                print(fname,year,monthno,nyears,nmonths)
                f1.write(fname+';'+str(year)+';'+str(monthno)+';'+str(smoothednino34index[i])+';'+str(elninoarr[i])+';'+str(laninaarr[i])+'\n')


    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()




#end def



#===NEW PROG=====================================
def plot_NINO34_patterns(exptname,fieldreq):
# this program will plot the patterns associated with an el nino


    filein='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/NINODATA/'+exptname+'_NINO3.4_index.txt'
    f=open(filein,'r')
    # read firstline and discard
    firstline=f.readline()

    arraysetup='n'
    for line in f:
        # extract filename
        linesplit=line.split(';')
        filename=linesplit[0]
        month=int(float(linesplit[2]))
        normval=float(linesplit[3])
        enind=float(linesplit[4])
        lnind=float(linesplit[5])

        # read in data from file
        if fieldreq == 'temp':
            ncdffile='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/temp_data/'+filename
        if fieldreq == 'precip_1':
            ncdffile='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/precip_data/'+filename[0:5]+'a@pd'+filename[9:15]+'precip.nc'


        f=Dataset(ncdffile,mode='r')
        f.dimensions
        f.variables

        filetemp = f.variables[fieldreq][:]
        filetemp=np.squeeze(filetemp)
        
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]

        f.close()

        xsize=len(lon)
        ysize=len(lat)
    
        # setup arrays nmonths nlats nlons
        if arraysetup=='n':
            enarr=np.zeros((12,ysize,xsize))
            lnarr=np.zeros((12,ysize,xsize))
            normarr=np.zeros((12,ysize,xsize))
            cten=np.zeros(12,dtype=np.int)
            ctln=np.zeros(12,dtype=np.int)
            ctnm=np.zeros(12,dtype=np.int)
            arraysetup='y'


        if enind != 0.0:
            enarr[month,:,:]=enarr[month,:,:]+filetemp
            cten[month]=cten[month] + 1
        elif lnind != 0.0:
            lnarr[month,:,:]=lnarr[month,:,:]+filetemp
            ctln[month]=ctln[month]+1
        else:
            normarr[month,:,:]=normarr[month,:,:]+filetemp
            ctnm[month]=ctnm[month]+1


   
            
    enarr=enarr/cten[:,np.newaxis, np.newaxis]
    lnarr=lnarr/ctln[:,np.newaxis, np.newaxis]
    normarr=normarr/ctnm[:,np.newaxis, np.newaxis]
        
   
       
    enmean=np.average(enarr,axis=0)
    lnmean=np.average(lnarr,axis=0)
    nmmean=np.average(normarr,axis=0)
    units='degC'
    if fieldreq == 'precip_1': 
        units='mm/month'
        enmean=enmean * 60. *60. * 24. * 30.
        lnmean=lnmean * 60. *60. * 24. * 30.
        nmmean=nmmean * 60. *60. * 24. * 30.


    if fieldreq == 'temp':
        # figures for temperature
        titlename=exptname+'elnino'
        plotdata(enmean,0,lon,lat,titlename,20.,30,1,0,'n',units)
        plotdata(lnmean,1,lon,lat,'lanina',20.,30,1,0,'n',units)
        enanom=enmean-nmmean
        plotdata(enanom,2,lon,lat,'elninoanom',-2.,2.1,0.1,0,'a',units)
        lnanom=lnmean-nmmean
        plotdata(lnanom,3,lon,lat,'lanina anom',-2.,2.1,0.1,0,'a',units)

    if fieldreq =='precip_1':
        # plots for precip
        # figures for temperature
        titlename=exptname+'elnino'
        plotdata(enmean,0,lon,lat,titlename,0.,300,10,0,'n',units)
        plotdata(lnmean,1,lon,lat,'lanina',0.,300,10,0,'n',units)
        enanom=enmean-nmmean
        plotdata(enanom,2,lon,lat,'elninoanom',-100.,110,10,0,'ar',units)
        lnanom=lnmean-nmmean
        plotdata(lnanom,3,lon,lat,'lanina anom',-100.,110,10.0,0,'ar',units)



    # get average temperature across pacific
    ix1=(lon >=100) & (lon <=300)
    ix2=(lat <= 20) & (lat >=-20) 
    lats_reg = lat[ix2]
    lons_reg = lon[ix1]
        
    avg_pac=nmmean[ix2]
    avg_pac=avg_pac[:,ix1]
        
    print(np.shape(avg_pac))
    for i in range(0,len(lons_reg)):
        for j in range(0,len(lats_reg)):
            if avg_pac[j,i] == 0:  # land point
                avg_pac[j,i]=float('NaN')
    
    WPminEP=np.nanmean(avg_pac[:,0:len(lons_reg)/2.0])- \
        np.nanmean(avg_pac[:,len(lons_reg)/2.0:len(lons_reg)])

    print('avg pac',avg_pac)
    
    print('mean is',np.nanmean(avg_pac),WPminEP)

    if fieldreq =='temp':
        titlename='normal W-E is:'+np.str(WPminEP)
        plotdata(nmmean-np.nanmean(avg_pac),4,lon,lat,titlename,-4.,4,0.1,0,'a','degC')





    # print out the average anomaly associated with elnino

    avg_enanom=enanom[ix2]
    avg_enanom=avg_enanom[:,ix1]
    avg_lnanom=lnanom[ix2]
    avg_lnanom=avg_lnanom[:,ix1]
   
    for i in range(0,len(lons_reg)):
        for j in range(0,len(lats_reg)):
            if avg_enanom[j,i] == 0 and avg_lnanom[j,i]==0:  # land point
                 avg_enanom[j,i]=float('NaN')
                 avg_lnanom[j,i]=float('NaN')
                 
    print('mean elnino signature',np.nanmean(avg_enanom))
    print('mean lanina signature',np.nanmean(avg_lnanom))
    print('mean absolute value elnino signature',np.nanmean(np.absolute(avg_enanom)))
    print('mean absolute value lanina signature',np.nanmean(np.absolute(avg_lnanom)))





    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_NINOindex/'+exptname+'_NINO3.4_patterns_'+fieldreq+'.eps'

    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()


    retdata=[enanom,lnanom,lon,lat]
    return retdata


#end def



#===NEW PROG=====================================
def plot_NINO34_freqspectra(exptnames):
# this program will plot the freqquency spectra of the NINO3.4 temperature


    nexpts=len(exptnames)

    allninotemp=[]
    allmonth=[]
    monthtemp=np.zeros(12)
    monthcount=np.zeros(12)
    titlename='power spectral density'



    for exptname in exptnames:

        filein='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/NINODATA/'+exptname+'_NINO3.4_temp.txt'
        f=open(filein,'r')
    # read first two lines and discard
        firstline=f.readline()
        firstline=f.readline()
        for line in f:
            # extract filename
            linesplit=line.split(';')
            filename=linesplit[0]
            month=int(float(linesplit[2]))
            nino34temp=float(linesplit[3])
            allninotemp.append(float(nino34temp))
            # stuff for removing average annual cycle
            allmonth.append(month)
            monthtemp[month]=monthtemp[month]+nino34temp
            monthcount[month]=monthcount[month]+1

        f.close()
        titlename=exptname+' '+titlename


    monthtemp=monthtemp/monthcount

 

    # remove annual cycle
    allninotemp=allninotemp-monthtemp[allmonth]

    titlename=titlename+' variance is ','%.2f' %np.var(allninotemp)


    # do a spectral analysis on the index
    Pxx_f, Pxx_den=sig.periodogram(allninotemp,1.0)
   
    Pxx_f=Pxx_f * 12. #  convert to cycles per year from cycles per month

    plt.subplot(2,1,1)
    nmonths=len(allninotemp)
    xarray=(np.arange(0,nmonths,1))/12.
    print('shapes')
    print(np.shape(xarray))
    print(np.shape(allninotemp))
    plt.plot(xarray,allninotemp) # raw data with restricted
    plt.title('NINO3.4 temp anomaly')
    plt.xlabel('year')
    plt.ylabel('nino3.4 temp')
    axes=plt.gca()
    axes.set_xlim([0,nmonths/12])
   


    plt.subplot(2,1,2)
    plt.plot(1.0/Pxx_f,Pxx_den) # plotting with period on y axis
    plt.title(titlename)
    axes=plt.gca()
    axes.set_xlim([0.5,150])
    axes.set_xscale('log')
    axes.set_xticks([1,2,5,10,20,40,80,160])
    axes.get_xaxis().set_major_formatter(mp.ticker.ScalarFormatter())
    plt.xlabel('period (years)')



    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_NINOindex/'+exptname+'_NINO3.4_freqspectra.eps'

    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()


    
   




#end def



################################
# main program

figureno=0
###############################
#get_NINO34_temperatures('xkvje')
#get_NINO34_temperatures('xkvjf')
#get_NINO34_temperatures('xkvjg')

##############################################
# get patterns associated with el nino for a chosen field
tempdata=plot_NINO34_patterns('xkvje','temp')
temp_enanom_pi=tempdata[0]
temp_lnanom_pi=tempdata[1]
lon=tempdata[2]
lat=tempdata[3]

tempdata=plot_NINO34_patterns('xkvjf','temp')
temp_enanom_plio=tempdata[0]
temp_lnanom_plio=tempdata[1]

tempdata=plot_NINO34_patterns('xkvjg','temp')
temp_enanom_pliop2=tempdata[0]
temp_lnanom_pliop2=tempdata[1]


# plot absolute change in temperature
plotdata(temp_enanom_pliop2-temp_enanom_pi,0,lon,lat,'diff in elnino temp (pliop2 and pi)',-1.,1.1,0.1,0,'ar','degC')
plotdata(temp_enanom_plio-temp_enanom_pi,1,lon,lat,'diff in elnino temp (plio and pi)',-1.,1.1,0.1,0,'ar','degC')
plotdata(temp_lnanom_pliop2-temp_lnanom_pi,2,lon,lat,'diff in lanina temp (pliop2 and pi)',-1.,1.1,0.1,0,'ar','degC')
plotdata(temp_lnanom_pliop2-temp_lnanom_pi,3,lon,lat,'diff in lanina temp (plio and pi)',-1.,1.1,0.1,0,'ar','degC')

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_NINOindex/xkvjg-xkvje_temp_patterns.eps'

plt.savefig(fileout, bbox_inches='tight')  
plt.close()



precipdata=plot_NINO34_patterns('xkvje','precip_1')
# extract fields for anomalies
precip_enanom_pi=precipdata[0]
precip_lnanom_pi=precipdata[1]
lon=precipdata[2]
lat=precipdata[3]

precipdata=plot_NINO34_patterns('xkvjf','precip_1')
precip_enanom_plio=precipdata[0]
precip_lnanom_plio=precipdata[1]

precipdata=plot_NINO34_patterns('xkvjg','precip_1')
precip_enanom_pliop2=precipdata[0]
precip_lnanom_pliop2=precipdata[1]

# plot absolute change in precipitation
plotdata(precip_enanom_pliop2-precip_enanom_pi,0,lon,lat,'diff in elnino precip (pliop2 and pi)',-100.,110,10.0,0,'ar','mm/month')
plotdata(precip_enanom_plio-precip_enanom_pi,1,lon,lat,'diff in elnino precip (plio and pi)',-100.,110,10.0,0,'ar','mm/month')
plotdata(precip_lnanom_pliop2-precip_lnanom_pi,2,lon,lat,'diff in lanina precip (pliop2 and pi)',-100.,110,10.0,0,'ar','mm/month')
plotdata(precip_lnanom_pliop2-precip_lnanom_pi,3,lon,lat,'diff in lanina precip (plio and pi)',-100.,110,10.0,0,'ar','mm/month')

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_NINOindex/xkvjg-xkvje_precip_patterns.eps'

plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# we can't do percentage change because the negative numbers make the plot really confusing






#########################################
#plot_NINO34_temperatures('xkvje')
# the temperature anomaly needs to be obtained from plot_NINO34_patterns

#tempanomaly=avg_pac_temp_PI-avg_pac_temp_plio # for normalising control index
#print(tempanomaly)
#plot_NINO34_temperatures('xkvjf','xkvje',tempanomaly)
#tempanomaly=avg_pac_temp_PI-avg_pac_temp_pliop2 # for normalising control index
#plot_NINO34_temperatures('xkvjg','xkvje',tempanomaly)


###########################
# plot freqency spectra to see what freqency of NINO
#plot_NINO34_freqspectra(['xkvja','xkvje'])
#plot_NINO34_freqspectra(['xkvjb','xkvjf'])
#plot_NINO34_freqspectra(['xkvjc','xkvjg'])


sys.exit(0)

####

