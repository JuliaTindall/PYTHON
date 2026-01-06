#!/usr/bin/env python2.7
#NAME
#    Nizar_plot_diags
#PURPOSE 
#    PLOT Diagnostics from Louise/Max timeslice experiments for Nizar
#
#
# Julia 31.1.2017


# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import sys
from netCDF4 import Dataset, MFDataset
from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid



# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,cbartitle,minval,maxval,diffval):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(2,2,fileno+1)
   # map=Basemap(projection='robin',resolution='l')

   # this may be a good region for middle east
    map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')   

    # this may be a good region for most of the globe
   # map=Basemap(llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=80,projection='mill')

    
   # map.drawmapboundary(fill_color='aqua')
    x, y = map(lons, lats)
    map.drawcoastlines()
    #V=np.arange(np.amin(plotdata),np.amax(plotdata),np.amin(plotdata)/10)
    V=np.arange(minval,maxval,diffval)
   # cs = map.contourf(x,y,plotdata,V)
    cs=map.contourf(x,y,plotdata,V,extend='both')
    plt.title(titlename)
  
    cbar=map.colorbar(cs,location='bottom',pad="5%")
    cbar.ax.tick_params(labelsize=6)
    cbar.ax.set_title(cbartitle)

def getKey(item):
    return item[0]

def plottimeseries(field,timeperiod,fieldname,axisname):



#   plot my data and overplot Nizars data
    plt.subplot(3,1,1)
    print(np.shape(timeperiod),timeperiod)
    print(np.shape(field),field)
    plt.plot(timeperiod,field)
    plt.xlabel('ka')
    plt.ylabel(axisname)
    plt.title(fieldname+' over Jordan')


    plt.tight_layout(pad=0.4,w_pad=0.5,h_pad=1.0)




#=====================================================
def extract_single_site_data(expt,filenames,expttime,fieldreq,lonname,latname):

    if fieldreq=='convP':
        extract338='convP'
        fieldreq='QCL'

    if fieldreq=='lsP':
        extract338='lsP'
        fieldreq='QCL'
    

    # 1.  Set up details.  Gridbox required and filename

    #Jordan
    lonreq=37.5
    latreq=30.0

    dirname='/home/users/jctindall/umoutput/BAS_timeslices/'+expt+'/'
    os.chdir(dirname)


    #2. Set up filename and extract the field we want


    print(filenames)
    f=MFDataset(filenames)
    f.dimensions
    f.variables

    lon = f.variables[lonname][:]
    lat = f.variables[latname][:]
    times = f.variables['t'][:]
    xsize=len(lon)
    ysize=len(lat)
    tsize=len(times)
    fieldval=f.variables[fieldreq][:]
    fieldval=np.squeeze(fieldval)




    if fieldreq=='QCL':
        if extract338=='convP':
            fieldval2=fieldval
            fieldval=fieldval2[:,9,:,:]+fieldval2[:,6,:,:]
            fieldval=fieldval * 60. * 60. * 24. * 30.
            fieldval2=0
            print(np.shape(fieldval))
        if extract338=='lsP':
            fieldval2=fieldval
            fieldval=fieldval2[:,0,:,:]+fieldval2[:,3,:,:]
            fieldval=fieldval * 60. * 60. * 24. * 30.
            fieldval2=0
            print(np.shape(fieldval))

    nt,ny,nx=np.shape(fieldval)

#############temporary code for testing


# extract Jordan value

#    ix1=(lon == lonreq)
#    ix2=(lat ==latreq)
    
#    temp=fieldval[:,ix2]
#    single_pt_avg=temp[:,:,ix1]
    
    #print('single point average is',single_pt_avg,expt)
#    if fieldreq=='precip':
#        single_pt_avg=single_pt_avg * 60. * 60. * 24. * 30.
#    datasize=len(single_pt_avg)
#    print(datasize)
#    for i in range(0,datasize):
#        print(fieldreq,i,single_pt_avg[i])
#    #print('mm month',single_pt_avg)
#    print('mm year',np.sum(single_pt_avg))

#############end temporary code for testing



    # check that nt is a multiple of 120 (so that we get a full year)

    print('j1')
    checkval=nt%12  # remainder after dividing by 120
    if checkval != 0:
        print('error there are some files missingin ',expt)
        print('you will not be able to get an accurate annual average')
        sys.exit()

    
    f.close()

    # we need to get annual average by weighting by precipitation amount

    avgfield=np.mean(fieldval,axis=0)
    print('j2')
    


    # extract Jordan value

    ix1=(lon == lonreq)
    ix2=(lat ==latreq)

    temp=avgfield[ix2]
    single_pt_avg=temp[:,ix1]
    
    print('single point average is',single_pt_avg,expt)
    

    #cbartitle='mm'
    #titlename='avg'
    #plotdata(tot_h2o*60.*60.*24.*30/tsize,0,lon,lat,titlename,cbartitle,0,100,5)

    

    retdata=[single_pt_avg]

    return retdata


#=================================================================
# MAIN PROGRAM STARTS xluba



# HERE 0ka
retdata=extract_single_site_data('xluba','xlubaa@pcr[6-9]*.nc','0ka','convP','longitude_1','latitude_1')
alldata=[['xluba',0.0,retdata]]
print(retdata)

# xlubb 1ka
retdata=extract_single_site_data('xlubb','xlubba@pcr[3-6]*.nc','1ka','convP','longitude_1','latitude_1')
alldata.append(['xlubb',-1.0,retdata])
print(retdata)

# xlubc 2ka
retdata=extract_single_site_data('xlubc','xlubca@pcr[3-6]*.nc','2ka','convP','longitude_1','latitude_1')
alldata.append(['xlubc',-2.0,retdata])

# xlubd 3ka
retdata=extract_single_site_data('xlubd','xlubda@pcr[3-6]*.nc','3ka','convP','longitude_1','latitude_1')
alldata.append(['xlubd',-3.0,retdata])


# xlube 4ka
retdata=extract_single_site_data('xlube','xlubea@pcr[3-6]*.nc','4ka','convP','longitude_1','latitude_1')
alldata.append(['xlube',-4.0,retdata])

# xlubf 5ka
retdata=extract_single_site_data('xlubf','xlubfa@pcr[3-6]*.nc','5ka','convP','longitude_1','latitude_1')
alldata.append(['xlubf',-5.0,retdata])

# xlubg 6ka
retdata=extract_single_site_data('xlubg','xlubga@pcr[3-6]*.nc','6ka','convP','longitude_1','latitude_1')
alldata.append(['xlubg',-6.0,retdata])

# xlubh 7ka
retdata=extract_single_site_data('xlubh','xlubha@pcr[3-6]*.nc','7ka','convP','longitude_1','latitude_1')
alldata.append(['xlubh',-7.0,retdata])

# xlubi 8ka
retdata=extract_single_site_data('xlubi','xlubia@pcr[6-9]*.nc','8ka','convP','longitude_1','latitude_1')
alldata.append(['xlubi',-8.0,retdata])

# xlubj 9ka
retdata=extract_single_site_data('xlubj','xlubja@pcr[3-6]*.nc','9ka','convP','longitude_1','latitude_1')
alldata.append(['xlubj',-9.0,retdata])

# xlubk 10ka
retdata=extract_single_site_data('xlubk','xlubka@pcr[3-6]*.nc','10ka','convP','longitude_1','latitude_1')
alldata.append(['xlubk',-10.0,retdata])

# xlubl 11ka
retdata=extract_single_site_data('xlubl','xlubla@pcr[3-6]*.nc','11ka','convP','longitude_1','latitude_1')
alldata.append(['xlubl',-11.0,retdata])

# xlubm 12ka
retdata=extract_single_site_data('xlubm','xlubma@pcr[3-6]*.nc','12ka','convP','longitude_1','latitude_1')
alldata.append(['xlubm',-12.0,retdata])

# xlubn 13ka
retdata=extract_single_site_data('xlubn','xlubna@pcr[3-6]*.nc','13ka','convP','longitude_1','latitude_1')
alldata.append(['xlubn',-13.0,retdata])

# xlubo 14ka
retdata=extract_single_site_data('xlubo','xluboa@pcr[3-6]*.nc','14ka','convP','longitude_1','latitude_1')
alldata.append(['xlubo',-14.0,retdata])

# xlubp 15ka
retdata=extract_single_site_data('xlubp','xlubpa@pcr[0-2]*.nc','15ka','convP','longitude_1','latitude_1')
alldata.append(['xlubp',-15.0,retdata])

# xlubq 16ka
retdata=extract_single_site_data('xlubq','xlubqa@pcr[0-2]*.nc','16ka','convP','longitude_1','latitude_1')
alldata.append(['xlubq',-16.0,retdata])

# xlubr 17ka
retdata=extract_single_site_data('xlubr','xlubra@pcr[6-9]*.nc','17ka','convP','longitude_1','latitude_1')
alldata.append(['xlubr',-17.0,retdata])

# xlubs 18ka
retdata=extract_single_site_data('xlubs','xlubsa@pcr[2-5]*.nc','18ka','convP','longitude_1','latitude_1')
alldata.append(['xlubs',-18.0,retdata])

# xlubt 19ka
retdata=extract_single_site_data('xlubt','xlubta@pcr[2-5]*.nc','19ka','convP','longitude_1','latitude_1')
alldata.append(['xlubt',-19.0,retdata])

# xlubu 20ka
retdata=extract_single_site_data('xlubu','xlubua@pcr[2-5]*.nc','20ka','convP','longitude_1','latitude_1')
alldata.append(['xlubu',-20.0,retdata])

# xlubv 21ka
retdata=extract_single_site_data('xlubv','xlubva@pcr[2-5]*.nc','21ka','convP','longitude_1','latitude_1')
alldata.append(['xlubv',-21.0,retdata])



# extract d18o

datalen=len(alldata)
convP=[alldata[i][2] for i in range(0,datalen)]
timeperiod=[alldata[i][1] for i in range(0,datalen)]


#TdegC=np.asarray(lsP)-273.15
convP=np.squeeze(convP)
print('convP',convP)
print('timeperiod',timeperiod)
plottimeseries(convP,timeperiod,'total convP','mm/month')
fileout='/home/users/jctindall/plots/python/Middleeast/Nizar_plot_diags/Jordan_convP_timeseries.eps' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/home/users/jctindall/plots/python/Middleeast/Nizar_plot_diags/Jordan_convP_timeseries.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
