#!/usr/bin/env python2.7
#NAME
#    PLOT_D18O
#PURPOSE 
#    PLOT D18O FROM LOUISE/MAX TIMESLICE EXPERIMENTS MULTIPLE FILES
#
#
# Julia 20.09.2016


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
    map=Basemap(llcrnrlon=-90.0,urcrnrlon=0.0,llcrnrlat=-50.0,urcrnrlat=0.0,projection='cyl',resolution='c')   

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

def plottimeseries(d18o,dD,timeperiod):




#   plot my data and overplot Nizars data
    plt.subplot(3,1,1)
    plt.plot(timeperiod,d18o)
    plt.xlabel('ka')
    plt.ylabel('permil')
    plt.title('d18o_p over paraiso')

    plt.subplot(3,1,2)
    plt.plot(timeperiod,13.65+(6.53*np.asarray(d18o)))
    plt.xlabel('ka')
    plt.ylabel('permil')
    plt.title('13.65-(6.53*d18o) over paraiso')


    plt.subplot(3,1,3)
    plt.plot(timeperiod,dD)
    plt.xlabel('ka')
    plt.ylabel('permil')
    plt.title('dD over paraiso')

    #plt.subplot(4,1,4)
    #dxs=dD-(8.0*np.asarray(d18o))
    #plt.plot(timeperiod,dxs)
    #plt.xlabel('ka')
    #plt.ylabel('permil')
    #plt.title('dxs over paraiso')

    plt.tight_layout(pad=0.4,w_pad=0.5,h_pad=1.0)




#=====================================================
def extract_single_site_data(expt,filenames,expttime):


    # 1.  Set up details.  Gridbox required and filename

    #Paraiso
    lonreq=303.75
    latreq=-5.0

    #Tigre perdido
    #lonreq=281.25
    #latreq=-5.0

    dirname='/home/users/jctindall/mholloway/'+expt+'/pcpd/'
    os.chdir(dirname)


    #2. Set up filename and extract stash code 338
    #   print d18o and dD for that file


    f=MFDataset(filenames)
    f.dimensions
    f.variables

    lon = f.variables['longitude_1'][:]
    lat = f.variables['latitude_1'][:]
    times = f.variables['t'][:]
    xsize=len(lon)
    ysize=len(lat)
    tsize=len(times)
    dD=f.variables['dD'][:]
    dD=np.squeeze(dD)
    d18o=f.variables['dO18'][:]
    d18o=np.squeeze(d18o)
    h2o=f.variables['h2o'][:]
    h2o=np.squeeze(h2o)
    
    
    f.close()

    # we need to get annual average by weighting by precipitation amount

    dD_weight=dD*h2o
    dD_weightsum=np.sum(dD_weight,axis=0)
    d18o_weight=d18o*h2o
    d18o_weightsum=np.sum(d18o_weight,axis=0)
    tot_h2o=np.sum(h2o,axis=0)
    dD_weightavg=dD_weightsum/tot_h2o
    d18o_weightavg=d18o_weightsum/tot_h2o
    #print("shape tpt h2o",tot_h2o.shape)
    #print("shape mean dDweight",dD_weightsum.shape)


    # extract Paraiso d18o and dD and total annual precipitation

    ix1=(lon == lonreq)
    ix2=(lat ==latreq)
   
    temp=d18o_weightavg[ix2]
    single_pt_d18o=np.mean(temp[:,ix1])
    temp=dD_weightavg[ix2]
    single_pt_dD=np.mean(temp[:,ix1])
    temp=tot_h2o[ix2]
    single_pt_h2o=np.mean(temp[:,ix1])
    
    

    cbartitle='mm'
    titlename='avg'

    lontemp=lon
    d18o_weightavg,lon=shiftgrid(180,d18o_weightavg,lon,start=False)
    lon=lontemp
    tot_h2o,lon=shiftgrid(180,tot_h2o,lon,start=False)
    lon=lontemp
    dD_weightavg,lon=shiftgrid(180,dD_weightavg,lon,start=False)
   
    plotdata(tot_h2o*60.*60.*24.*30/tsize,0,lon,lat,titlename,cbartitle,0,200,10)

    

    #weighted dD
    titlename='dD'
    cbartitle='permil '
    plotdata(dD_weightavg,1,lon,lat,titlename,cbartitle,-100,200,20)

    #weighted d18o
    titlename='d18o'
    cbartitle='permil '
    plotdata(d18o_weightavg,2,lon,lat,titlename,cbartitle,-12,0,1)

    #weighted dxs
    titlename='dxs'
    cbartitle='permil '
    dxs=dD_weightavg - (8.0 * d18o_weightavg)
    plotdata(dxs,3,lon,lat,titlename,cbartitle,0,100,5)

    plt.tight_layout(pad=0.4,w_pad=0.5,h_pad=1.0)

    fileout='/home/users/jctindall/plots/python/Amazon/plot_d18o/'+expt+'_'+expttime+'_d18o_spatial.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    fileout='/home/users/jctindall/plots/python/Amazon/plot_d18o/'+expt+'_'+expttime+'_d18o_spatial.png' 
    plt.savefig(fileout, bbox_inches='tight')  

  
    plt.close()
   

    retdata=[single_pt_d18o,single_pt_dD,single_pt_h2o]

    return retdata


#=================================================================
# MAIN PROGRAM STARTS HERE


# xluba 0ka
retdata=extract_single_site_data('xluba','iso_xlubaa@pcr[3-6]*.nc','0ka')
alldata=[['xluba',0.0,retdata]]

# xlubb 1ka
retdata=extract_single_site_data('xlubb','iso_xlubba@pcr[3-6]*.nc','1ka')
alldata.append(['xlubb',-1.0,retdata])

# xlubc 2ka
retdata=extract_single_site_data('xlubc','iso_xlubca@pcr[3-6]*.nc','2ka')
alldata.append(['xlubc',-2.0,retdata])

# xlubd 3ka
retdata=extract_single_site_data('xlubd','iso_xlubda@pcr[3-6]*.nc','3ka')
alldata.append(['xlubd',-3.0,retdata])

# xlube 4ka
retdata=extract_single_site_data('xlube','iso_xlubea@pcr[3-6]*.nc','4ka')
alldata.append(['xlube',-4.0,retdata])

# xlubf 5ka
retdata=extract_single_site_data('xlubf','iso_xlubfa@pcr[3-6]*.nc','5ka')
alldata.append(['xlubf',-5.0,retdata])

# xlubg 6ka
retdata=extract_single_site_data('xlubg','iso_xlubga@pcr[3-6]*.nc','6ka')
alldata.append(['xlubg',-6.0,retdata])

# xlubh 7ka
retdata=extract_single_site_data('xlubh','iso_xlubha@pcr[3-6]*.nc','7ka')
alldata.append(['xlubh',-7.0,retdata])

# xlubi 8ka
retdata=extract_single_site_data('xlubi','iso_xlubia@pcr[3-6]*.nc','8ka')
alldata.append(['xlubi',-8.0,retdata])

# xlubj 9ka
retdata=extract_single_site_data('xlubj','iso_xlubja@pcr[3-6]*.nc','9ka')
alldata.append(['xlubj',-9.0,retdata])

# xlubk 10ka
retdata=extract_single_site_data('xlubk','iso_xlubka@pcr[3-6]*.nc','10ka')
alldata.append(['xlubk',-10.0,retdata])

# xlubl 11ka
retdata=extract_single_site_data('xlubl','iso_xlubla@pcr[3-6]*.nc','11ka')
alldata.append(['xlubl',-11.0,retdata])

# xlubm 12ka
retdata=extract_single_site_data('xlubm','iso_xlubma@pcr[3-6]*.nc','12ka')
alldata.append(['xlubm',-12.0,retdata])

# xlubn 13ka
retdata=extract_single_site_data('xlubn','iso_xlubna@pcr[3-6]*.nc','13ka')
alldata.append(['xlubn',-13.0,retdata])

# xlubo 14ka
retdata=extract_single_site_data('xlubo','iso_xluboa@pcr[3-6]*.nc','14ka')
alldata.append(['xlubo',-14.0,retdata])

# xlubp 15ka
retdata=extract_single_site_data('xlubp','iso_xlubpa@pcr[0-2]*.nc','15ka')
alldata.append(['xlubp',-15.0,retdata])

# xlubq 16ka
retdata=extract_single_site_data('xlubq','iso_xlubqa@pcr[0-2]*.nc','16ka')
alldata.append(['xlubq',-16.0,retdata])

# xlubr 17ka
retdata=extract_single_site_data('xlubr','iso_xlubra@pcr[0-2]*.nc','17ka')
alldata.append(['xlubr',-17.0,retdata])

# xlubs 18ka
retdata=extract_single_site_data('xlubs','iso_xlubsa@pcr[2-5]*.nc','18ka')
alldata.append(['xlubs',-18.0,retdata])

# xlubt 19ka
retdata=extract_single_site_data('xlubt','iso_xlubta@pcr[2-5]*.nc','19ka')
alldata.append(['xlubt',-19.0,retdata])

# xlubu 20ka
retdata=extract_single_site_data('xlubu','iso_xlubua@pcr[2-5]*.nc','20ka')
alldata.append(['xlubu',-20.0,retdata])

# xlubv 21ka
retdata=extract_single_site_data('xlubv','iso_xlubva@pcr[2-5]*.nc','21ka')
alldata.append(['xlubv',-21.0,retdata])



# extract d18o

datalen=len(alldata)
d18o=[alldata[i][2][0] for i in range(0,datalen)]
dD=[alldata[i][2][1] for i in range(0,datalen)]
timeperiod=[alldata[i][1] for i in range(0,datalen)]


plottimeseries(d18o,dD,timeperiod)
fileout='/home/users/jctindall/plots/python/Amazon/plot_d18o/Paraiso_d18o_timeseries.eps' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/home/users/jctindall/plots/python/Amazon/plot_d18o/Paraiso_d18o_timeseries.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
