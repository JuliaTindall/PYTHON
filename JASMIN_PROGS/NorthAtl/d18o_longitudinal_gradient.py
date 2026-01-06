#!/usr/bin/env python2.7
#NAME
#    D18o longitudinal gradient
#PURPOSE 
#    
#    Jonathan Holmes wants to see if the longitudinal gradient of d18o in the
#    model agrees with the data.
#    This will be done for the 9-11ka - 7-8ka period
#
#    He actually wants the longitudinal gradient just at 
#    the sites that he has
#
#
#
# Julia 30.09.20181


# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import sys
from netCDF4 import Dataset, MFDataset
from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid
from scipy import stats


# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,cbartitle,minval,maxval,diffval):
    lons, lats = np.meshgrid(lon,lat)
    if fileno <= 90:
        plt.subplot(2,2,fileno+1)
   # map=Basemap(projection='robin',resolution='l')

   # this may be a good region for the North Atlantic region
    map=Basemap(llcrnrlon=-20.0,urcrnrlon=40.0,llcrnrlat=30.0,urcrnrlat=70.0,projection='cyl',resolution='c')   

    # this may be a good region for most of the globe
    #map=Basemap(llcrnrlon=0,llcrnrlat=-80,urcrnrlon=360,urcrnrlat=80,projection='mill')

    
   # map.drawmapboundary(fill_color='aqua')
    x, y = map(lons, lats)
    map.drawcoastlines()
    #V=np.arange(np.amin(plotdata),np.amax(plotdata),np.amin(plotdata)/10)
    V=np.arange(minval,maxval,diffval)
   # cs = map.contourf(x,y,plotdata,V)
    cs=map.contourf(x,y,plotdata,V,extend='both',cmap='RdBu_r')
    #cs=map.contourf(x,y,plotdata,vmin=minval,vmax=maxval)
    plt.title(titlename)
  
    cbar=map.colorbar(cs,location='bottom',pad="5%")
    cbar.ax.tick_params(labelsize=10)
    cbar.ax.set_title(cbartitle)


def oplotdata(anom_lakes_d18o,datalons,datalats,minval,maxval,diffval,lake_spel):
# plot filled circles of the data over the map
    print(anom_lakes_d18o)
    V=np.arange(minval,maxval,diffval)
    print(minval,maxval)
    if lake_spel=='l': # lake
        plt.scatter(datalons,datalats,color='black',marker='o',s=110)
        plt.scatter(datalons,datalats,vmin=minval,vmax=maxval,c=anom_lakes_d18o,marker='o',s=70,cmap='RdBu_r')
    if lake_spel=='s': # speleothem
        plt.scatter(datalons,datalats,color='black',marker='^',s=110)
        plt.scatter(datalons,datalats,vmin=minval,vmax=maxval,c=anom_lakes_d18o,marker='^',s=70,cmap='RdBu_r')
    #plt.colorbar()
    


def getKey(item):
    return item[0]



#=====================================================
def extract_data(expt,filenames,expttime):


    # 1.  Set up details.  Gridbox required and filename


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

   # titlename=expttime
   # cbartitle='mm'
    lontemp=lon
    tot_h2o,lon=shiftgrid(180.,tot_h2o,lon,start=False,cyclic=360)
   # plotdata(tot_h2o*60.*60.*24.*30/tsize,0,lon,lat,titlename,cbartitle,0,150,10)
   # cbartitle='permil'
    lon=lontemp
    d18o_weightavg,lon=shiftgrid(180.,d18o_weightavg,lon,start=False,cyclic=360)
#    plotdata(d18o_weightavg,2,lon,lat,titlename,cbartitle,-40,10,5)
 #   plt.show()

    retdata=[d18o_weightavg,tot_h2o/tsize,lon,lat]

    return retdata





#=====================================================
def extract_data_seas(expt,filestart,expttime,monthnames):

    # 1.  Set up details.  Gridbox required and filename

    nmonths=len(monthnames)
    seasname=''  # get seasonname by using first letter of each month
    for mon in monthnames:
        seasname=seasname+mon[0]

    dirname='/home/users/jctindall/umoutput/BAS_timeslices/'+expt
    os.chdir(dirname)


    #2. Set up filename and extract stash code 338
    #   print d18o and dD for that file


    for monthno in range (0,nmonths):
        filenames=filestart+monthnames[monthno]+'.nc'
        print(dirname+filenames)
        f=MFDataset(filenames)
        f.dimensions
        f.variables

        lon = f.variables['longitude_1'][:]
        lat = f.variables['latitude_1'][:]
        times = f.variables['t'][:]
        xsize=len(lon)
        ysize=len(lat)
        tsize=len(times)
        print(filenames,tsize)

        alldata=f.variables['QCL'][:]
        h2o=alldata[:,0,:,:]+alldata[:,3,:,:]+alldata[:,6,:,:]+alldata[:,9,:,:]
        tot18o=alldata[:,1,:,:]+alldata[:,4,:,:]+alldata[:,7,:,:]+alldata[:,10,:,:]
    
        if monthno ==0: 
            # we need to get annual average by weighting by precipitation amount
            all_tot18o=tot18o
            all_h2o=h2o
            tsizesave=tsize
        else:
            if tsizesave != tsize:
                print('you have not got the same number of files for each month')
                print(tsizesave,tsize)

                sys.exit()

            all_tot18o=all_tot18o+tot18o
            all_h2o=all_h2o+h2o
    
        f.close()


    tot18o_avg=np.sum(all_tot18o,axis=0)  
    toth2o_avg=np.sum(all_h2o,axis=0)

    lontemp=lon
    toth2o_avg,lon=shiftgrid(180.,toth2o_avg,lon,start=False,cyclic=360)
    lon=lontemp
    tot18o_avg,lon=shiftgrid(180.,tot18o_avg,lon,start=False,cyclic=360)
    d18o_avg=((tot18o_avg/toth2o_avg)-2005.2E-6)/2005.2E-9


    retdata=[d18o_avg,toth2o_avg/tsize,lon,lat]

    return retdata

#====================================================================
def get_Jonathan_data(filename):
    
    f1=open(filename,'r')
    #discard titleline
    textline=f1.readline()

    datalons=[]
    datalats=[]
    data_5_7ka=[]
    data_9_11ka=[]

    for line in f1:
        linesplit=line.split(',') # the data in the file is split by comma
        datalats.append(np.float(linesplit[2]))
        datalons.append(np.float(linesplit[3]))
        data_5_7ka.append(np.float(linesplit[5]))
        data_9_11ka.append(np.float(linesplit[7]))

   
    retdata=[datalons,datalats,data_5_7ka,data_9_11ka]
    return retdata

#======================================================
def full_lon_plot(lon,lat,d18o_anom):
# this will plot the longitudinal gradient averaged over
# a range of latitudes just from the model

# make reduced array
    latmin=40.
    latmax=60.
    count=0
    for i in range(0,len(lon)):
        if lon[i] >= -10. and lon[i] <= 30.:
            count=count+1
    lonredu=np.zeros(count)
    count=0
    for j in range(0,len(lat)):
        if lat[j] >= latmin and lat[j] <= latmax:
            count=count+1
    latredu=np.zeros(count)

    d18oanom_redu=np.zeros((len(latredu),len(lonredu)))
    loncount=0
    for i in range(0,len(lon)):
        if lon[i] >= -10. and lon[i] <= 30.:
            lonredu[loncount]=lon[i]
            latcount=0
            for j in range(0,len(lat)):
                if lat[j] >= latmin and lat[j] <= latmax:
                    if loncount==0:
                        latredu[latcount]=lat[j]
                    d18oanom_redu[latcount,loncount]=d18o_anom[j,i]
                    latcount=latcount+1
            loncount=loncount+1

    meand18oanom=np.mean(d18oanom_redu,axis=0)
    maxd18oanom=np.max(d18oanom_redu,axis=0)
    mind18oanom=np.min(d18oanom_redu,axis=0)

    plt.subplot(211)
    plt.plot(lonredu,meand18oanom,label='meand18o')
    plt.xlabel('longitude')
    plt.ylabel('permille')
    plt.title('longitudinal average'+str(latmin)+'N-'+str(latmax)+'N')
    plt.legend()

    plt.subplot(212)
    plt.plot(lonredu,maxd18oanom,label='maximum d18o')
    plt.plot(lonredu,mind18oanom,label='minimum d18o')
    plt.xlabel('longitude')
    plt.ylabel('permille')
    plt.legend()

   
    fileout='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_average_d18o.png' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()


#end def full_lon_plot

#########################################################
def site_lon_plot(lon,lat,d18oanom,lakedatalons,lakedatalats,lakedata_anom,speldatalons,speldatalats,speldata_anom):

    # mark lake outliers (on Jonathans plot these have 
    # longitude > 15 and d18o > 1)

    nlakes=len(lakedatalons)
    outliers=np.zeros(nlakes)
    for i in range(0,nlakes):
        if lakedatalons[i] > 15. and lakedata_anom[i] > 1.0:
            outliers[i]=1

        if lakedata_anom[i] < -1.5:
            outliers[i]=1
            print('found outlier')

    # temporary change so there are no outliers
    print('j2',lakedatalons)
    outliers[:]=0

    # set up model data at lake locations
    d18omod_lakeloc=np.zeros(nlakes)
    d18omod_lake_outliers=np.zeros(nlakes)
    d18omod_lake_nooutliers=np.zeros(nlakes)
    lonmod_l=np.zeros(nlakes)
    latmod_l=np.zeros(nlakes)

    # get nearest longitude and latitude
   
    lonss=np.zeros(nlakes)
    latss=np.zeros(nlakes)
    for lake in range(0,nlakes):
        lonss[lake]=(np.abs(lon-lakedatalons[lake])).argmin()
  
    
    for lake in range(0,nlakes):
        latss[lake]=(np.abs(lat-lakedatalats[lake])).argmin()
  
    lonss=lonss.astype(int) # change from float to int
    latss=latss.astype(int) # change from float to int
   
     
    for lake in range(0,nlakes):
        lonmod_l[lake]=lon[lonss[lake]]
        latmod_l[lake]=lat[latss[lake]]
        d18omod_lakeloc[lake]=d18oanom[latss[lake],lonss[lake]]

        if outliers[lake]==1:
            d18omod_lake_outliers[lake]=d18omod_lakeloc[lake]
            d18omod_lake_nooutliers[lake]=np.nan
        else:
            d18omod_lake_outliers[lake]=np.nan
            d18omod_lake_nooutliers[lake]=d18omod_lakeloc[lake]

    
       

    # set up model data at spel locations
    nspels=len(speldatalons)
    d18omod_spelloc=np.zeros(nspels)
    lonmod_s=np.zeros(nspels)
    latmod_s=np.zeros(nspels)

    # get nearest longitude and latitude
   
    lonss=np.zeros(nspels)
    latss=np.zeros(nspels)
    for spel in range(0,nspels):
        lonss[spel]=(np.abs(lon-speldatalons[spel])).argmin()
  
    
    for spel in range(0,nspels):
        latss[spel]=(np.abs(lat-speldatalats[spel])).argmin()
  
    lonss=lonss.astype(int) # change from float to int
    latss=latss.astype(int) # change from float to int
   
     
    for spel in range(0,nspels):
        lonmod_s[spel]=lon[lonss[spel]]
        latmod_s[spel]=lat[latss[spel]]
        d18omod_spelloc[spel]=d18oanom[latss[spel],lonss[spel]]


    # linear regression for speleothem
    slope_s,intercept_s,r_value,p_value,std_err=stats.linregress(lonmod_s,d18omod_spelloc)
    # linear regression for lake
    finiteYmask=np.isfinite(d18omod_lake_nooutliers)
    print('j1',lonmod_l)
    print('j1',d18omod_lake_nooutliers)
    slope_l,intercept_l,r_value,p_value,std_err=stats.linregress(lonmod_l[finiteYmask],d18omod_lake_nooutliers[finiteYmask])
      


    #plt.plot(lonmod_l,d18omod_lakeloc,'b^',label='lake')
    plt.plot(lonmod_l,d18omod_lake_outliers,'g^',markersize=10,label='lake outliers')
    plt.plot(lonmod_l,d18omod_lake_nooutliers,'bo',label='lake')
    plt.plot(lonmod_l,intercept_l+slope_l * lonmod_l,'b')
   
    plt.plot(lonmod_s,d18omod_spelloc,'ro',label='speleothem')
    plt.plot(lonmod_s,intercept_s+slope_s * lonmod_s,'r')
    plt.legend()
    plt.xlabel('longitude')
    plt.ylabel('permille')
    plt.title('linear regression at sites of data')
    fileout='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_linregress_at_sites.png' 
    plt.savefig(fileout, bbox_inches='tight')  
    fileout='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_linregress_at_sites.eps' 
    plt.savefig(fileout, bbox_inches='tight')  



    # print out lat, lon and d18odata and d18omodel to a text file
    filewrite='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_data_model.txt'
    f2=open(filewrite,'w+')
    f2.write('speleothem\n')
    f2.write('lon   lat   data(d18oanom)  model(d18oanom)\n')
    for spel in range(0,nspels):
       f2.write(str(speldatalons[spel])+' '+str(speldatalats[spel])+' '+str(round(speldata_anom[spel],2))+' '+str(round(d18omod_spelloc[spel],2))+'\n')
 
    f2.write('lakes \n')
    f2.write('lon   lat   data(d18oanom)  model(d18oanom)\n')
    for lake in range(0,nlakes):
        f2.write(str(lakedatalons[lake])+' '+str(lakedatalats[lake])+' '+str(round(lakedata_anom[lake],2))+' '+str(round(d18omod_lakeloc[lake],2))+'\n')
 
    f2.close()
#end def site_lon_plot


#=================================================================
# MAIN PROGRAM STARTS HERE

seasname='jfmamjjasond'
#seasname='son'
if seasname =='djf':
    monthnames=['dc','ja','fb']
if seasname =='jja':
    monthnames=['jn','jl','ag']
if seasname =='mam':
    monthnames=['mr','ar','my']
if seasname =='son':
    monthnames=['sp','ot','nv']
if seasname =='jfmamjjasond':
    monthnames=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']

#retdata=extract_data('xlubh','iso_xlubha@pcr[3-6]*','7ka')


# xlubh 7ka
retdata=extract_data_seas('xlubh','xlubha@pcr[3-6]*','7ka',monthnames)
d18o_7ka=retdata[0]+0.04 # constant is ice volume correction
h2o_7ka=retdata[1]
lon=retdata[2]
lat=retdata[3]

# xlubg 6ka
retdata=extract_data_seas('xlubg','xlubga@pcr[3-6]*','6ka',monthnames)
d18o_6ka=retdata[0]+0.04
h2o_6ka=retdata[1]
lon=retdata[2]
lat=retdata[3]

# xlubf 5ka
retdata=extract_data_seas('xlubf','xlubfa@pcr[3-6]*','7ka',monthnames)
d18o_5ka=retdata[0]+0.04
h2o_5ka=retdata[1]
lon=retdata[2]
lat=retdata[3]




# xlubj 9ka
retdata=extract_data_seas('xlubj','xlubja@pcr[3-6]*','9ka',monthnames)
print('julia',retdata[0])
d18o_9ka=retdata[0]+0.25 # constant is ice volume correction
h2o_9ka=retdata[1]


# xlubk 10ka
retdata=extract_data_seas('xlubk','xlubka@pcr[3-6]*','10ka',monthnames)
d18o_10ka=retdata[0]+0.33
h2o_10ka=retdata[1]


# xlubl 11ka
retdata=extract_data_seas('xlubl','xlubla@pcr[3-6]*','11ka',monthnames)
d18o_11ka=retdata[0]+0.42
h2o_11ka=retdata[1]


# get average 7ka d18o - weighted by amount

d18o_5_7ka=((d18o_7ka * h2o_7ka) +  (d18o_6ka * h2o_6ka)+  (d18o_5ka * h2o_5ka))/(h2o_7ka+h2o_6ka+h2o_5ka)
d18o_9_11ka=((d18o_9ka * h2o_9ka) +  (d18o_10ka * h2o_10ka) + (d18o_11ka * h2o_11ka))/(h2o_9ka+h2o_10ka+h2o_11ka)

d18o_anom=d18o_9_11ka - d18o_5_7ka

print(np.shape(d18o_anom))


# do a full longitudinal plot at all model gridboxes
# averaged over a range  

full_lon_plot(lon,lat,d18o_anom)



# get Jonathans data
#filename='/home/users/jctindall/programs/PYTHON/NorthAtl/data/lakes_data.csv'
filename='/home/users/jctindall/programs/PYTHON/NorthAtl/data/Jonathan_new_data.csv'
lakedata=get_Jonathan_data(filename)
lakedatalons=lakedata[0]
lakedatalats=lakedata[1]
lakedata_5_7ka=lakedata[2]
lakedata_9_11ka=lakedata[3]

filename='/home/users/jctindall/programs/PYTHON/NorthAtl/data/speleothem_data.csv'
speleodata=get_Jonathan_data(filename)
speleodatalons=speleodata[0]
speleodatalats=speleodata[1]
speleodata_5_7ka=speleodata[2]
speleodata_9_11ka=speleodata[3]


lakedata_anom=np.asarray(lakedata_9_11ka) - np.asarray(lakedata_5_7ka)

speleodata_anom=np.asarray(speleodata_9_11ka) - np.asarray(speleodata_5_7ka)


# do a plot at the sites
site_lon_plot(lon,lat,d18o_anom,lakedatalons,lakedatalats,lakedata_anom,speleodatalons,speleodatalats,speleodata_anom)



h2o_5_7ka=(h2o_7ka+h2o_6ka+h2o_5ka)/3.0
h2o_9_11ka=(h2o_9ka+h2o_10ka+h2o_11ka)/3.0


plotdata(d18o_9_11ka,0,lon,lat,'9-11ka d18o','permille',-12,-6,0.5)
plotdata(h2o_9_11ka*60.*60.*24.*30.,1,lon,lat,'9-11ka h2o','mm/month',0,100,10)

anom_h2o=h2o_9_11ka - h2o_5_7ka
titlename='precip 9-11ka - 5-7ka'
cbartitle='mm/month'
plotdata(anom_h2o*60.*60.*24.*30,3,lon,lat,titlename,cbartitle,-10,10,1)

print(np.shape(d18o_5_7ka), np.shape(d18o_9_11ka))

anom_d18o=d18o_9_11ka - d18o_5_7ka

print(np.shape(lakedata_5_7ka), np.shape(lakedata_9_11ka))
anom_lakes_d18o=np.asarray(lakedata_9_11ka) - np.asarray(lakedata_5_7ka)
anom_speleo_d18o=np.asarray(speleodata_9_11ka) - np.asarray(speleodata_5_7ka)

titlename='d18o 9-11ka - 5-7ka '+seasname
cbartitle='permille'
#plotdata(anom_d18o,2,lon,lat,titlename,cbartitle,-6,2,0.5)
plotdata(anom_d18o,2,lon,lat,titlename,cbartitle,-2,1,0.2)


fileout='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_9-11ka_5-7ka_'+seasname+'.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_9-11ka_5-7ka_'+seasname+'.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

plotdata(anom_d18o,99,lon,lat,titlename,cbartitle,-2.0,2.0,0.25)
oplotdata(anom_lakes_d18o,lakedatalons,lakedatalats,-2.0,2.0,0.25,'l')
oplotdata(anom_speleo_d18o,speleodatalons,speleodatalats,-2.0,2.0,0.25,'s')

print('lake',lakedatalons)
print('speleo',speleodatalons)


fileout='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_9-11ka_5-7ka_d18o_data'+seasname+'.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/home/users/jctindall/plots/python/NorthAtl/plot_d18o/North_atl_9-11ka_5-7ka_d18o_data'+seasname+'.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()



sys.exit()
