#NAME
#    omni_deglac_temp.py
#PURPOSE 
#    Lauren wants to put the temperature from the BBC runs on the omniglobe
#    This program will produce suitable png files for this on Annie
#
#
# Julia 8.2.2017


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
def plotdata(plotdata,lon,lat,minval,maxval,diffval):


    # set image size for omniglobe
#    h_inches = 30.0  # 3000 by 1500 i.e. OmniGlobe at 100dpi
#    v_inches = 15.0
#    plt.figure(figsize=(h_inches,v_inches))
  
    # make canvas as big as figure area (for OmniGlobe PNGs)
#    plt.subplots_adjust(0,0,1,1)



    # Cater for lack of 360.0 information (periodic)
    nlats = len(lat)
    nlons = len(lon)
    lon_p = np.zeros(nlons+1)
    plotdata_p=np.zeros((nlats,nlons+1))
    lon_p[:-1] = lon
    lon_p[nlons] = 360.
    plotdata_p[:,:-1]=plotdata
    plotdata_p[:,nlons]=plotdata[:,0]


    lons, lats = np.meshgrid(lon_p,lat)
    map=Basemap(resolution='c', llcrnrlon=0.0,llcrnrlat=-90, urcrnrlon=360.0,urcrnrlat=90)

    x, y = map(lons, lats)
    #map.drawcoaslines()
    map.drawcoastlines(linewidth=0.25)
    V=np.arange(minval,maxval,diffval)
    cs=map.contourf(x,y,plotdata_p,V,extend='both')
    
  



#=====================================================
def extract_SAT(expt):
# this function will extract the SAT from the file


    filename='/nfs/annie/earljg/data/BBC_runs/'+expt+'a.pdcljja.nc'

    print(filename)
    f=Dataset(filename)
    f.dimensions
    f.variables

    lon = f.variables['longitude'][:]
    lat = f.variables['latitude'][:]
    xsize=len(lon)
    ysize=len(lat)
    fieldval=f.variables['temp_mm_1_5m'][:]
    fieldval=np.squeeze(fieldval)


    retdata=[lon,lat,fieldval]

    return retdata


#=================================================================
# MAIN PROGRAM STARTS HERE

minyear=-120
maxyear=0.5
exptinfo=[]
exptinfo.append(["tdwz9",-120.,"120k"])
exptinfo.append(["tdwz8",-116.,"116k"])
exptinfo.append(["tdwz7",-112.,"112k"])
exptinfo.append(["tdwz6",-108.,"108k"])
exptinfo.append(["tdwz5",-104.,"104k"])
exptinfo.append(["tdwz4",-100.,"100k"])
exptinfo.append(["tdwz3",-96.,"96k"])
exptinfo.append(["tdwz2",-92.,"92k"])
exptinfo.append(["tdwz1",-88.,"88k"])
exptinfo.append(["tdwz0",-84.,"84k"])
exptinfo.append(["tdwzZ",-80.,"80k"])
exptinfo.append(["tdwzY",-78.,"78k"])
exptinfo.append(["tdwzX",-76.,"76"])
exptinfo.append(["tdwzW",-74.,"74"])
exptinfo.append(["tdwzV",-72.,"72k"])
exptinfo.append(["tdwzU",-70.,"70k"])
exptinfo.append(["tdwzT",-68.,"68k"])
exptinfo.append(["tdwzS",-66.,"66k"])
exptinfo.append(["tdwzR",-64.,"64k"])
exptinfo.append(["tdwzQ",-62.,"62k"])
exptinfo.append(["tdwzP",-60.,"60k"])
exptinfo.append(["tdwzO",-58,"58k"])
exptinfo.append(["tdwzN",-56.,"56k"])
exptinfo.append(["tdwzM",-54.,"54k"])
exptinfo.append(["tdwzL",-52.,"52k"])
exptinfo.append(["tdwzK",-50.,"50k"])
exptinfo.append(["tdwzJ",-48.,"48k"])
exptinfo.append(["tdwzI",-46.,"46k"])
exptinfo.append(["tdwzH",-44.,"44k"])
exptinfo.append(["tdwzG",-42.,"42k"])
exptinfo.append(["tdwzF",-40.,"40k"])
exptinfo.append(["tdwzE",-38.,"38k"])
exptinfo.append(["tdwzD",-36.,"36k"])
exptinfo.append(["tdwzC",-34.,"34k"])
exptinfo.append(["tdwzB",-32.,"32k"])
exptinfo.append(["tdwzA",-30.,"30k"])
exptinfo.append(["tdwzz",-28.,"28k"])
exptinfo.append(["tdwzy",-26.,"26k"])
exptinfo.append(["tdwzx",-24.,"24k"])
exptinfo.append(["tdwzw",-22.,"22k"])
exptinfo.append(["tdwzv",-21.,"21k"])
exptinfo.append(["tdwzu",-20.,"20k"])
exptinfo.append(["tdwzt",-19.,"19k"])
exptinfo.append(["tdwzs",-18.,"18k"])
exptinfo.append(["tdwzr",-17.,"17k"])
exptinfo.append(["tdwzq",-16.,"16k"])
exptinfo.append(["tdwzp",-15.,"15k"])
exptinfo.append(["tdwzo",-14.,"14k"])
exptinfo.append(["tdwzn",-13.,"13k"])
exptinfo.append(["tdwzm",-12.,"12k"])
exptinfo.append(["tdwzl",-11.,"11k"])
exptinfo.append(["tdwzk",-10.,"10k"])
exptinfo.append(["tdwzj",-9.,"9k"])
exptinfo.append(["tdwzi",-8.,"8k"])
exptinfo.append(["tdwzh",-7.,"7k"])
exptinfo.append(["tdwzg",-6.,"6k"])
exptinfo.append(["tdwzf",-5.,"5k"])
exptinfo.append(["tdwze",-4.,"4k"])
exptinfo.append(["tdwzd",-3.,"3k"])
exptinfo.append(["tdwzc",-2.,"2k"])
exptinfo.append(["tdwzb",-1.,"1k"])
exptinfo.append(["tdwza",-0.,"0k"])


nexpts=len(exptinfo)
print(nexpts)

# get the data from the files
for n in range(0,nexpts):
    indinfo=exptinfo[n]
    exptname=indinfo[0]
    timeperiod=indinfo[1]
    timetitle=indinfo[2]
    retdata=extract_SAT(exptname)
    lon=retdata[0]
    lat=retdata[1]
    SAT=retdata[2]
    
    
    if n==0: 
        alltime=np.zeros(nexpts)
        allSAT_nonint=np.zeros((nexpts,len(lat),len(lon)))
        print('shape',np.shape(allSAT_nonint))
        
    allSAT_nonint[n,:,:]=SAT
    alltime[n]=timeperiod


# interpolate onto a 500 year array

years=np.arange(minyear,maxyear,0.5)
print(years)
nreq=len(years)
print(years,nreq)

SAT_int=np.zeros((nreq,len(lat),len(lon))) # to store 500 years of data
for treq in range(0,nreq):
    for tfile in range(0,len(alltime)):

        # if there is a file for the required time use that
        if years[treq] == alltime[tfile]:
            SAT_int[treq,:,:]=allSAT_nonint[tfile,:,:]
            
        else:
        # if not interpolate
            if years[treq] > alltime[tfile] and  \
                    years[treq] < alltime[tfile+1]:
                diffupp=alltime[tfile+1]-years[treq]
                difflow=years[treq]-alltime[tfile]
                difffull=alltime[tfile+1]-alltime[tfile]
                SAT_int[treq,:,:]=(((difffull-difflow) * allSAT_nonint[tfile,:,:]) +  ((difffull-diffupp) * allSAT_nonint[tfile+1,:,:]))/difffull
                
                
              
# plot interpolated SAT to a file

SAT_int=SAT_int-273.15 # convert to celcius
for treq in range(0,nreq):

    plotdata(SAT_int[treq,:,:],lon,lat,-40,40,1.0)
    fileyear=120+years[treq]
    if fileyear < 10:
        fileyear='00'+str(fileyear)
    else:
        if fileyear < 100:
            fileyear='0'+str(fileyear)
        else:
            fileyear=str(fileyear)
    
    fname_out = '/nfs/annie/earjcti/CEMAC/omniplots/deglac_SAT/SAT_120k-'+fileyear+'K.eps'
    print(fname_out)  
    plt.savefig(fname_out,bbox_inches="tight",pad_inches=0.0)
    plt.close()

    # check by plotting difference in consecutive files
    
#    plotdata(SAT_int[treq+1,:,:]-SAT_int[treq,:,:],lon,lat,-2,2,0.1)
#    fname_out = '/nfs/annie/earjcti/CEMAC/omniplots/deglac_SAT/SAT_diff'+np.str#(years[treq+1])+'_'+np.str(years[treq])+'.eps'
#    print(fname_out)  
#    plt.savefig(fname_out,bbox_inches="tight",pad_inches=0.0)
#    plt.close()

    

sys.exit()
