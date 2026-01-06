#NAME
#    Database_temperature_annual.py
#PURPOSE 
#
#  This program will create database like files for annual averages.
#  it is actually used for Dan Hills energy balance calculation
#
#  The program will also get 
#         ?????_Annual_Average_a@pd_TotalPrecipitationRate.nc
#         ?????_Annaul_Average_a@pd_TotalCloud.nc
#
#
# Julia 8.2.2017
# Julia 20.10.2018 ; included the ability to create database HadCM3 files
# Julia 17.01.2025 ; copied from the CEMAC directory and adapted for pliomip3

# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import sys
from netCDF4 import Dataset, MFDataset
#from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid




#=====================================================
def get_monthly_data(expt,startyear,endyear,HadCM3):

# this function will extract the data and write out to a file

    if HadCM3 == 'orig':
        century=np.floor(startyear/100.)
        choices = {10 : 'a', 11: 'b',  12: 'c',  13: 'd', 14: 'e', 
                       15: 'f', 16: 'g', 17: 'h', 18: 'i', 19: 'j',  
                       20: 'k', 21: 'l', 22: 'm', 23: 'n', 24: 'o',
                       25: 'p', 26: 'q', 27: 'r', 28: 's', 29: 't',
                       30: 'u', 31: 'v', 32: 'w', 33: 'x', 34: 'y', 35:'z'} 

        extra=choices.get(century, np.str(np.int(century))) 
          
        sep = '@'
        startyearuse=np.int(startyear - (century * 100))
        endyearuse = np.int(endyear - (century * 100))
    else:
        sep = '#'
        staryearuse=startyear
        endyearuse=endyear
        extra=''


    outdir='/nfs/hera1/earjcti/um/'+expt+'/database_averages/'+expt+'_Monthly_Average_winds_'
    infile='/nfs/hera1/earjcti/um/'+expt+'/pcpd/'+expt+'a'+sep+'pc'
 

    # loop over all years
    months=['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
  
    for year in range(startyearuse,endyearuse):
        for i,month in enumerate(months):
            if HadCM3 == 'y':
                filename=infile+str(year).zfill(9)+month+'+.nc'
            if HadCM3 == 'orig':
                filename= infile + extra + str(year) + month + '.nc'
            print(filename)
            f=MFDataset(filename)
            #print(f.variables)
       
            if year == startyearuse and i==0:
                latin = f.variables['latitude'][:]
                latsize=len(latin)
                lonin = f.variables['longitude'][:]
                lonsize=len(lonin)
                pressin = f.variables['p'][:]
                zsize=len(pressin)
                allvar_u=np.ma.zeros((endyear-startyear,12,zsize,latsize,lonsize))
                allvar_v=np.ma.zeros((endyear-startyear,12,zsize,latsize,lonsize))
           
            uvar=f.variables['u'][:]
            vvar=f.variables['v'][:]
            allvar_u[year-startyearuse,i,:,:,:]=np.squeeze(uvar)
            allvar_v[year-startyearuse,i,:,:,:]=np.squeeze(vvar)
           


    avgvar_u=np.ma.mean(allvar_u,axis=0)
    avgvar_v=np.ma.mean(allvar_v,axis=0)
    #avgvar_u=np.ma.where(avgvar_u.mask == 1.0, -99999., avgvar_u)
    #avgvar_v=np.ma.where(avgvar_v.mask == 1.0, -99999., avgvar_v)
  
    # write average variable out to a netcdf file

    # set up filename

    fout=(outdir+'_#pc_uv_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
 
    f2=Dataset(fout,mode='w',format='NETCDF3_CLASSIC')
    # create dimensions
    lon=f2.createDimension('longitude',lonsize)
    lat=f2.createDimension('latitude',latsize)
    level=f2.createDimension('p',zsize)
    time=f2.createDimension('time',12)
    # create variables
    lons=f2.createVariable('longitude',np.float32,('longitude',))
    lats=f2.createVariable('latitude',np.float32,('latitude',))
    levels=f2.createVariable('p',np.float32,('p',))
    times=f2.createVariable('time',np.float32,('time',))

    varfield_u=f2.createVariable('u',np.float32,
                                   ('time','p','latitude','longitude'))
    longname_u=u"U COMPNT OF WIND ON PRESSURE LEVELS"
    unitsname_u=u"m s-1"

    varfield_v=f2.createVariable('v',np.float32,
                                   ('time','p','latitude','longitude'))
    longname_v=u"V COMPNT OF WIND ON PRESSURE LEVELS"
    unitsname_v=u"m s-1"

   
    # create variable attributes
    lons.setncatts({'units':u"degrees_east"}) 
    lats.setncatts({'units':u"degrees_north"})
    levels.setncatts({'units':u"m"})
    times.setncatts({'units':u"days since 0000-01-01 00.00",\
                         'calendar':"360_day"})
    varfield_u.setncatts({'long_name': longname_u,\
                                'units':unitsname_u})
    varfield_v.setncatts({'long_name': longname_v,\
                                'units':unitsname_v})

    # assign data to variables
    lons[:]=lonin
    lats[:]=latin
    levels[:]=pressin
    times[:]=(np.arange(0,360,30))+15
    varfield_u[:,:,:,:]=avgvar_u
    varfield_v[:,:,:,:]=avgvar_v
    
        
    f2.close()
        



    return 


#=================================================================
# MAIN PROGRAM STARTS HERE

#expt='xkvje'
#startyear=2301
#nyears=100
#HadCM3='n'

#expt='xogzl'
#startyear=2980
#nyears=50
#HadCM3='y'

#expt='xqbwc'
#expts=['xqbwc','xqbwd','xqbwe','xqbwg','xqbwh','xqbwi','xqbwj','xqbwk','xqbwl',
#       'xqbwm','xqbwr','xqbws','xqbwt']
#expts=['xqbwh','xqbwi','xqbwj','xqbwk','xqbwl',
#       'xqbwm','xqbwr','xqbws','xqbwt']
#expts=['xqbwn']
#startyear=3900
#endyear=4000
#HadCM3='y'

HadCM3='orig' # y / orig /n
expt='xozzf'
startyear=2450
endyear=2500

get_monthly_data(expt,startyear,endyear,HadCM3)
