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


def get_glob_avg(varreq, latin, field):
    """
    gets the global average of the fields
    """
    zonalmean = np.mean(varreq, axis=(0,1,3))
    latrad = latin * 2. * np.pi / 360.
    coslatrad = np.cos(latrad)
    
    globavg = np.sum(zonalmean * coslatrad) / np.sum(coslatrad)
    if field == 'temp' or field =='temp_1':
        globavg = globavg-273.15
    if field == 'precip_1':
        globavg = globavg * 60. * 60. *24.
   
    return globavg


#=====================================================
def get_annual_data(expt,startyear,endyear,field):

# this function will extract the data and write out to a file

    outdir='/Uolstore/Research/a/hera1/earjcti/um/'+expt+'/database_averages/'+expt+'_Annual_Average'
    if field == 'SST':
        infile='/Uolstore/Research/a/hera1/earjcti/um/'+expt+'/pf/'+expt+'o#pf'
    elif field == 'SSS' or field == 'salinity':
        infile='/Uolstore/Research/a/hera1/earjcti/um/'+expt+'/pg/'+expt+'o#pg'
    else:
        infile='/Uolstore/Research/a/hera1/earjcti/um/'+expt+'/pcpd/'+expt+'a#pd'
    txtfile='/Uolstore/Research/a/hera1/earjcti/um/'+expt+'/database_averages/'+expt+'_'+ field + str(startyear) + '_' + str(endyear)  +'_global_avg.txt'
    ftxt=open(txtfile,'w')

    # loop over all years
    for year in range(startyear,endyear):
        filename=infile+str(year).zfill(9)+'*.nc'
        print('j1',filename,year)
        f=MFDataset(filename)
       
        if year == startyear:
            latin = f.variables['latitude'][:]
            latsize=len(latin)
            lonin = f.variables['longitude'][:]
            lonsize=len(lonin)
            if field == 'salinity':
                depthin = f.variables['depth_1'][:]
        
            if field == 'salinity':
                allvar=np.ma.zeros((endyear-startyear,20,latsize,lonsize))
            else:
                allvar=np.ma.zeros((endyear-startyear,latsize,lonsize))
            timeseries=np.zeros(endyear-startyear)


        if field == 'SST':
            print(filename)
            varreq=f.variables['temp'][:]
            allvar[year-startyear,:,:]=np.ma.mean(varreq,axis=0)
        elif field == 'SSS':
            print(filename)
            varreq=f.variables['salinity'][:]
            salin=np.ma.mean(varreq,axis=0)
            allvar[year-startyear,:,:]=(salin[0,:,:] * 1000.)+35.0
        elif field == 'salinity':
            print(filename)
            varreq=f.variables['salinity'][:]
            salin=np.ma.mean(varreq,axis=0)
            allvar[year-startyear,:,:,:]=(salin * 1000.)+35.0
        else:
            varreq=f.variables[field][:]
            allvar[year-startyear,:,:]=np.mean(varreq,axis=0)
      
        globavg = get_glob_avg(varreq, latin,field)
        
        timeseries[year-startyear]=globavg
        ftxt.write(np.str(year) + ',' + np.str(globavg) + '\n')
   
        f.close()
    ftxt.close()

    print('averaging now', np.shape(allvar))
    avgvar=np.ma.mean(allvar,axis=0)
    avgvar=np.ma.where(avgvar.mask == 1.0, -99999., avgvar)
  
    # write average variable out to a netcdf file

    # set up filename

    fout=(outdir+'_#pd_'+field+'_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
    if ((field == 'temp_1' and HadCM3 !='y')
         or (field == 'temp' and HadCM3 =='y')):
        fout=(outdir+'_#pd_Temperature_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
    if field == 'precip_1':
        fout=(outdir+'_#pd_TotalPrecipitationRate_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
    if field == 'field30':
        fout=(outdir+'_#pd_TotalCloud_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
    if field == 'SST':
        fout=(outdir+'_#pf_SST_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
    if field == 'SSS':
        fout=(outdir+'_#pf_SSS_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
    if field == 'salinity':
        fout=(outdir+'_#pg_salinity_' + str(startyear) + 
          '_' + str(endyear) + '.nc')
    print(fout)
 
    f2=Dataset(fout,mode='w',format='NETCDF3_CLASSIC')
    # create dimensions
    lon=f2.createDimension('longitude',lonsize)
    lat=f2.createDimension('latitude',latsize)
    level=f2.createDimension('ht',1)
    time=f2.createDimension('time',1)
    # create variables
    lons=f2.createVariable('longitude',np.float32,('longitude',))
    lats=f2.createVariable('latitude',np.float32,('latitude',))
    levels=f2.createVariable('ht',np.float32,('ht',))
    times=f2.createVariable('time',np.float32,('time',))
    if ((field == 'temp_1' and HadCM3 !='y')
         or (field == 'temp' and HadCM3 =='y')):
        varfield=f2.createVariable('temp',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"TEMPERATURE AT 1.5M"
        unitsname=u"K"
    if field == 'precip_1':
        varfield=f2.createVariable('precip',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"TOTAL PRECIPITATION RATE KG/M2/S"
        unitsname=u"kg m-2"
    if field == 'field30':
        varfield=f2.createVariable('field30',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"TOTAL CLOUD AMOUNT - RANDOM OVERLAP"
        unitsname=u"0-1"
    if field == 'field201':
        varfield=f2.createVariable('field201',np.float32,
                              ('time','ht','latitude','longitude'))
        longname=u"OUTGOING SW RAD FLUX (TOA)"
        unitsname=u"W m-2"
    if field == 'field207_1':
        varfield=f2.createVariable('field207_1',np.float32,
                              ('time','ht','latitude','longitude'))
        longname=u"CLEAR-SKY (II) UP SURFACE SW FLUX"
        unitsname=u"W m-2"
    if field == 'field200':
        varfield=f2.createVariable('field200',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"INCOMING SW RAD FLUX (TOA): ALL TSS"
        unitsname=u"W m-2"
    if field == 'field207':
        varfield=f2.createVariable('field207',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"CLEAR-SKY (II) UPWARD SW FLUX (TOA)"
        unitsname=u"W m-2"
    if field == 'ilr':
        varfield=f2.createVariable('ilr',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"DOWNWARD LW RAD FLUX: SURFACE"
        unitsname=u"W m-2"
    if field == 'olr':
        varfield=f2.createVariable('olr',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"OUTGOING LW RAD FLUX (TOA)"
        unitsname=u"W m-2"
    if field == 'p_1':
        varfield=f2.createVariable('p_1',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"PSTAR AFTER TIMESTEP"
        unitsname=u"Pa"
    if field == 'csolr':
        varfield=f2.createVariable('field207',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"CLEAR-SKY (II) UPWARD LW FLUX (TOA)"
        unitsname=u"W m-2"
    if field == 'longwave':
        varfield=f2.createVariable('longwave',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"NET DOWN SURFACE LW RAD FLUX"
        unitsname=u"W m-2"
    if field == 'solar':
        varfield=f2.createVariable('solar',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"NET DOWN SURFACE SW FLUX: SW TS ONLY"
        unitsname=u"W m-2"
    if field == 'field203':
        varfield=f2.createVariable('field203',np.float32,
                                   ('time','ht','latitude','longitude'))
        longname=u"TOTAL DOWNWARD SURFACE SW FLUX"
        unitsname=u"W m-2"
    if field == 'field208':
        varfield=f2.createVariable('field208',np.float32,
                              ('time','ht','latitude','longitude'))
        longname=u"CLEAR-SKY (II) DOWN SURFACE SW FLUX"
        unitsname=u"W m-2"
    if field == 'sh':
        varfield=f2.createVariable('sh',np.float32,
                              ('time','ht','latitude','longitude'))
        longname=u"SURFACE & B.LAYER HEAT FLUXES"
        unitsname=u"W m-2"
    if field == 'lh':
        varfield=f2.createVariable('lh',np.float32,
                              ('time','ht','latitude','longitude'))
        longname=u"SURFACE LATENT HEAT FLUX"
        unitsname=u"W m-2"

    if field == 'SST':
        varfield=f2.createVariable('SST',np.float32,
                              ('time','ht','latitude','longitude'),
                               fill_value = -99999.)
        longname=u"OCN TOP-LEVEL TEMPERATURE"
        unitsname=u"K"

    if field == 'SSS':
        varfield=f2.createVariable('SSS',np.float32,
                              ('time','ht','latitude','longitude'),
                               fill_value = -99999.)
        longname=u"OCN TOP-LEVEL SALINITY"
        unitsname=u"psu"

    if field == 'salinity':
        depth=f2.createDimension('depth',20)
        depth=f2.createVariable('depth',np.float32,('depth',))
        depth[:]=depthin
        varfield=f2.createVariable('salinity',np.float32,
                              ('time','depth','latitude','longitude'),
                               fill_value = -99999.)
        longname=u"SALINITY"
        unitsname=u"psu"

   
    # create variable attributes
    lons.setncatts({'units':u"degrees_east"}) 
    lats.setncatts({'units':u"degrees_north"})
    levels.setncatts({'units':u"m"})
    times.setncatts({'units':u"days since 0000-01-01 00.00",\
                         'calendar':"360_day"})
    varfield.setncatts({'long_name': longname,\
                                'units':unitsname})

    # assign data to variables
    lons[:]=lonin
    lats[:]=latin
    levels[:]=-1.0
    times[:]=0.0
    varfield[0,:,:,:]=avgvar
    
        
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

expts=['xpsie','xpsig','xpsid']
#expts=['xqbwd','xqbwe','xqbwg','xqbwh','xqbwi','xqbwj','xqbwk','xqbwl',
#       'xqbwm','xqbwn','xqbwo','xqbwp','xqbwq','xqbwr','xqbws','xqbwt']
startyear=1000
endyear=1200
HadCM3='y'

#HadCM3='y'
#expt='xibos'
#startyear=3500
#endyear=50

#field='p_1'
#field = 'temp'
#get_annual_data(expt,startyear,endyear,field)
#sys.exit(0)

for expt in expts:
    field='salinity'
    get_annual_data(expt,startyear,endyear,field)
sys.exit(0)

field='sh'
get_annual_data(expt,startyear,endyear,field)

field='lh'
get_annual_data(expt,startyear,endyear,field)



field='olr'
get_annual_data(expt,startyear,endyear,field)

field='field200'
get_annual_data(expt,startyear,endyear,field)

field='field201'
get_annual_data(expt,startyear,endyear,field)

#field='precip_1'
#get_annual_data(expt,startyear,endyear,field)


#temperature at 1.5m
if HadCM3=='y':
    field='temp'   # temperature at 1.5m in HadCM3
else:
    field='temp_1' # temperature at 1.5m in HadGEM
get_annual_data(expt,startyear,endyear,field)



field='field207'
get_annual_data(expt,startyear,endyear,field)

field='ilr'
get_annual_data(expt,startyear,endyear,field)
field='csolr'
get_annual_data(expt,startyear,endyear,field)

field='longwave'
get_annual_data(expt,startyear,endyear,field)



field='field203'
get_annual_data(expt,startyear,endyear,field)
field='solar'
get_annual_data(expt,startyear,endyear,field)

field='field30'
get_annual_data(expt,startyear,endyear,field)

field='field208'
get_annual_data(expt,startyear,endyear,field)

field='field207_1'
get_annual_data(expt,startyear,endyear,field)
    

sys.exit()
