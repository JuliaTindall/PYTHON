#!/usr/bin/env python2.7
#NAME
#    TOA_rad_inbalance
#PURPOSE
#    This program was based on PYTHON/PROGRAMS/HadGEM2/plot_radiation.py
#
#  It should do the same thing but might have been updated for PlioMIP3.
#
# search for 'main program' to find end of functions
# Julia 22/11/2016



import os
import numpy as np
#import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
#import cartopy.crs as ccrs
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
#import subprocess


#functions are:
#  def annmean

def get_rad_inbal(year):
    # we will get  a) incoming sw ra flux (toa) field200
    #              b) outgoing sw ra flux (toa) field201
    #              c) outgoing lw rad flux (toa) olr
    # other things if the budgets don't balance

    if HadCM3 == 'y':
        datasetname='/uolstore/Research/a/hera1/earjcti/um/'+exptname+'/pcpd/'+exptname+'a#pd' + str(year).zfill(9) + '??+.nc'
        print(datasetname)
        #sys.exit(0)
        f=MFDataset(datasetname)
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        lat = f.variables['latitude'][:]
        lon = f.variables['longitude'][:]
        insw=f.variables['field200'][:]
        outsw=f.variables['field201'][:]
        outlw=f.variables['olr'][:]
        netlw= (-1.0) *outlw
        netsw=insw-outsw


    netsw=np.squeeze(netsw)
    netlw=np.squeeze(netlw)
    ntimes,ny,nx=np.shape(netsw)
    print(ntimes,ny,nx)
    
#average across the time dimension
    netsw_ann=np.mean(netsw,axis=0)
    netlw_ann=np.mean(netlw,axis=0)
    
    
    # calculate residual and mean residual weighted by cos latitude

    residual=netsw_ann+ netlw_ann
  
    weights=np.cos(np.radians(lat))
    resid_zon=np.average(residual,axis=0,weights=weights)
    average_residual=np.average(resid_zon)


    return(average_residual)


################################
# main program

# annual mean
figureno=0


HadCM3='y'
startyear=1991
endyear=3001

#endings=['i','j','k','l','m','n','o','p','q','r','s','t']
endings=['j']
for ending in endings:
    #exptname = 'xqbw' + ending
    exptname = 'xpsi' + ending
    ann_avg_ts = np.zeros(endyear-startyear+1)
    years=np.arange(startyear, endyear+1)

    try:
        for year in range(startyear,endyear+1):
            print(year)
            ann_avg = get_rad_inbal(year)
            ann_avg_ts[year-startyear]=ann_avg

        print('got data')
       # get running mean (30 year)

        weights = np.ones(30)/30
        #print('sum weights',np.sum(weights))
        run_mean = np.convolve(ann_avg_ts, weights, "valid")
        years_run_mean = np.arange(startyear+14, endyear-14)

#print(run_mean)
#print(np.shape(run_mean))
#print(np.shape(years_run_mean))
#sys.exit(0)

        print('about to plot')
        plt.plot(years,ann_avg_ts)
        print('j1')
        plt.plot(years_run_mean, run_mean)
        print('j2')
        plt.hlines(y=0,xmin=startyear, xmax=endyear)
        plt.title('TOA radiation inbalance:'+exptname)
        plt.xlabel('year')
        plt.ylabel('W/m2')
        fileout = '/uolstore/Research/a/hera1/earjcti/um/' + exptname + '/spinup/' + exptname + '_TOAinbal.eps'
        plt.savefig(fileout)
        plt.close()


        datafile = '/uolstore/Research/a/hera1/earjcti/um/' + exptname + '/spinup/' + exptname + '_TOAinbal.tex'
        print('about to open',datafile)
        with open(datafile, 'w') as f:
            f.write('% year  toa_imbalance_Wm2\n')
            for yr, val in zip(years, ann_avg_ts):
                f.write(f'{yr:d} {val:.6f}\n')
      

    except:
        print('failure on',exptname)



sys.exit(0)

####

