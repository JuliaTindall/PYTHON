#!/usr/bin/env python2.7
#NAME
#
# This program will plot the density at different levels
#
#
import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import iris
import iris.quickplot as qplt
import cartopy.crs as ccrs
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid
import subprocess
import pandas as pd



def get_avg(year, latstart):
    """
    gets average salinity for this year between latstart and -90
    """
 
    yearuse = str(year).zfill(9)
    filename=('/nfs/hera1/earjcti/um/'+exptname+'/pg/'+exptname+'o#pg'
              + yearuse + 'c1+.nc')
    cube_alllevs = iris.load_cube(filename,
                                  'SALINITY (OCEAN)       (PSU-35)/1000')
    cube = cube_alllevs[0,0,:,:]  # this is surface salinity
    cube=cube * 1000. + 35.0
    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    weights = iris.analysis.cartography.area_weights(cube) # this is in m2
  
    for j,lat in enumerate(cube.coord('latitude').points):
        if lat > latstart:
            weights[j,:]=0.0
      
    SH_sal_mean = cube.collapsed(['latitude','longitude'],
                                 iris.analysis.MEAN,weights=weights)

  
    return (SH_sal_mean.data)

#####################################################################
def plotdrifts(salinity,startyear,endyear,latstart):
    """
    plots the timeseries of seaice area
    file
    """

    plt.subplot(1,1,1)
    plt.plot(salinity)
    plt.ylim(33.0,34.5)
    plt.title('Salinity ' + P3name.get(exptname))
    plt.ylabel('psu')
    plt.xlabel('year')
    if latstart < 0:
        latstartuse = str(latstart * -1.0) + 'S'
    else:
        latstartuse = str(latstart) + 'N'


    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/salinity_'+exptname+'_' + latstartuse + '-90S_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.eps') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    fileout=('/nfs/hera1/earjcti/um/'+exptname+'/spinup/salinity_'+exptname+'_' + latstartuse + '-90S_' + np.str(np.int(startyear)) + '_'+ np.str(np.int(endyear)) +'.png') 
    print('savingfig',fileout)
    plt.savefig(fileout, bbox_inches='tight')  

    
    plt.close()


   


#######################################################
def get_SHsalin(HadCM3,exptname,startyear,endyear,latstart):
    """
    reads in the salinity and extracts for the SH
    """

    # arrays for storing seaice
    salinity_SH  = np.zeros(endyear-startyear+1)

    # obtain means for each year and store in the arrays
    
    for year in range(startyear,endyear+1):
        print(year)
        salin = get_avg(year,latstart)
        salinity_SH[year-startyear]=salin

    # plot and save
    plotdrifts(salinity_SH,startyear,endyear,latstart)
  


################################
# main program

# annual mean
figureno=0


P3name = {'xpsie' : 'EP400',    'xpsig':'EP490' ,'xpsig':'LP'}

latstart=-70   # will plot from latstart to -90.0
HadCM3='y'
exptname='xpsig'
startyear=12  # can't start before year 12 because we aren't outputting d18o
endyear=1999

filename=('/home/earjcti/um/'+ exptname + '/basin_diagnostics/'
          + exptname + '_Pacific_density' + str(startyear) + '_' +
          str(endyear) + '_' + str(latstart) + '.txt')


df = pd.read_csv(filename)
print(df.columns.tolist())

# Access columns
years = df["Year"].tolist()
dens_surf = df[" dens_surf"].tolist()
dens_0_100 = df["dens_0_100"].tolist()
dens_0_200 = df["dens_0_200"].tolist()
dens_0_1000 = df["dens_0_1000"].tolist()
dens_1000_3000 = df["dens_1000_3000"].tolist()


# print mean value from year 500-1000

mask = (df["Year"] >= 500) & (df["Year"] <= 1000)
df_subset = df.loc[mask]

if df_subset.empty:
    print("No rows with Year between 100 and 500.")
else:
    # Select numeric columns excluding Year
    numeric_cols = df_subset.select_dtypes(include=[np.number]).columns.tolist()
  
    # Compute means
    means = df_subset[numeric_cols].mean()

    # Print to display (pretty)
    print("Means for Year in [500,1000]:")
    for col, val in means.items():
        print(f"{col}: {val:.6f}")


        
# plot timeseries
fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,figsize=(11,8),
                              sharex=True,gridspec_kw=dict(hspace=0.25))

ax1.plot(years,dens_surf,label='surface density')
ax1.plot(years,dens_0_100,label='0-100m')
ax1.plot(years,dens_0_200,label='0-200m')
ax1.plot(years,dens_0_1000,label='0-1000m')
ax1.plot(years,dens_1000_3000,label='1000-3000m')

ax1.set_ylim(35.5,37.0)
ax1.set_xlabel("Year")
ax1.set_ylabel("density")
ax1.set_title("density at different levels: " + exptname)
ax1.grid(True, alpha=0.3)
ax1.legend()  # adjust ncols if many lines


ax2.plot(years,[b-a for b,a in zip(dens_0_100,dens_surf)],
         label='anom 0-100 minus surface density')
ax2.plot(years,[b-a for b,a in zip(dens_0_200,dens_0_100)],
         label='anom 0-200m - 0-100')
ax2.plot(years,[b-a for b,a in zip(dens_0_1000,dens_0_200)],
         label='anom 0-1000 - 0-200m')
ax2.plot(years,[b-a for b,a in zip(dens_1000_3000,dens_0_1000)],
         label='anom 1000-3000 - 0-1000m')

ax2.set_xlabel("Year")
ax2.set_ylabel("density")
ax2.set_title("1000-3000m density minus density at different levels: " + exptname)
ax2.set_ylim(0.,0.4)
#ax2.set_xlim(500,750)
ax2.grid(True, alpha=0.3)
ax2.legend(loc='lower left')  # adjust ncols if many lines

plt.tight_layout()
plt.show()
####

