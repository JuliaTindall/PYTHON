#NAME
#    plot_SAL  # salinity
#PURPOSE 
# This program will plot the salinity over the southern ocean to compare
# with MLD plots

import xarray as xr
#import xesmf as xe # for regridding
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import cartopy.crs as ccrs
import sys

# === CONFIG ===
#years='0000039[7-9]'
#expts = ['xqbwg','xqbwd','xqbwe','xqbwg']

years='00000202'  # second 10 years
expts = ['xpsid','xpsio']

salvar = "salinity"
lat_name = "latitude"
lon_name = "longitude"


period = {'xqbwc':'PI', 'xqbwd':'LP','xqbwe':'EP400','xqbwg':'EP',
           'xqbwi':'LP280','xqbwj':'LP490'}

def recenter_lon(lon):
    lon = np.where(lon > 180, lon - 360, lon)
    return lon



#####################################################################


all_sal = []
for expt in expts:
    exptsal = []
    allfiles =  sorted(Path('/home/earjcti/um/'+expt+'/pg/').
                     glob(expt + '*pg'+years+'*.nc'))

    for file in allfiles:
        ds_file = xr.open_dataset(file)
        sal = ds_file[salvar].isel(depth_1=0).squeeze()
        sal = (sal * 1000.) + 35.
        lon = ds_file[lon_name]
        lat = ds_file[lat_name]
        exptsal.append(sal)

    sal_mean = xr.concat(exptsal, dim="file").mean(dim="file")
    
        
    #print(sal_mean.values.shape)
    #print(sal_mean[lon_name].values)
    #sys.exit(0)
    all_sal.append(sal_mean)

#==================================
# plot mixed layer depth
plt.figure(figsize=(8, 12))

# set up for colorbar
#vals = [-200, -100, -50, -25, -10, -5, -2, 0, 2, 5, 10, 25, 50, 100, 200]
#vals = [val/100 for val in vals]
vals = np.arange(32.5,36.0,0.5)

cmap_base = plt.get_cmap('viridis', len(vals) - 1) 
colors = [cmap_base(i) for i in range(cmap_base.N)]
#mid_index = vals.index(0.0) - 1  # bin just below zero
#colors[mid_index] = (1, 1, 1, 1)  # pure white RGBA
#colors[mid_index+1] = (1, 1, 1, 1)  # pure white RGBA
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

lon2d, lat2d = np.meshgrid(recenter_lon(sal[lon_name].values),
                           sal[lat_name].values)


for i,expt in enumerate(expts):
    sal=all_sal[i]
    ax = plt.subplot(221+i,projection=ccrs.SouthPolarStereo())
    cf = ax.pcolormesh(lon2d, lat2d, sal.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
    ax.coastlines()
    
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.5,
        color='black',
        alpha=0.8,
        linestyle='--'
    )
    
    
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    plt.colorbar(cf, orientation="horizontal", label=f"psu",extend='both')
    plt.title(expt + '-' + period.get(expt,expt))

    
if len(expts) == 2:
    vals = np.arange(-1,1.1,0.1)

    cmap_base = plt.get_cmap('RdBu_r', len(vals) - 1) 
    colors = [cmap_base(i) for i in range(cmap_base.N)]
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

    print('doing difference between experiments')
    sal_1 = all_sal[0]
    sal_2 = all_sal[1]
    diff = sal_2 - sal_1
    ax = plt.subplot(223,projection=ccrs.SouthPolarStereo())
    cf = ax.pcolormesh(lon2d, lat2d, diff.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
    ax.coastlines()  
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.5,
        color='black',
        alpha=0.8,
        linestyle='--' )
   
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    plt.colorbar(cf, orientation="horizontal", label=f"psu",extend='both')
    plt.title(expts[1] + '-' +  expts[0])
    
    
fileout = "sal_"+years + "*.png"
plt.savefig(fileout)
plt.close()

