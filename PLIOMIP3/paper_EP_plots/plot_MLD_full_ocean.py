#NAME
#    plot_MLD  # mixed layer depth
#PURPOSE 
# This program will plot the mixed layer depth over the southern ocean to
# assess mixing between the surface and deeper layers

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
expts = ['xqbwj','xqbwg']

years= '0000039'  # 1400-1500 and 2000-210
fig = ['a)','b)']
#expt = 'xpsie'

mldvar = "field653"
lat_name = "latitude"
lon_name = "longitude"


period = {'xqbwc':'PI', 'xqbwd':'LP','xqbwe':'EP400','xqbwg':'EP',
           'xqbwi':'LP280','xqbwj':'LP$_{490}$'}

def recenter_lon(lon):
    lon = np.where(lon > 180, lon - 360, lon)
    return lon



#####################################################################


all_mld = []
for expt in expts:
    exptmld = []
    allfiles =  sorted(Path('/home/earjcti/um/'+expt+'/pg/').
                     glob(expt + '*pg'+years+'*.nc'))

    print(years,len(allfiles))
    for file in allfiles:
        ds_file = xr.open_dataset(file)
        mld = ds_file[mldvar].squeeze()
        lon = ds_file[lon_name]
        lat = ds_file[lat_name]
        exptmld.append(mld)

    mld_mean = xr.concat(exptmld, dim="file").mean(dim="file")
    
        
    #print(mld_mean.values.shape)
    #print(mld_mean[lon_name].values)
    #sys.exit(0)
    all_mld.append(mld_mean)

#==================================
# plot mixed layer depth
plt.figure(figsize=(8, 8))

# set up for colorbar
#vals = [-200, -100, -50, -25, -10, -5, -2, 0, 2, 5, 10, 25, 50, 100, 200]
#vals = [val/100 for val in vals]
vals = np.arange(0,220,20)

cmap_base = plt.get_cmap('viridis', len(vals) - 1) 
colors = [cmap_base(i) for i in range(cmap_base.N)]
#mid_index = vals.index(0.0) - 1  # bin just below zero
#colors[mid_index] = (1, 1, 1, 1)  # pure white RGBA
#colors[mid_index+1] = (1, 1, 1, 1)  # pure white RGBA
cmap = mcolors.ListedColormap(colors)
cmap.set_bad('lightgrey')
norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

lon2d, lat2d = np.meshgrid(lon.values,lat.values)

#lon2d, lat2d = np.meshgrid(recenter_lon(mld[lon_name].values),
#                           mld[lat_name].values)


for i,expt in enumerate(expts):
    mld=all_mld[i].where(all_mld[i] > 0)
    #mld=all_mld[i]
    ax = plt.subplot(221+i,projection=ccrs.PlateCarree())
    #cf = ax.pcolormesh(lon2d, lat2d, mld.values,
    #               transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
    vals = np.arange(0,210,10)
    cf = plt.contourf(lon2d, lat2d, mld.values,
                      transform=ccrs.PlateCarree(), cmap=cmap, levels=vals,
                      extend='max')
    
    ax.coastlines()
    #ax.set_global()
    
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.5,
        color='black',
        alpha=0.8,
        linestyle='--'
    )
    
    
    #ax.set_extent([0, 360, -90, 90], ccrs.PlateCarree())
    cbar=plt.colorbar(cf, orientation="horizontal",extend='max')
    cbar.ax.tick_params(labelsize=12)
    cbar.ax.set_title('metres', fontsize=14)
    plt.title(fig[i] + ' MLD : ' +  period.get(expt),fontsize=16)


# plot anomaly
vals = np.arange(-35,40,10)

cmap_base = plt.get_cmap('RdBu_r', len(vals) - 1) 
colors = [cmap_base(i) for i in range(cmap_base.N)]
mid = len(colors)//2
#colors[mid-1] = (1,1,1,1)
colors[mid] = (1,1,1,1)
cmap = mcolors.ListedColormap(colors)
cmap.set_bad('lightgrey')
norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

print('doing difference between experiments')
mld_1 = all_mld[0].where(mld > 0)
mld_2 = all_mld[1].where(mld > 0)
diff = mld_2 - mld_1
ax = plt.subplot(223,projection=ccrs.PlateCarree())
#cf = ax.pcolormesh(lon2d, lat2d, diff.values,
#                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
cf = plt.contourf(lon2d, lat2d, diff.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,levels=vals,
                   extend='both')

ax.coastlines()  
gl = ax.gridlines(
    crs=ccrs.PlateCarree(),
    draw_labels=False,
    linewidth=0.5,
    color='black',
    alpha=0.8,
    linestyle='--' )

#ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
cbar=plt.colorbar(cf, orientation="horizontal",extend='both')
cbar.ax.tick_params(labelsize=12)
cbar.ax.set_title('metres', fontsize=14)
plt.title('c) MLD :' + period.get(expts[1]) + ' - ' + period.get(expts[0]),fontsize=16)

fileout = 'mld_full.png'
plt.savefig(fileout)
plt.close()

