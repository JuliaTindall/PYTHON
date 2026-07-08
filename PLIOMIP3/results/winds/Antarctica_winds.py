#NAME
#    Antarctica_winds
#PURPOSE 
# This program will plot the winds around Antarctica and the difference
# between them

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

years='00000301'  # second 10 years
#expts = ['xpsid','xpsio']
expts = ['xqbwd','xqbwo']

uvar = "u"
vvar = "v"
lat_name = "latitude"
lon_name = "longitude"


period = {'xqbwc':'PI', 'xqbwd':'LP','xqbwe':'EP400','xqbwg':'EP',
           'xqbwi':'LP280','xqbwj':'LP490'}

def recenter_lon(lon):
    lon = np.where(lon > 180, lon - 360, lon)
    return lon



#####################################################################


all_u = []
all_v = []
for expt in expts:
    exptu = []
    exptv = []
    allfiles =  sorted(Path('/home/earjcti/um/'+expt+'/pcpd/').
                     glob(expt + '*pc'+years+'*.nc'))

    #print(years)
    #print(expt)
    #print(allfiles)
    #sys.exit(0)

    for file in allfiles:
        ds_file = xr.open_dataset(file)
        u2 = ds_file[uvar].squeeze()
        v2 = ds_file[vvar].squeeze()
        u=u2.isel(p=0)
        v=v2.isel(p=0)
        lon = ds_file[lon_name]
        lat = ds_file[lat_name]
        exptu.append(u)
        exptv.append(v)

    u_mean = xr.concat(exptu, dim="file").mean(dim="file")
    v_mean = xr.concat(exptv, dim="file").mean(dim="file")
    
        
    #print(u_mean.values.shape)
    #print(u_mean[lon_name].values)
    #sys.exit(0)
    all_u.append(u_mean)
    all_v.append(v_mean)

#==================================
# plot mixed layer depth
plt.figure(figsize=(8, 12))

# set up for colorbar
#vals = [-200, -100, -50, -25, -10, -5, -2, 0, 2, 5, 10, 25, 50, 100, 200]
#vals = [val/100 for val in vals]
vals = np.arange(0,11,1)

cmap_base = plt.get_cmap('YlGnBu', len(vals) - 1) 
colors = [cmap_base(i) for i in range(cmap_base.N)]
#mid_index = vals.index(0.0) - 1  # bin just below zero
#colors[mid_index] = (1, 1, 1, 1)  # pure white RGBA
#colors[mid_index+1] = (1, 1, 1, 1)  # pure white RGBA
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

lon2d, lat2d = np.meshgrid(recenter_lon(u[lon_name].values),
                           u[lat_name].values)


for i,expt in enumerate(expts):
    u=all_u[i]
    v=all_v[i]
    wind_speed = np.sqrt(u.values**2 + v.values**2)

    #print('ustart',u.values[:3,:3],'uend')
    #print('vstart',v.values[:3,:3],'vend')
    #sys.exit(0)
    ax = plt.subplot(221+i,projection=ccrs.SouthPolarStereo())
    cf = ax.pcolormesh(lon2d, lat2d, wind_speed,
                       transform=ccrs.PlateCarree(), cmap=cmap,norm=norm,
                       alpha=0.5)

    skipx = 3
    skipy = 1
    q = ax.quiver(lon2d[::skipy,::skipx], lat2d[::skipy,::skipx],
              u.values[::skipy,::skipx], v.values[::skipy,::skipx],
                  transform=ccrs.PlateCarree(),scale=100,
                  width=0.003,headwidth=4,headlength=5)

    ax.coastlines()
    ax.quiverkey(q, 0.9, -0.1, 10, "10 m/s")
    
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.5,
        color='black',
        alpha=0.8,
        linestyle='--'
    )
    
    
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    plt.colorbar(cf, orientation="horizontal", label=f"m/s",extend='max')
    plt.title(expt + '-' + period.get(expt,expt))

#plt.show()
#sys.exit(0)

if len(expts) == 2:
    vals = np.arange(-4.0,4.4,0.4)

    cmap_base = plt.get_cmap('RdBu_r', len(vals) - 1) 
    colors = [cmap_base(i) for i in range(cmap_base.N)]
    #mid_index = vals.index(0.0) - 1  # bin just below zero
    #colors[mid_index] = (1, 1, 1, 1)  # pure white RGBA
    #colors[mid_index+1] = (1, 1, 1, 1)  # pure white RGBA
    cmap = mcolors.ListedColormap(colors)
    norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)


    print('doing difference between experiments')
    u_1 = all_u[0]
    u_2 = all_u[1]
    diff_u = u_2 - u_1
    v_1 = all_v[0]
    v_2 = all_v[1]
    diff_v = v_2 - v_1
    wind_speed =  (np.sqrt(u_1.values**2 + v_1.values**2)-np.sqrt(u_2.values**2 + v_2.values**2))

    ax = plt.subplot(223,projection=ccrs.SouthPolarStereo())
    cf = ax.pcolormesh(lon2d, lat2d, wind_speed,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
    skipx = 3
    skipy = 1
    q = ax.quiver(lon2d[::skipy,::skipx], lat2d[::skipy,::skipx],
              diff_u.values[::skipy,::skipx], diff_v.values[::skipy,::skipx],
                  transform=ccrs.PlateCarree(),scale=20)
    
    ax.coastlines()  
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=False,
        linewidth=0.5,
        color='black',
        alpha=0.8,
        linestyle='--' )
   
    ax.set_extent([-180, 180, -90, -60], ccrs.PlateCarree())
    plt.colorbar(cf, orientation="horizontal", label=f"Sv")
    ax.quiverkey(q, 0.9, -0.1, 1, "1 m/s")
   
    plt.title(expts[1] + '-' +  expts[0])

plt.show()
sys.exit(0)


fileout = "u_v_"+years + "*.png"
plt.savefig(fileout)
plt.close()

