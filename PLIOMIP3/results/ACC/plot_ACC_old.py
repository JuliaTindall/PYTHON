#NAME
#    plot_ACC
#PURPOSE 
# This program calculates Antarctic Circumpolar Current (ACC) transport through
# a defined Drake Passage section.  
# It reads annual NetCDF files for both an experiment and control run 
# It produces a time series of transport for each case and their difference in Sverdrups.  
# It also creates a polar stereographic map of the 100‑year mean surface zonal velocity.  

# This program will plot the Antarctic circumpolar current
# Import necessary libraries

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from pathlib import Path
import cartopy.crs as ccrs
import sys



# === CONFIG ===
years='0000039[7-9]'
xqbwd_files = sorted(Path("/home/earjcti/um/xqbwd/pg/").
                     glob("xqbwd*pg"+years+"*.nc"))
xqbwc_files = sorted(Path("/home/earjcti/um/xqbwc/pg/").
                     glob("xqbwc*pg"+years+"*.nc"))
xqbwe_files = sorted(Path("/home/earjcti/um/xqbwe/pg/").
                     glob("xqbwe*pg"+years+"*.nc"))
xqbwg_files = sorted(Path("/home/earjcti/um/xqbwg/pg/").
                     glob("xqbwg*pg"+years+"*.nc"))

uvar = "field703"
lat_name = "latitude_1"
lon_name = "longitude_1"
lev_name = "depth_1"

lon0 = -68.0      # section longitude
lat_min = -75.0
lat_max = -52.0

def recenter_lon(lon):
    lon = np.where(lon > 180, lon - 360, lon)
    return lon

def layer_thickness(depth):
    z = depth.values
    edges = np.zeros(z.size + 1)
    edges[1:-1] = 0.5 * (z[1:] + z[:-1])
    edges[0] = z[0] - (edges[1] - z[0])
    edges[-1] = z[-1] + (z[-1] - edges[-2])
    return np.diff(edges)

def section_transport(u, lat, depth, lon, lon0):
    """
    Computes the total volume transport across a defined ocean section 
    using model velocity data.  
    Integrates zonal velocity over depth and width to return transport in 
    Sverdrups for each dataset.  
    """
    usq=u.squeeze()
    usq = usq / 100.  # convert from cm/s to m/s

    # Find nearest longitude
    lon = recenter_lon(lon)
    # index of required longitude
    i = np.argmin(np.abs(lon - lon0))

    # get section we care about (ie between latmin and latmax at lon0)
    u_sec = usq.isel({lon_name: i}).sel({lat_name: slice(lat_min, lat_max)})

    dy = np.abs(np.gradient(np.deg2rad(u_sec[lat_name].values))) * 6371000.0
    dz = layer_thickness(depth)

    term = u_sec * dz[:, None] * dy[None, :]
    #print(u_sec.shape)
    #print(dz.shape)
    #print(dz[:,None])
    #print(dy[None,:])
    #print(dy.shape)
    #print(u_sec.data[3,10],dz[3,None],dy[None,10])
    #print('term is',term.data[3,10],term.data[3,10]/1.0E6)
    #print(term.data/1E6)
    #print(np.nansum(term.data)/1E6)
    trans = term.sum(dim=(lev_name, lat_name),skipna=True) / 1e6  # Sv
    return float(trans)

def mean_depth_integrated_zonal_transport(files, uvar, lev_name, lat_name,
                                          lon_name):
    """
    mean depth integrated zonal transport from u
    returns trans which is transport
    """

    acc = None
    acc_surf = None
    n = 0
    for f in files:
        ds = xr.open_dataset(f)
        u = ds[uvar] / 100.  # convert to m/s
        
        # get section we care about (ie between latmin and latmax)
        u_sec = u.squeeze().sel({lat_name: slice(-90.0, -30.0)})

        # Depth-integrate u: sum_z u * dz  -> units m^2/s
        dz = layer_thickness(ds[lev_name])
        uint = (u_sec * dz[:,None,None]).sum(dim=lev_name)
        uint_surf = (u_sec[0,:,:] * dz[0,None,None])
    
        # Convert to transport through a meridional wall segment by multiplying by dy(lat)
        dy = np.abs(np.gradient(np.deg2rad(u_sec[lat_name].values))) * 6371000.0
    
        trans = (uint * dy[:,None]).rename("zonal_transport_Sv") / 1e6  # Sv
        trans_surf = (uint_surf * dy[:,None]).rename("zonal_transport_Sv") / 1e6  # Sv

        # Accumulate
        acc = trans if acc is None else acc + trans
        acc_surf = trans_surf if acc_surf is None else acc_surf + trans_surf
        n += 1

    mean_trans = acc / n
    mean_trans_surf = acc_surf / n
    
    lon = mean_trans[lon_name]
    lon = xr.where(lon > 180, lon - 360, lon) 
    mean_trans = mean_trans.assign_coords({lon_name: lon}).sortby(lon_name)
    mean_trans_surf = mean_trans_surf.assign_coords({lon_name: lon}).sortby(lon_name)
    return mean_trans,mean_trans_surf


#####################################################################
# === Loop over years ===

trans_xqbwd, trans_xqbwc, trans_xqbwg, trans_xqbwe = [], [], [], []

for fxqbwd, fxqbwc, fxqbwg,fxqbwe in zip(xqbwd_files,
                                         xqbwc_files,
                                         xqbwg_files,
                                         xqbwe_files):
    ds_xqbwd = xr.open_dataset(fxqbwd)
    ds_xqbwc = xr.open_dataset(fxqbwc)
    ds_xqbwg = xr.open_dataset(fxqbwg)
    ds_xqbwe = xr.open_dataset(fxqbwe)
    t_xqbwd = section_transport(ds_xqbwd[uvar], ds_xqbwd[lat_name], ds_xqbwd[lev_name], ds_xqbwd[lon_name], lon0)
    t_xqbwc = section_transport(ds_xqbwc[uvar], ds_xqbwc[lat_name], ds_xqbwc[lev_name], ds_xqbwc[lon_name], lon0)
    t_xqbwe = section_transport(ds_xqbwe[uvar], ds_xqbwe[lat_name], ds_xqbwe[lev_name], ds_xqbwe[lon_name], lon0)
    t_xqbwg = section_transport(ds_xqbwg[uvar], ds_xqbwg[lat_name], ds_xqbwg[lev_name], ds_xqbwg[lon_name], lon0)
    trans_xqbwg.append(t_xqbwg)
    trans_xqbwe.append(t_xqbwe) 
    trans_xqbwd.append(t_xqbwd)
    trans_xqbwc.append(t_xqbwc)  # for the prindustrial Ai says this should be
                             # 130-150Sv

trans_xqbwd = np.array(trans_xqbwd)
trans_xqbwc = np.array(trans_xqbwc)
trans_xqbwg = np.array(trans_xqbwg)
trans_xqbwe = np.array(trans_xqbwe)
diff = trans_xqbwd - trans_xqbwc

print(f"Mean ACC transport (xqbwd): {trans_xqbwd.mean():.2f} Sv")
print(f"Mean ACC transport (xqbwc): {trans_xqbwc.mean():.2f} Sv")
print(f"Mean ACC transport (xqbwg): {trans_xqbwg.mean():.2f} Sv")
print(f"Mean ACC transport (xqbwe): {trans_xqbwe.mean():.2f} Sv")
print(f"Mean difference: {diff.mean():.2f} Sv")

# === Time series plot ===

plt.figure()
plt.plot(trans_xqbwg, label="xqbwg-EP")
plt.plot(trans_xqbwe, label="xqbwe-EP400")
plt.plot(trans_xqbwd, label="xqbwd-LP")
plt.plot(trans_xqbwc, label="xqbwc_PI")
#plt.plot(diff, label="Xqbwd - xqbwc", ls="--")
plt.xlabel("Year")
plt.ylabel("Transport (Sv)")
plt.title("ACC Transport at Drake Passage")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("acc_transport_timeseries.png", dpi=150)
plt.close()

# === get ACC zonally integrated and at the surface===
(transport_xqbwd,
surftran_xqbwd) = mean_depth_integrated_zonal_transport(xqbwd_files,uvar,
                                                       lev_name, lat_name,
                                                       lon_name)
(transport_xqbwc,
 surftran_xqbwc)= mean_depth_integrated_zonal_transport(xqbwc_files,uvar,
                                                       lev_name, lat_name,
                                                       lon_name)

(transport_xqbwe,
 surftran_xqbwe)= mean_depth_integrated_zonal_transport(xqbwe_files,uvar,
                                                       lev_name, lat_name,
                                                       lon_name)
(transport_xqbwg,
 surftran_xqbwg)= mean_depth_integrated_zonal_transport(xqbwg_files,uvar,
                                                       lev_name, lat_name,
                                                       lon_name)
lon2d, lat2d = np.meshgrid(recenter_lon(transport_xqbwd[lon_name].values),
                           transport_xqbwd[lat_name].values)


#==================================
# plot depth integrated field
plt.figure(figsize=(8, 8))

# set up for colorbar
vals = [-200, -100, -50, -25, -10, -5, -2, 0, 2, 5, 10, 25, 50, 100, 200]
#vals = [val/100 for val in vals]
cmap_base = plt.get_cmap('RdBu_r', len(vals) - 1) 
colors = [cmap_base(i) for i in range(cmap_base.N)]
mid_index = vals.index(0) - 1  # bin just below zero
colors[mid_index] = (1, 1, 1, 1)  # pure white RGBA
colors[mid_index+1] = (1, 1, 1, 1)  # pure white RGBA
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

ax = plt.subplot(221,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, transport_xqbwc.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwc - PI")

ax = plt.subplot(222,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, transport_xqbwd.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwd LP")

ax = plt.subplot(223,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, transport_xqbwe.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwe  EP400")

ax = plt.subplot(224,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, transport_xqbwg.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwg EP")


fileout = "depth_integrated_zonal_transport_"+years + "*.png"
plt.savefig(fileout)
plt.close()



#==================================
# plot surface integrated field
plt.figure(figsize=(8, 12))

# set up for colorbar
#vals = [-200, -100, -50, -25, -10, -5, -2, 0, 2, 5, 10, 25, 50, 100, 200]
#vals = [val/100 for val in vals]
vals = np.arange(-0.5,0.55,0.05)

cmap_base = plt.get_cmap('RdBu_r', len(vals) - 1) 
colors = [cmap_base(i) for i in range(cmap_base.N)]
#mid_index = vals.index(0.0) - 1  # bin just below zero
#colors[mid_index] = (1, 1, 1, 1)  # pure white RGBA
#colors[mid_index+1] = (1, 1, 1, 1)  # pure white RGBA
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

ax = plt.subplot(321,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, surftran_xqbwc.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwc - PI")

ax = plt.subplot(322,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, surftran_xqbwd.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwd LP")

ax = plt.subplot(323,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, surftran_xqbwe.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwe  EP400")

ax = plt.subplot(324,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, surftran_xqbwg.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwg EP")


vals = np.arange(-0.1,0.11,0.01)

cmap_base = plt.get_cmap('RdBu_r', len(vals) - 1) 
colors = [cmap_base(i) for i in range(cmap_base.N)]
#mid_index = vals.index(0.0) - 1  # bin just below zero
#colors[mid_index] = (1, 1, 1, 1)  # pure white RGBA
#colors[mid_index+1] = (1, 1, 1, 1)  # pure white RGBA
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(boundaries=vals, ncolors=cmap.N)

ax = plt.subplot(325,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, surftran_xqbwg.values - surftran_xqbwe.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwg - xqbwe EP-EP400")

ax = plt.subplot(326,projection=ccrs.SouthPolarStereo())
cf = ax.pcolormesh(lon2d, lat2d, surftran_xqbwg.values - surftran_xqbwc.values,
                   transform=ccrs.PlateCarree(), cmap=cmap,norm=norm)
ax.coastlines()
ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
plt.colorbar(cf, orientation="horizontal", label=f"Sv")
plt.title("xqbwg - xqbwc EP-PI")


fileout = "top_layer_zonal_transport_"+years + "*.png"
plt.savefig(fileout)
plt.close()

