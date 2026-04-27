# 
# plots density poleward of 60S with depth.  We need to already have produced
# basin_diagnostics/mean_xpsid_Pacific1400_1499.nc using program
# plot_density_by_basin.py
#

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# --- File setup ---
filestart = '/home/earjcti/um/xpsid/basin_diagnostics/mean_xpsid_Pacific'
filename1 = filestart + '1400_1499.nc'
filename2 = filestart.replace('xpsid', 'xpsie') + '1400_1499.nc'
filename3 = filestart.replace('xpsid', 'xpsie') + '1600_1699.nc'

filenames = [filename1, filename2, filename3]
legends   = ['LP_1400', 'EP400_1400', 'EP400_1600']

depth_name = 'depth_1'
lat_name   = 'latitude'
dens_name  = 'insitu_T_0'

ranges = [(0, 200), (200, 1000), (1000, 4000)]

fig, axes = plt.subplots(3, 1, figsize=(6, 9), sharex=True)

for i, file_path in enumerate(filenames):
    ds = xr.open_dataset(file_path)
    density = ds[dens_name]
    lat     = ds[lat_name]
    depth   = ds[depth_name]

    # Select latitudes south of 60S
    south_mask = lat <= -60
    density_south = density.sel({lat_name: south_mask})
    lat_south     = lat.sel({lat_name: south_mask})

    # Area weights
    weights = np.cos(np.deg2rad(lat_south))

    # Weighted average over latitude
    dens_w = density_south.weighted(weights).mean(dim=lat_name)

    # Plot each depth band in its own panel
    for ax, (dmin, dmax) in zip(axes, ranges):
        dcoord = dens_w[depth_name]
        mask   = (dcoord >= dmin) & (dcoord <= dmax)
        dens_band  = dens_w.where(mask, drop=True)
        depth_band = dcoord.where(mask, drop=True)

        ax.plot(dens_band.values, depth_band.values, label=legends[i], linewidth=2)
        #ax.set_title(f'{dmin}-{dmax} m', fontsize=12)
        ax.grid(True)
        ax.set_ylim(dmax, dmin)  # invert y-axis

# Add y-labels for each subplot
for ax in axes:
    ax.set_ylabel('Depth (m)', fontsize=10)

# Shared x-label
fig.text(0.5, 0.04, 'Density (kg/m³)', ha='center', fontsize=12)

# Legend in bottom panel
axes[-1].legend(title='Experiment', fontsize=10)

plt.tight_layout(rect=[0, 0.06, 1, 1])
plt.show()
