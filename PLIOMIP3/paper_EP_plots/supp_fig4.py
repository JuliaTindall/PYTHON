#python program
#
# all plots for figure 4 of the paper
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import iris
#import iris.quickplot as qplt
import iris.plot as iplt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys


##################################################################
# set up figure and subplot
fig = plt.figure(figsize=[10.5,9])

spec = gridspec.GridSpec(ncols=3, nrows=2,
                         width_ratios=[1, 1,1], wspace=0.5,
                         hspace=0.5, height_ratios=[1,1])


period = {'xqbwd':'LP','xqwbj':'LP$_{490}','xqbwe':'EP$_{400}$',
          'xqbwg':'EP','xqbwc':'PI'}
filestart = '/home/earjcti/um/'

# read in data
f_ep400 = (filestart + 'xpsie/basin_diagnostics/mean_xpsie_Pacific2_51.nc')
f_ep = (filestart + 'xpsig/basin_diagnostics/mean_xpsig_Pacific2_51.nc')
f_lp = (filestart + 'xpsid/basin_diagnostics/mean_xpsid_Pacific2_51.nc')

#density EP400-LP

ax0=fig.add_subplot(spec[0])
EP400_cube_a = iris.load_cube(f_ep400,'density basin')
constraint1 = iris.Constraint(latitude = lambda cell: -80 < cell <  -60)
constraint2 = iris.Constraint(depth_1 = lambda cell: cell < 2200)
EP400_cube_b = EP400_cube_a.extract(constraint1)
EP400_cube = EP400_cube_b.extract(constraint2)
       
LP_cube_a = iris.load_cube(f_lp,'density basin')
LP_cube_b = LP_cube_a.extract(constraint1)
LP_cube =  LP_cube_b.extract(constraint2)

anom_cube = EP400_cube - LP_cube

vals = np.arange(-0.1,0.11,0.01)
cs=iplt.contourf(anom_cube,levels=vals,extend='both',cmap='RdBu_r')
#plt.xlabel('Latitude (degrees)',fontsize=14)
plt.ylabel('Depth (m)',fontsize=10)
ax0.tick_params(labelsize=10)
plt.title('a) Density: EP$_{400}$ - LP',fontsize=14)
    
cb=plt.colorbar(cs,orientation='vertical')
cb.ax.set_xlabel('kg m$^{-3}$',fontsize=12)
cb.ax.tick_params(labelsize=12)
cb.ax.xaxis.set_label_position('bottom')
cb.ax.xaxis.labelpad=15


#salinity EP400-LP

ax0=fig.add_subplot(spec[1])
EP400_cube_a = iris.load_cube(f_ep400,'salinity basin')
constraint1 = iris.Constraint(latitude = lambda cell: -80 < cell <  -60)
constraint2 = iris.Constraint(depth_1 = lambda cell: cell < 2200)
EP400_cube_b = EP400_cube_a.extract(constraint1)
EP400_cube = EP400_cube_b.extract(constraint2)
       
LP_cube_a = iris.load_cube(f_lp,'salinity basin')
LP_cube_b = LP_cube_a.extract(constraint1)
LP_cube =  LP_cube_b.extract(constraint2)

anom_cube = EP400_cube - LP_cube

vals = np.arange(-0.2,0.22,0.02)
cs=iplt.contourf(anom_cube,levels=vals,extend='both',cmap='RdBu_r')
#plt.xlabel('Latitude (degrees)',fontsize=14)
#plt.ylabel('Depth (m)',fontsize=10)
ax0.set_yticks([])
plt.title('b) Salinity: EP$_{400}$ - LP',fontsize=14)
    
cb=plt.colorbar(cs,orientation='vertical')
cb.ax.set_xlabel('psu',fontsize=12)
cb.ax.tick_params(labelsize=12)
cb.ax.xaxis.set_label_position('bottom')
cb.ax.xaxis.labelpad=15

#temperature EP400-LP

ax0=fig.add_subplot(spec[2])
EP400_cube_a = iris.load_cube(f_ep400,'temperature basin')
constraint1 = iris.Constraint(latitude = lambda cell: -80 < cell <  -60)
constraint2 = iris.Constraint(depth_1 = lambda cell: cell < 2200)
EP400_cube_b = EP400_cube_a.extract(constraint1)
EP400_cube = EP400_cube_b.extract(constraint2)
       
LP_cube_a = iris.load_cube(f_lp,'temperature basin')
LP_cube_b = LP_cube_a.extract(constraint1)
LP_cube =  LP_cube_b.extract(constraint2)

anom_cube = EP400_cube - LP_cube

vals = np.arange(-0.5,0.55,0.05)
cs=iplt.contourf(anom_cube,levels=vals,extend='both',cmap='RdBu_r')
#plt.xlabel('Latitude (degrees)',fontsize=14)
#plt.ylabel('Depth (m)',fontsize=10)
ax0.tick_params(labelsize=10)
plt.title('c) Temperature: EP$_{400}$ - LP',fontsize=14)
ax0.set_yticks([])

cb=plt.colorbar(cs,orientation='vertical')
cb.ax.set_xlabel('$^\circ$C',fontsize=12)
cb.ax.tick_params(labelsize=12)
cb.ax.xaxis.set_label_position('bottom')
cb.ax.xaxis.labelpad=15

# now do EP-LP

#density EP-LP

ax0=fig.add_subplot(spec[3])
EP_cube_a = iris.load_cube(f_ep,'density basin')
constraint1 = iris.Constraint(latitude = lambda cell: -80 < cell <  -60)
constraint2 = iris.Constraint(depth_1 = lambda cell: cell < 2200)
EP_cube_b = EP_cube_a.extract(constraint1)
EP_cube = EP_cube_b.extract(constraint2)
       
LP_cube_a = iris.load_cube(f_lp,'density basin')
LP_cube_b = LP_cube_a.extract(constraint1)
LP_cube =  LP_cube_b.extract(constraint2)

anom_cube = EP_cube - LP_cube

vals = np.arange(-0.1,0.11,0.01)
cs=iplt.contourf(anom_cube,levels=vals,extend='both',cmap='RdBu_r')
plt.xlabel('Latitude ($^\circ$)',fontsize=10)
plt.ylabel('Depth (m)',fontsize=10)
ax0.tick_params(labelsize=10)
plt.title('d) Density: EP - LP',fontsize=14)
    
cb=plt.colorbar(cs,orientation='vertical')
cb.ax.set_xlabel('kg m$^{-3}$',fontsize=12)
cb.ax.tick_params(labelsize=12)
cb.ax.xaxis.set_label_position('bottom')
cb.ax.xaxis.labelpad=15


#salinity EP-LP

ax0=fig.add_subplot(spec[4])
EP_cube_a = iris.load_cube(f_ep,'salinity basin')
constraint1 = iris.Constraint(latitude = lambda cell: -80 < cell <  -60)
constraint2 = iris.Constraint(depth_1 = lambda cell: cell < 2200)
EP_cube_b = EP_cube_a.extract(constraint1)
EP_cube = EP_cube_b.extract(constraint2)
       
LP_cube_a = iris.load_cube(f_lp,'salinity basin')
LP_cube_b = LP_cube_a.extract(constraint1)
LP_cube =  LP_cube_b.extract(constraint2)

anom_cube = EP_cube - LP_cube

vals = np.arange(-0.2,0.22,0.02)
cs=iplt.contourf(anom_cube,levels=vals,extend='both',cmap='RdBu_r')
plt.xlabel('Latitude ($^\circ$)',fontsize=10)
#plt.ylabel('Depth (m)',fontsize=10)
ax0.set_yticks([])
plt.title('e) Salinity: EP - LP',fontsize=14)
    
cb=plt.colorbar(cs,orientation='vertical')
cb.ax.set_xlabel('psu',fontsize=12)
cb.ax.tick_params(labelsize=12)
cb.ax.xaxis.set_label_position('bottom')
cb.ax.xaxis.labelpad=15

#temperature EP-LP

ax0=fig.add_subplot(spec[5])
EP_cube_a = iris.load_cube(f_ep,'temperature basin')
constraint1 = iris.Constraint(latitude = lambda cell: -80 < cell <  -60)
constraint2 = iris.Constraint(depth_1 = lambda cell: cell < 2200)
EP_cube_b = EP_cube_a.extract(constraint1)
EP_cube = EP_cube_b.extract(constraint2)
       
LP_cube_a = iris.load_cube(f_lp,'temperature basin')
LP_cube_b = LP_cube_a.extract(constraint1)
LP_cube =  LP_cube_b.extract(constraint2)

anom_cube = EP_cube - LP_cube

vals = np.arange(-0.5,0.55,0.05)
cs=iplt.contourf(anom_cube,levels=vals,extend='both',cmap='RdBu_r')
plt.xlabel('Latitude ($^\circ$)',fontsize=10)
#plt.ylabel('Depth (m)',fontsize=10)
ax0.tick_params(labelsize=10)
plt.title('f) Temperature: EP - LP',fontsize=14)
ax0.set_yticks([])

cb=plt.colorbar(cs,orientation='vertical')
cb.ax.set_xlabel('$^\circ$C',fontsize=12)
cb.ax.tick_params(labelsize=12)
cb.ax.xaxis.set_label_position('bottom')
cb.ax.xaxis.labelpad=15


plt.tight_layout()
plt.savefig('Supp_Figure4.png')


# Save the plot to a file
#plt.savefig('seaice/sea-ice_30yr_2_2999.png')
