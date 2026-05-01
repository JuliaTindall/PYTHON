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

spec = gridspec.GridSpec(ncols=2, nrows=3,
                         width_ratios=[1, 1], wspace=0.5,
                         hspace=0.5, height_ratios=[1, 1,0.3])


period = {'xqbwd':'LP','xqwbj':'LP$_{490}','xqbwe':'EP$_{400}$',
          'xqbwg':'EP','xqbwc':'PI'}
filestart = '/home/earjcti/um/'

# read in data
fileg_Atl = (filestart +
             'xqbwg/basin_diagnostics/mean_xqbwg_Atlantic3900_3999.nc')
filej_Atl = (filestart +
             'xqbwj/basin_diagnostics/mean_xqbwj_Atlantic3900_3999.nc')
fileg_Pac = (filestart +
             'xqbwg/basin_diagnostics/mean_xqbwg_Pacific3900_3999.nc')
filej_Pac = (filestart +
             'xqbwj/basin_diagnostics/mean_xqbwj_Pacific3900_3999.nc')


#temperature xqbwg-xqbwj Atlantic
ax0=fig.add_subplot(spec[0])
xqbwg_atltemp_cube = iris.load_cube(fileg_Atl,'temperature basin')
xqbwj_atltemp_cube = iris.load_cube(filej_Atl,'temperature basin')

Tanom_cube = xqbwg_atltemp_cube - xqbwj_atltemp_cube

cs=iplt.contourf(Tanom_cube,cmap='RdBu_r',
              levels=np.arange(-2.0,2.2,0.2),extend='both')

ax0.set_title('a) Atlantic temperature anomaly: EP-LP$_{400}$',fontsize=15)
plt.ylabel('depth (m)',fontsize=14)
#plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax0.set_xticks([-60,-30,0,30,60])


#colorbar
cax_holder=fig.add_subplot(spec[4])
cax_holder.axis("off")
cax = inset_axes(
    cax_holder,
    width="100%",   # fraction of the subplot width
    height="30%",  # fraction of the subplot height
    loc="upper center"
)
cbar=plt.colorbar(cs, cax=cax, orientation="horizontal")
cbar.set_label("$^\circ$ C", fontsize=15)
cbar.ax.tick_params(labelsize=14,labelrotation=45)

cax_holder.spines[:].set_visible(False)
cax_holder.tick_params(axis='both', which='both', bottom=False, top=False,
                  right=False, left=False,labelbottom=False,
                  labeltop=False,labelleft=False,labelright=False)
                      




# now do atlantic salinity anomaly
ax1=fig.add_subplot(spec[1])
xqbwg_atl_cube = iris.load_cube(fileg_Atl,'salinity basin')
xqbwj_atl_cube = iris.load_cube(filej_Atl,'salinity basin')

Sanom_cube = xqbwg_atl_cube - xqbwj_atl_cube

cs2=iplt.contourf(Sanom_cube,cmap='RdBu_r',
              levels=np.arange(-0.5,0.55,0.05),extend='both')

ax1.set_title('b) Atlantic salinity anomaly: EP-LP$_{400}$',fontsize=15)
plt.ylabel('depth (m)',fontsize=14)
#plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax1.set_xticks([-60,-30,0,30,60])

#colorbar
cax_holder=fig.add_subplot(spec[5])
cax_holder.axis("off")
cax = inset_axes(
    cax_holder,
    width="100%",   # fraction of the subplot width
    height="30%",  # fraction of the subplot height
    loc="upper center"
)
cbar=plt.colorbar(cs2, cax=cax, orientation="horizontal")
cbar.set_label("psu", fontsize=15)
cbar.ax.tick_params(labelsize=14,labelrotation=45)

cax_holder.spines[:].set_visible(False)
cax_holder.tick_params(axis='both', which='both', bottom=False, top=False,
                  right=False, left=False,labelbottom=False,
                  labeltop=False,labelleft=False,labelright=False)
                      


#temperature xqbwg-xqbwj Pacific
ax2=fig.add_subplot(spec[2])
xqbwg_pactemp_cube = iris.load_cube(fileg_Pac,'temperature basin')
xqbwj_pactemp_cube = iris.load_cube(filej_Pac,'temperature basin')

Tanom_cube = xqbwg_pactemp_cube - xqbwj_pactemp_cube

iplt.contourf(Tanom_cube,cmap='RdBu_r',
              levels=np.arange(-2.0,2.2,0.2),extend='both')

ax2.set_title('c) Pacific temperature anomaly: EP-LP$_{400}$',fontsize=15)
plt.ylabel('depth (m)',fontsize=14)
plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax2.set_xticks([-60,-30,0,30,60])

# now do pacific salinity anomaly
ax3=fig.add_subplot(spec[3])
xqbwg_pac_cube = iris.load_cube(fileg_Pac,'salinity basin')
xqbwj_pac_cube = iris.load_cube(filej_Pac,'salinity basin')

Sanom_cube = xqbwg_pac_cube - xqbwj_pac_cube

iplt.contourf(Sanom_cube,cmap='RdBu_r',
              levels=np.arange(-0.5,0.55,0.05),extend='both')

ax3.set_title('d) Pacific salinity anomaly: EP-LP$_{400}$',fontsize=15)
plt.ylabel('depth (m)',fontsize=14)
plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax3.set_xticks([-60,-30,0,30,60])



plt.tight_layout()
plt.savefig('Figure4.png')


# Save the plot to a file
#plt.savefig('seaice/sea-ice_30yr_2_2999.png')
