#python program
#
# all plots for figure 4 of the paper
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import numpy as np
import iris
from iris.analysis.cartography import area_weights
import iris.quickplot as qplt
import iris.plot as iplt
import cartopy.crs as ccrs
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys


##################################################################
# set up figure and subplot
fig = plt.figure(figsize=[10.5,6.5])

spec = gridspec.GridSpec(ncols=5, nrows=3,hspace=0.1,
                         wspace=0.1, width_ratios=[1,0.08, 1,0.08,1],
                         height_ratios=[0.2,1,1.0])


period = {'xqbwd':'LP','xqwbj':'LP$_{490}','xqbwe':'EP$_{400}$',
          'xqbwg':'EP','xqbwc':'PI'}
filestart = '/home/earjcti/um/'
fileend = '_Monthly_Average_#pd_Temperature_3900_4000.nc'

# read in data
file_LP = (filestart + 'xqbwd/database_averages/xqbwd' + fileend)
file_LP490 = (filestart + 'xqbwj/database_averages/xqbwj' + fileend)
file_EP = (filestart + 'xqbwg/database_averages/xqbwg' + fileend)
file_EP400 = (filestart + 'xqbwe/database_averages/xqbwe' + fileend)

field = 'TEMPERATURE AT 1.5M'
# read in monthly data

LPall_cube = iris.load_cube(file_LP,field)
EPall_cube = iris.load_cube(file_EP,field)
LP490all_cube = iris.load_cube(file_LP490,field)
EP400all_cube = iris.load_cube(file_EP400,field)

# get annual mean
LP_cube = iris.util.squeeze(LPall_cube.collapsed('time',iris.analysis.MEAN))
EP_cube = iris.util.squeeze(EPall_cube.collapsed('time',iris.analysis.MEAN))
LP490_cube = iris.util.squeeze(LP490all_cube.collapsed('time',
                                                       iris.analysis.MEAN))
EP400_cube = iris.util.squeeze(EP400all_cube.collapsed('time',
                                                       iris.analysis.MEAN))

# calculate difference and get mean
EP_LP_cube = EP_cube - LP_cube
EP_LP_cube.coord('latitude').guess_bounds()
EP_LP_cube.coord('longitude').guess_bounds()
grid_areas = area_weights(EP_LP_cube)
EP_LP_mean = EP_LP_cube.collapsed(['latitude','longitude'],iris.analysis.MEAN,
                                  weights=grid_areas)

LP490_LP_cube = LP490_cube - LP_cube
LP490_LP_cube.coord('latitude').guess_bounds()
LP490_LP_cube.coord('longitude').guess_bounds()
grid_areas = area_weights(LP490_LP_cube)
LP490_LP_mean = LP490_LP_cube.collapsed(['latitude','longitude'],
                                        iris.analysis.MEAN,
                                        weights=grid_areas)

EP400_LP_cube = EP400_cube - LP_cube
EP400_LP_cube.coord('latitude').guess_bounds()
EP400_LP_cube.coord('longitude').guess_bounds()
grid_areas = area_weights(EP400_LP_cube)
EP400_LP_mean = EP400_LP_cube.collapsed(['latitude','longitude'],
                                        iris.analysis.MEAN,
                                        weights=grid_areas)

EP_LP490_cube = EP_cube - LP490_cube
EP_LP490_cube.coord('latitude').guess_bounds()
EP_LP490_cube.coord('longitude').guess_bounds()
grid_areas = area_weights(EP_LP490_cube)
EP_LP490_mean = EP_LP490_cube.collapsed(['latitude','longitude'],
                                        iris.analysis.MEAN,
                                        weights=grid_areas)

EP_EP400_cube = EP_cube - EP400_cube
EP_EP400_cube.coord('latitude').guess_bounds()
EP_EP400_cube.coord('longitude').guess_bounds()
grid_areas = area_weights(EP_EP400_cube)
EP_EP400_mean = EP_EP400_cube.collapsed(['latitude','longitude'],
                                        iris.analysis.MEAN,
                                        weights=grid_areas)

#EP_LP490_cube = EP_cube - LP490_cube
#EP_EP400_cube = EP_cube - EP400_cube
#EP400_LP_cube = EP400_cube - LP_cube
#LP490_LP_cube = LP490_cube - LP_cube

#######################
# plot temperature
##########################
# EP - LP
ax0=fig.add_subplot(spec[5],projection=ccrs.PlateCarree())

vals = np.arange(-4.0,4.5,0.5)

cs=iplt.contourf(EP_LP_cube,cmap='RdBu_r',levels=vals,extend='both')

titlename = 'a) EP - LP = ' + f'{EP_LP_mean.data:.2f}' + '$^\circ$C'
ax0.set_title(titlename,fontsize=13)
plt.gca().coastlines()


# EP - LP490
ax1=fig.add_subplot(spec[7],projection=ccrs.PlateCarree())

cs=iplt.contourf(EP_LP490_cube,cmap='RdBu_r',levels=vals,extend='both')

titlename = 'b) EP - LP$_{490}$ = ' + f'{EP_LP490_mean.data:.2f}' + '$^\circ$C'
ax1.set_title(titlename,fontsize=13)
plt.gca().coastlines()

# LP490 - LP
ax2=fig.add_subplot(spec[9],projection=ccrs.PlateCarree())

cs=iplt.contourf(LP490_LP_cube,cmap='RdBu_r',levels=vals,extend='both')

titlename = 'c) LP$_{490}$ - LP = ' + f'{LP490_LP_mean.data:.2f}' + '$^\circ$C'
ax2.set_title(titlename,fontsize=13)
plt.gca().coastlines()

# EP400 - LP
ax4=fig.add_subplot(spec[12],projection=ccrs.PlateCarree())

cs=iplt.contourf(EP_LP_cube,cmap='RdBu_r',levels=vals,extend='both')

titlename = 'd) EP$_{400}$ - LP = ' + f'{EP400_LP_mean.data:.2f}' + '$^\circ$C'
ax4.set_title(titlename,fontsize=13)
plt.gca().coastlines()


# EP - EP400
ax5=fig.add_subplot(spec[14],projection=ccrs.PlateCarree())

cs=iplt.contourf(EP_EP400_cube,cmap='RdBu_r',levels=vals,extend='both')

titlename = 'e) EP - EP$_{400}$ = ' + f'{EP_EP400_mean.data:.2f}' + '$^\circ$C'
ax5.set_title(titlename,fontsize=13)
plt.gca().coastlines()


# now add colorbar
cax = fig.add_subplot(spec[1:3,0])
pos = cax.get_position()
cax.set_position([pos.x0,pos.y0+0.15,pos.width,pos.height*0.1])
cbar = fig.colorbar(cs, cax=cax,orientation='horizontal')
cbar.ax.tick_params(labelsize=14)
cbar.set_label('°C',fontsize=14)


# now draw black box around b and d
axt = fig.add_subplot(spec[2])
axt.set_axis_off()
axt.text(0.5, 1.0,'Effects of CAS opening',fontsize=14,va='top',
           fontweight='bold',ha='center')

positions = [ax1.get_position(),ax4.get_position(),axt.get_position()]

left   = min(p.x0 for p in positions)
right  = max(p.x1 for p in positions)
bottom = min(p.y0 for p in positions)
top    = max(p.y1 for p in positions)
pad = 0.01
title_space = 0.06

rect = Rectangle(
    (left-pad, bottom-pad),
    (right-left)+2*pad,
    (top-bottom)+2*pad + title_space,
    fill=False,
    edgecolor='black',
    linewidth=2,
    transform=fig.transFigure,
    zorder=1000
)

fig.add_artist(rect)

# now draw black box around a and title
axt3 = fig.add_subplot(spec[0])
axt3.set_axis_off()
axt3.text(0.5, 1.0,'Total EP - LP warming',fontsize=14,va='top',
          fontweight='bold',ha='center')
positions = [ax0.get_position(),axt3.get_position()]

left   = min(p.x0 for p in positions)
right  = max(p.x1 for p in positions)
bottom = min(p.y0 for p in positions)
top    = max(p.y1 for p in positions)
pad = 0.01
title_space = 0.06

rect = Rectangle(
    (left-pad, bottom-pad),
    (right-left)+2*pad,
    (top-bottom)+2*pad + title_space,
    fill=False,
    edgecolor='black',
    linewidth=2,
    transform=fig.transFigure,
    zorder=1000)
fig.add_artist(rect)

# draw black box around c and e
axt2 = fig.add_subplot(spec[4])
axt2.set_axis_off()
axt2.text(0.5, 1.0,'Effects of CO$_2$ increase',fontsize=14,va='top',
           fontweight='bold',ha='center')
positions = [ax2.get_position(),ax5.get_position(),axt2.get_position()]

left   = min(p.x0 for p in positions)
right  = max(p.x1 for p in positions)
bottom = min(p.y0 for p in positions)
top    = max(p.y1 for p in positions)
print(left,right,bottom,top)
pad = 0.01
title_space = 0.06

rect = Rectangle(
    (left-pad, bottom-pad),
    (right-left)+2*pad,
    (top-bottom)+2*pad + title_space,
    fill=False,
    edgecolor='black',
    linewidth=2,
    transform=fig.transFigure,
    zorder=1000
)

fig.add_artist(rect)


#plt.tight_layout()
#plt.show()
#sys.exit(0)
plt.savefig('Figure1.png')
plt.close()

