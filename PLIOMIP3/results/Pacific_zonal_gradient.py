"""
#NAME
#    pacific_zonal_gradient.py
#PURPOSE 
#
#  This program will plot the gradient across the pacific for each experiment
#  
"""

# Import necessary libraries

import os
import numpy as np
import math
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.quickplot as qplt
import iris.plot as iplt
#from netCDF4 import Dataset, MFDataset
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def get_Pacific_temps(expt):
    """
    gets the mean 20-20 temperature by longitude across the Pacific
    """
    tempfile = ('/nfs/hera1/earjcti/um/' + expt + '/database_averages/' + 
                expt + '_Annual_Average_#pf_SST_3900_4000.nc')
    print(tempfile)
    temp_cube = iris.load_cube(tempfile)
    temp_cube=iris.util.squeeze(temp_cube)
    temp_cube.coord('latitude').guess_bounds()
    temp_cube.coord('longitude').guess_bounds()
    grid_areas=iris.analysis.cartography.area_weights(temp_cube)
    for j,lat in enumerate(temp_cube.coord('latitude').points):
        if (lat > TROPNORTH or lat < TROPSOUTH):
            grid_areas[j,:]=0.0

   
    # 20N-20S mean
    temp_mean_cube = temp_cube.collapsed(['latitude'],
                                         iris.analysis.MEAN,
                                         weights=grid_areas)
    pac_constraint = iris.Constraint(longitude=lambda cell:
                                     150.0 < cell < 260.0)                     
    trop_pac_mean = temp_mean_cube.extract(pac_constraint)
    trop_pac_mean.data = trop_pac_mean.data
    trop_pac_mean.units='degC'
    
    return (trop_pac_mean)
                                 





###############################################################
CONTROL='xqbwc'
TROPNORTH = 0.0
TROPSOUTH = -5.0  # the tropical range 5=5N-5S
EXPTS=['xqbwc','xqbwd','xqbwg','xqbwr','xqbwl','xqbwm',
       'xqbwi','xqbwj','xqbwk','xqbwn','xqbwo','xqbwt','xqbws']
EXPTNAMES={'xqbwc':'PI', 'xqbwd':'LP','xqbwg':'EP',
           'xqbwr':'LP_alt','xqbwl':'PI400','xqbwm':'PI560',
           'xqbwi':'LP280','xqbwj':'LP490','xqbwk':'LP560',
           'xqbwn':'LP_highNH_orb','xqbwo':'LP_lowNH_orb',
           'xqbwt':'PI_dyn-veg','xqbws':'LP_dyn-veg'}
#EXPTS=['xqbwc','xqbwd','xqbwg']

pac_temp_cubelist = CubeList([])
gradient=[]
exptnames_used=[]
for expt in EXPTS:
    pacific_temperature_cube = get_Pacific_temps(expt)
    pac_temp_cubelist.append(pacific_temperature_cube)
    gradient.append(np.max(pacific_temperature_cube.data) - 
                    np.min(pacific_temperature_cube.data))
    exptnames_used.append(EXPTNAMES.get(expt))


# plot all the pacific temperatures
for i,cube in enumerate(pac_temp_cubelist):
   qplt.plot(cube,label=EXPTNAMES.get(EXPTS[i]))
plt.subplots_adjust(bottom=0.3)
plt.legend(loc='lower left',ncol=4,bbox_to_anchor=(-0.05, -0.4),
           borderaxespad=0)
print(gradient)
plt.savefig('plots/Pacific_zonal_gradient1_'+str(TROPNORTH) + 'N-' + 
            str(np.abs(TROPSOUTH))+'S.eps') 
plt.close()

# plot as a bar chart with NH and SH next to each other
x=np.arange(len(EXPTS))
width=0.8 # the width of the bars
bar_colors=['tab:blue']
fig,ax=plt.subplots(layout='constrained')
#plot NH
for i,amp in enumerate(gradient):
    if i == 0:
        rects=ax.bar(x,gradient,width,
                    color=bar_colors[0])
    else:
        rects=ax.bar(x,gradient,width,color=bar_colors[0])
#    ax.bar_label(rects,EXPTS[i],padding=3)

# labels
ax.set_title('Tropical pacific gradient (max value - min value)')
ax.set_ylabel('Temperature difference (degC)')
ax.set_xticks(x,exptnames_used,rotation=90)
ax.set_ylim(0.0,4.0)
ax.set_xlim(-0.5,12.5)

# put some colors dividing the experiments
plt.axvspan(-0.5,3.5,facecolor=(1.0,0.9,0.9))
plt.axvspan(3.5,8.5,facecolor=(1.0,1.0,0.9))
plt.axvspan(8.5,10.5,facecolor=(0.9,0.9,0.9))
plt.axvspan(10.5,12.5,facecolor=(0.9,1.0,0.9))
plt.axvline(x=3.5,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=8.5,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=10.5,ymin=0.0,ymax=1.0,color='black')
plt.text(-0.25,3.75,'core and extension')
plt.text(5.0, 3.75,'CO2 sensitivity')
plt.text(9.0,3.75,'orbital')
plt.text(10.75,3.75,'vegetaton')

#plot the data again because of the background
for i,amp in enumerate(gradient):
    rects=ax.bar(x,gradient,width,color=bar_colors[0])

plt.savefig('plots/Pacific_zonal_gradient2_'+str(TROPNORTH) + 'N-' + 
            str(np.abs(TROPSOUTH))+'S.eps') 
plt.close()
           
           
