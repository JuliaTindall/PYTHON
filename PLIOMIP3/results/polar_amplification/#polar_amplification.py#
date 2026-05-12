"""
#NAME
#    polar_amplification.py
#PURPOSE 
#
#  This program will show the polar amplification for each experiment
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
import cartopy.crs as ccrs
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import iris.quickplot as qplt
import iris.plot as iplt
from netCDF4 import Dataset, MFDataset
import sys
#from netCDF4 import Dataset, MFDataset

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def get_temperature_bands(expt):
    """
    gets the polar amplification for the NH and the SH
    """
    tempfile = ('/nfs/hera1/earjcti/um/' + expt + '/database_averages/' + 
                expt + '_Annual_Average_#pd_Temperature_3900_4000.nc')
    temp_cube = iris.load_cube(tempfile)
    temp_cube=iris.util.squeeze(temp_cube)
    temp_cube.coord('latitude').guess_bounds()
    temp_cube.coord('longitude').guess_bounds()
    grid_areas=iris.analysis.cartography.area_weights(temp_cube)

    # global mean
    temp_mean = temp_cube.collapsed(['longitude','latitude'],
                                    iris.analysis.MEAN,
                                    weights=grid_areas)
  
    # 60-90N
    grid_areas_NH = np.copy(grid_areas)
    grid_areas_SH = np.copy(grid_areas)
    for j,lat in enumerate(temp_cube.coord('latitude').points):
        if lat < 60:
            grid_areas_NH[j,:]=0.0
        if lat > -60:
            grid_areas_SH[j,:]=0.0
    NH_6090_mean=temp_cube.collapsed(['longitude','latitude'],
                                    iris.analysis.MEAN,
                                    weights=grid_areas_NH)
    SH_6090_mean=temp_cube.collapsed(['longitude','latitude'],
                                    iris.analysis.MEAN,
                                    weights=grid_areas_SH)

    print(expt,EXPTNAMES.get(expt),temp_mean.data-273.15,
          NH_6090_mean.data-273.15,SH_6090_mean.data-273.15)

    
    return (temp_mean.data-273.15,
          NH_6090_mean.data-273.15,SH_6090_mean.data-273.15)
                                 






CONTROL='xqbwc'
EXPTS=['xqbwd','xqbwg','xqbwr','xqbwl','xqbwm',
       'xqbwi','xqbwj','xqbwk','xqbwn','xqbwo','xqbwt','xqbws']
EXPTNAMES={'xqbwc':'PI', 'xqbwd':'LP','xqbwg':'EP',
           'xqbwr':'LP_alt','xqbwl':'PI400','xqbwm':'PI560',
           'xqbwi':'LP280','xqbwj':'LP490','xqbwk':'LP560',
           'xqbwn':'LP_highNH_orb','xqbwo':'LP_lowNH_orb',
           'xqbwt':'PI_dyn-veg','xqbws':'LP_dyn-veg'}
#EXPTS=['xqbwd','xqbwg','xqbwr']
(controlT, controlNH_60N90N, controlSH_60S90S)=get_temperature_bands(CONTROL)

NH_amplification=[]
SH_amplification=[]
exptnames_used=[]
for expt in EXPTS:
    (globalmeanT,NH_60N90N_T,SH_60S90S_T) = get_temperature_bands(expt)
    glob_warm=globalmeanT-controlT
    NH_warm=NH_60N90N_T - controlNH_60N90N
    SH_warm=SH_60S90S_T - controlSH_60S90S

    NH_amplification.append(NH_warm / glob_warm)
    SH_amplification.append(SH_warm / glob_warm)
    exptnames_used.append(EXPTNAMES.get(expt))

# plot as a bar chart with NH and SH next to each other
x=np.arange(len(EXPTS))
width=0.35 # the width of the bars
bar_colors=['tab:blue','tab:orange']
fig,ax=plt.subplots(layout='constrained')
#plot NH
for i,amp in enumerate(NH_amplification):
    if i == 0:
        rects=ax.bar(x,NH_amplification,width,label='Northern Hemisphere',
                    color=bar_colors[0])
    else:
        rects=ax.bar(x,NH_amplification,width,color=bar_colors[0])
#    ax.bar_label(rects,EXPTS[i],padding=3)
for i,amp in enumerate(SH_amplification):
    if i == 0:
        rects=ax.bar(x+width,SH_amplification,width,label='Southern Hemisphere',
                     color=bar_colors[1])
    else:
        rects=ax.bar(x+width,SH_amplification,width,color=bar_colors[1])
# labels
ax.set_title('Polar Amplification')
ax.set_ylabel('Amplification Factor')
ax.set_xticks(x+0.2,exptnames_used,rotation=90)
ax.legend(loc='lower left',ncol=2,bbox_to_anchor=(-0.05, -0.3))
ax.set_ylim(1.0,4.0)
ax.set_xlim(-0.5,11.75)

# put some colors dividing the experiments
plt.axvspan(-0.5,2.75,facecolor='pink',alpha=0.2)
plt.axvspan(2.75,7.75,facecolor='yellow',alpha=0.2)
plt.axvspan(7.75,9.75,facecolor='purple',alpha=0.2)
plt.axvspan(9.75,11.75,facecolor='green',alpha=0.2)
plt.axvline(x=2.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=7.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=9.75,ymin=0.0,ymax=1.0,color='black')
plt.text(-0.25,3.75,'core and extension')
plt.text(4.0, 3.75,'CO2 sensitivity')
plt.text(8.0,3.75,'orbital')
plt.text(10.0,3.75,'vegetaton')

#plot the data again because of the background
for i,amp in enumerate(NH_amplification):
    rects=ax.bar(x,NH_amplification,width,color=bar_colors[0])
for i,amp in enumerate(SH_amplification):
    rects=ax.bar(x+width,SH_amplification,width,color=bar_colors[1])

plt.show()
           
           
