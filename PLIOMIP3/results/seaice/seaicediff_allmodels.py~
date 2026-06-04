"""
#NAME
#    seaice_allmodels.py
#PURPOSE 
#
#  This program will show the seaice area in each model
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


def get_seaice(expt):
    """
    gets the seaice for the NH and the SH
    """
    useann = {'xqbwn':'y','xqbwo':'y','xqbwd':'y'}
    filestart = '/home/earjcti/um/'
    if useann.get(expt,'n') == 'y':
        fileend = '_Annual_Average_#pg_seaice_3900_4000.nc'
    else:
        fileend = '_Monthly_Average_#pd_SeaIceConc_3900_4000.nc'
    
    #filestart = '/uolstore/Research/a/hera1/earjcti/um'

    tempfile = (filestart + expt + '/database_averages/' + expt + fileend)
    temp_cube = iris.load_cube(tempfile)
    temp_cube=iris.util.squeeze(temp_cube)
    temp_cube.data = np.where(temp_cube.data < 0.0,  0.0, temp_cube.data)

    # get north and south
    SH_constraint = iris.Constraint(latitude = lambda cell: cell < 0)
    NH_constraint = iris.Constraint(latitude = lambda cell: cell > 0)
    SH_cube = temp_cube.extract(SH_constraint)
    NH_cube = temp_cube.extract(NH_constraint)

    if useann.get(expt,'n') == 'y':
        NH_ann_cube = NH_cube.copy()
        SH_ann_cube = SH_cube.copy()
    else:
        NH_ann_cube = NH_cube.collapsed(['time'],iris.analysis.MEAN)
        SH_ann_cube = SH_cube.collapsed(['time'],iris.analysis.MEAN)
  

    # get N/S mean
    NH_ann_cube.coord('latitude').guess_bounds()
    NH_ann_cube.coord('longitude').guess_bounds()
    NH_grid_areas=iris.analysis.cartography.area_weights(NH_ann_cube)
    SH_ann_cube.coord('latitude').guess_bounds()
    SH_ann_cube.coord('longitude').guess_bounds()
    SH_grid_areas=iris.analysis.cartography.area_weights(SH_ann_cube)
    # set land to zero
   
   
    print(NH_ann_cube)
    print(NH_grid_areas.shape)
    # NH / SH mean
    NH_total_cube = NH_ann_cube.collapsed(['longitude','latitude'],
                                    iris.analysis.SUM,
                                    weights=NH_grid_areas)
    SH_total_cube = SH_ann_cube.collapsed(['longitude','latitude','time'],
                                    iris.analysis.SUM,
                                    weights=SH_grid_areas)
    # divide by 1E6 to convert to km2
    NH_ann_ice = NH_total_cube.data / 1.0E6
    SH_ann_ice = SH_total_cube.data / 1.0E6
  
    #print(NH_ann_ice)
    #print(SH_ann_ice)
    
    return (NH_ann_ice,SH_ann_ice)
                                 






EXPTS=['xqbwc','xqbwl','xqbwm','xqbwi','xqbwd','xqbwj','xqbwk','xqbwr','xqbwg',
       'xqbwn','xqbwo',
       'xqbwt','xqbws']
EXPTNAMES={'xqbwc':'PI', 'xqbwd':'LP','xqbwg':'EP',
           'xqbwr':'LP_alt','xqbwl':'PI400','xqbwm':'PI560',
           'xqbwi':'LP280','xqbwj':'LP490','xqbwk':'LP560',
           'xqbwn':'LP_highNH_orb','xqbwo':'LP_lowNH_orb',
           'xqbwt':'PI_dyn-veg','xqbws':'LP_dyn-veg'}

NH_ice_all=[]
SH_ice_all=[]
exptnames_used=[]
for expt in EXPTS:
    (NH_ann_ice, SH_ann_ice) = get_seaice(expt)
    exptnames_used.append(EXPTNAMES.get(expt))
    NH_ice_all.append(NH_ann_ice)
    SH_ice_all.append(SH_ann_ice)

# plot as a bar chart with NH and SH next to each other
x=np.arange(len(EXPTS))
width=0.35 # the width of the bars
bar_colors=['tab:blue','tab:orange']
fig,ax=plt.subplots(layout='constrained',figsize=[12.,10.])
#plot NH
rects=ax.bar(x,NH_ice_all,width,label='Northern Hemisphere',
                    color=bar_colors[0])
rects=ax.bar(x+width,SH_ice_all,width,label='Southern Hemisphere',
                     color=bar_colors[1])

# labels
ax.set_title('Sea ice',fontsize=16)
ax.set_ylabel('Area',fontsize=16)
ax.set_xticks(x+0.2,exptnames_used,rotation=90,fontsize=16)
#ax.set_yticklabels(['1','1.5','2','2.5','3','3.5','4'],fontsize=16)
ax.legend(loc='lower left',ncol=2,bbox_to_anchor=(-0.05, -0.3),fontsize=16)
#ax.set_ylim(1.0,4.0)
ax.set_xlim(-0.5,12.75)
# put some colors dividing the experiments
plt.axvspan(-0.5,2.75,facecolor='pink',alpha=0.2)
plt.axvspan(2.75,8.75,facecolor='yellow',alpha=0.2)
plt.axvspan(8.75,10.75,facecolor='purple',alpha=0.2)
plt.axvspan(10.75,12.75,facecolor='green',alpha=0.2)
plt.axvline(x=2.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=8.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=6.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=10.75,ymin=0.0,ymax=1.0,color='black')
plt.text(1.,1.0E7,'PI',fontsize=16)
plt.text(5.0, 1.0E7,'LP',fontsize=16)
plt.text(9.5,1.0E7,'orbital',fontsize=16)
plt.text(11.0,1.0E7,'vegetaton',fontsize=16)

#plot the data again because of the background
for i,amp in enumerate(NH_ice_all):
    rects=ax.bar(x,NH_ice_all,width,color=bar_colors[0])
for i,amp in enumerate(SH_ice_all):
    rects=ax.bar(x+width,SH_ice_all,width,color=bar_colors[1])

plt.savefig('seaice_allexpts.png')
plt.show()
           
           
