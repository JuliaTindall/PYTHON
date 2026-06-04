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
                                 



EXPTNAMES={'xqbwc':'PI', 'xqbwd':'LP','xqbwg':'EP',
           'xqbwr':'LP_alt','xqbwl':'PI$_{400}$','xqbwm':'PI$_{560}$',
           'xqbwi':'LP$_{280}$','xqbwj':'LP$_{490}$','xqbwk':'LP$_{560}$',
           'xqbwn':'LP$_{highNHorb}$','xqbwo':'LP$_{lowNHorb}$',
           'xqbwt':'PI$_{dynveg}$','xqbws':'LP$_{dynveg}$'}

CONTROLNAMES = {'xqbwi':'xqbwc', # LP280, PI
                'xqbwl':'xqbwc', # PI400, PI
                'xqbwm':'xqbwc', # PI560, PI
                'xqbwd':'xqbwi', # LP, LP280
                'xqbwj':'xqbwi', # LP490, LP280
                'xqbwk':'xqbwi', # lP560, LP280
                'xqbwn':'xqbwd', # LP_high orbit, LP
                'xqbwo':'xqbwd', # LP_low_orbit, LP
                'xqbwt':'xqbwc', # PI dyn veg, PI
                'xqbws':'xqbwd', # Plio dyn veg, PI
                'xqbwg':'xqbwj', # EP - LP490
                'xqbwr':'xqbwd'}  # LP_alt - LP
                
EXPTS=['xqbwi','xqbwl','xqbwm','xqbwd','xqbwj','xqbwk','xqbwn',
       'xqbwo','xqbwt','xqbws','xqbwg','xqbwr']





NH_ice_all=[]
SH_ice_all=[]
exptnames_used=[]
for expt in EXPTS:
    (NH_ann_ice, SH_ann_ice) = get_seaice(expt)
    cntrl = CONTROLNAMES.get(expt)
    (NH_ann_ice_cntrl, SH_ann_ice_cntrl) = get_seaice(cntrl)

    
    NH_ice_all.append((NH_ann_ice - NH_ann_ice_cntrl)/ NH_ann_ice)
    SH_ice_all.append((SH_ann_ice - SH_ann_ice_cntrl) / SH_ann_ice)

    exptnames_used.append(EXPTNAMES.get(expt) +  ' - ' + EXPTNAMES.get(cntrl))

# plot as a bar chart with NH and SH next to each other
x=np.arange(len(EXPTS))
width=0.35 # the width of the bars
bar_colors=['tab:blue','tab:orange']
fig,ax=plt.subplots(figsize=[12.,10.])
#plot NH and SH
rects=ax.bar(x,NH_ice_all,width,label='Northern Hemisphere',
                    color=bar_colors[0])
rects=ax.bar(x+width,SH_ice_all,width,label='Southern Hemisphere',
                     color=bar_colors[1])

# labels
ax.set_title('Sea ice',fontsize=16)
ax.set_ylabel('Fraction lost',fontsize=16)
ax.set_xticks(x+0.2,exptnames_used,rotation=90,fontsize=16)
#ax.set_yticklabels(['1','1.5','2','2.5','3','3.5','4'],fontsize=16)
ax.legend(loc='lower left',ncol=2,bbox_to_anchor=(-0.05, -0.4),fontsize=16)
#ax.set_ylim(1.0,4.0)
ax.set_xlim(-0.5,12.75)
# put some colors dividing the experiments
plt.axvspan(0.75,5.75,facecolor='yellow',alpha=0.2)
plt.axvspan(5.75,7.75,facecolor='pink',alpha=0.2)
plt.axvspan(7.75,9.75,facecolor='green',alpha=0.2)
plt.axvspan(9.75,11.75,facecolor='purple',alpha=0.2)
plt.axvline(x=0.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=5.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=7.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=9.75,ymin=0.0,ymax=1.0,color='black')
plt.text(1.,1.0E7,'PI',fontsize=16)
plt.text(5.0, 1.0E7,'LP',fontsize=16)
plt.text(9.5,1.0E7,'orbital',fontsize=16)
plt.text(11.0,1.0E7,'vegetaton',fontsize=16)


#plot the data again because of the background
for i,amp in enumerate(NH_ice_all):
    rects=ax.bar(x,NH_ice_all,width,color=bar_colors[0])
for i,amp in enumerate(SH_ice_all):
    rects=ax.bar(x+width,SH_ice_all,width,color=bar_colors[1])

plt.subplots_adjust(bottom=0.3)
plt.show()
plt.savefig('seaice_allexpts.png')
plt.close()           
           
