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
    filestart = '/home/earjcti/um/'
    #filestart = '/uolstore/Research/a/hera1/earjcti/um/'
    tempfile = (filestart + expt + '/database_averages/' + 
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

    #print(expt,EXPTNAMES.get(expt),temp_mean.data-273.15,
    #      NH_6090_mean.data-273.15,SH_6090_mean.data-273.15)

    
    return (temp_mean.data-273.15,
          NH_6090_mean.data-273.15,SH_6090_mean.data-273.15)
                                 






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
#EXPTS=['xqbwn']

NH_amplification=[]
SH_amplification=[]
globalTanom = []
exptnames_used=[]
for expt in EXPTS:
    (globalmeanT,NH_60N90N_T,SH_60S90S_T) = get_temperature_bands(expt)
    
    cntrl = CONTROLNAMES.get(expt)
    (controlT,controlNH_60N90N,controlSH_60S90S) = get_temperature_bands(cntrl)
  
    glob_warm=globalmeanT-controlT
    globalTanom.append(str(np.round(np.abs(glob_warm),1)))
    NH_warm=NH_60N90N_T - controlNH_60N90N
    SH_warm=SH_60S90S_T - controlSH_60S90S

    NH_amplification.append(NH_warm / glob_warm)
    SH_amplification.append(SH_warm / glob_warm)
    exptnames_used.append(EXPTNAMES.get(expt) +  ' - ' + EXPTNAMES.get(cntrl))

    print('j',expt,glob_warm, NH_warm,SH_warm)

sys.exit(0)
# plot as a bar chart with NH and SH next to each other
x=np.arange(len(EXPTS))
width=0.35 # the width of the bars
bar_colors=['tab:blue','tab:orange']
fig,ax=plt.subplots(layout='constrained',figsize=[12,10])
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
ax.set_title('Polar Amplification',fontsize=14)
ax.set_ylabel('Amplification Factor',fontsize=14)
ax.set_xticks(x+0.2,exptnames_used,rotation=90,fontsize=16)
for i,amp in enumerate(NH_amplification):
    print(i,amp,globalTanom[i])
    plt.text(x[i],max(amp, SH_amplification[i],0) + 0.15,globalTanom[i],fontsize=12)  
#ax.legend(loc='upper left',ncol=2,bbox_to_anchor=(-0.05, -0.3))
ax.set_ylim(0.1,10.0)
ax.set_xlim(-0.5,11.75)
ax.set_yscale('symlog',linthresh=2,linscale=0.5)
ax.set_yticks([0.1,0.25,0.5,1,2,4,8,10])
ax.get_yaxis().set_major_formatter(mp.ticker.ScalarFormatter())
ax.set_yticklabels(['0.1','0.25','0.5','1','2','4','8','10'],fontsize=16)

# put some colors dividing the experiments
#plt.axvspan(-0.5,0.75,facecolor='pink',alpha=0.2)
plt.axvspan(0.75,5.75,facecolor='yellow',alpha=0.2)
plt.axvspan(5.75,7.75,facecolor='pink',alpha=0.2)
plt.axvspan(7.75,9.75,facecolor='green',alpha=0.2)
plt.axvspan(9.75,11.75,facecolor='purple',alpha=0.2)
plt.axvline(x=0.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=5.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=7.75,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=9.75,ymin=0.0,ymax=1.0,color='black')
plt.text(2.5, 6.00,'CO$_2$ SENSITIVITY',fontsize=14)
plt.text(6.2, 6.00,'ORBITAL',fontsize=14)
plt.text(8.0, 6.00,'VEGETATION',fontsize=14)
plt.text(10.5, 6.00,'EXT',fontsize=14)
#plt.text(8.0,3.75,'orbital')
#plt.text(10.0,3.75,'vegetaton')
#ax.legend(loc='upper center',fontsize=14)
ax.legend(loc='lower left',ncol=2,bbox_to_anchor=(-0.05, -0.3),fontsize=18)

plt.axhline(y=1.0,color='black')

#plot the data again because of the background
for i,amp in enumerate(NH_amplification):
    rects=ax.bar(x,NH_amplification,width,color=bar_colors[0])
for i,amp in enumerate(SH_amplification):
    rects=ax.bar(x+width,SH_amplification,width,color=bar_colors[1])


plt.savefig('polar_amp_bcs.png')
plt.savefig('polar_amp_bcs.eps')
plt.close()
#plt.show()
           
           
