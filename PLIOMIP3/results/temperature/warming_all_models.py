"""
#NAME
#   warming_all_models.py  
#PURPOSE 
#
#   This will show the mean / winter / summer warming for all models
#
#  also has the option of splitting into land and ocean
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
from netCDF4 import Dataset, MFDataset
import sys
#from netCDF4 import Dataset, MFDataset

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")


def get_mean_val(expt,cube):
    """
    gets the mean of the cube over the globe, land or ocean as appropriate
    """
    direct={'xqbwc':'preind2', 'xqbwd':'P4_enh','xqbwg':'EP',
           'xqbwr':'LP_MIN_LSM_Change','xqbwl':'preind2','xqbwm':'preind2',
           'xqbwi':'P4_enh','xqbwj':'P4_enh','xqbwk':'P4_enh',
           'xqbwn':'P4_enh','xqbwo':'P4_enh',
           'xqbwt':'preind2','xqbws':'P4_enh'}

    filestart={'xqbwc':'qrparm.', 'xqbwd':'P4_enh_qrparm.','xqbwg':'EP_',
           'xqbwr':'LPLSMmin_','xqbwl':'qrparm.','xqbwm':'qrparm.',
           'xqbwi':'P4_enh_qrparm.','xqbwj':'P4_enh_qrparm.',
           'xqbwk':'P4_enh_qrparm.',
           'xqbwn':'P4_enh_qrparm.','xqbwo':'P4_enh_qrparm.',
           'xqbwt':'qrparm.','xqbws':'P4_enh_qrparm.'}


    cube=iris.util.squeeze(cube)
    if cube.coord('latitude').has_bounds():
        pass
    else:
        cube.coord('latitude').guess_bounds()
    if cube.coord('longitude').has_bounds():
        pass
    else:
        cube.coord('longitude').guess_bounds()

    if LAND_OCN_IND == 'Globe':
        grid_areas=iris.analysis.cartography.area_weights(cube)
    if LAND_OCN_IND == 'NH':
        grid_temp=iris.analysis.cartography.area_weights(cube)
        grid_areas=np.copy(grid_temp)
        for j,lat in enumerate(cube.coord('latitude').points):
            if lat < 0.0:
                grid_areas[j,:]=0.0
    if LAND_OCN_IND == 'SH':
        grid_temp=iris.analysis.cartography.area_weights(cube)
        grid_areas=np.copy(grid_temp)
        for j,lat in enumerate(cube.coord('latitude').points):
            if lat > 0.0:
                grid_areas[j,:]=0.0
  
    filename = ('/uolstore/Research/a/hera1/earjcti/ancil/'+direct.get(expt)+'/' + 
                filestart.get(expt) + 'mask.nc')
    if LAND_OCN_IND == 'Land':
        lsm_cube = iris.load_cube(filename,'LAND MASK (LOGICAL: LAND=TRUE)')
        lsm_cube = iris.util.squeeze(lsm_cube)
        grid_areas=iris.analysis.cartography.area_weights(cube)
        grid_areas=grid_areas * lsm_cube.data

    if LAND_OCN_IND == 'Ocean':
        lsm_cube = iris.load_cube(filename,'LAND MASK (LOGICAL: LAND=TRUE)')
        lsm_cube = iris.util.squeeze(lsm_cube)
        grid_areas=iris.analysis.cartography.area_weights(cube)
        grid_areas=grid_areas * (lsm_cube.data * (-1.0) + 1.0)

    # global mean
    cube_mean = cube.collapsed(['longitude','latitude'],
                                    iris.analysis.MEAN,
                                    weights=grid_areas)
 
    return cube_mean

def get_globmean_temperature(expt):
    """
    gets the globalmean temperature for the annual average 
    the JJA average and the DJF average
    """
    tempfile = ('/uolstore/Research/a/hera1/earjcti/um/' + expt + '/database_averages/' + 
                expt + '_Annual_Average_#pd_Temperature_3900_4000.nc')
    ann_temp_cube = iris.load_cube(tempfile)
    ann_temp_cube=iris.util.squeeze(ann_temp_cube)

  
    # get the values for the seasons
    tempfile = ('/uolstore/Research/a/hera1/earjcti/um/' + expt + '/database_averages/' + 
                expt + '_Monthly_Average_#pd_Temperature_3900_4000.nc')
    temp_cube = iris.load_cube(tempfile)
    temp_cube=iris.util.squeeze(temp_cube)
    temp_alt_cube = temp_cube.collapsed(['time'],iris.analysis.MEAN)
   
    djf_cubelist = CubeList([])
    djf_cubelist.append(temp_cube[0,:,:,:])
    djf_cubelist.append(temp_cube[1,:,:,:])
    djf_cubelist.append(temp_cube[11,:,:,:])
    djf_3cube = djf_cubelist.merge_cube()
    djf_cube = djf_3cube.collapsed(['time'],iris.analysis.MEAN)
   
    jja_cubelist = CubeList([])
    jja_cubelist.append(temp_cube[5,:,:,:])
    jja_cubelist.append(temp_cube[6,:,:,:])
    jja_cubelist.append(temp_cube[7,:,:,:])
    jja_3cube = jja_cubelist.merge_cube()
    jja_cube = jja_3cube.collapsed(['time'],iris.analysis.MEAN)
   
  
    
    # get the mean values (this will return global / land /ocn mean as reqd
    ann_temp_mean = get_mean_val(expt,ann_temp_cube)
    ann_temp_mean_alt = get_mean_val(expt,temp_alt_cube)
    djf_temp_mean = get_mean_val(expt,djf_cube)
    jja_temp_mean = get_mean_val(expt,jja_cube)

   
    return (ann_temp_mean.data-273.15,
          jja_temp_mean.data-273.15,djf_temp_mean.data-273.15)
                                 

############################################################




CONTROL='xqbwc'
EXPTS=['xqbwd','xqbwg','xqbwr','xqbwl','xqbwm',
       'xqbwi','xqbwj','xqbwk',
       'xqbwn','xqbwo',
       'xqbwt','xqbws']
EXPTNAMES={'xqbwc':'PI', 'xqbwd':'LP','xqbwg':'EP',
           'xqbwr':'LP_alt','xqbwl':'PI400','xqbwm':'PI560',
           'xqbwi':'LP280','xqbwj':'LP490','xqbwk':'LP560',
           'xqbwn':'LP_highNH_orb','xqbwo':'LP_lowNH_orb',
           'xqbwt':'PI_dyn-veg','xqbws':'LP_dyn-veg'}

#CONTROL = 'xqbwd'
#EXPTS = ['xqbwj','xqbwg']
LAND_OCN_IND='NH'  # options are Land; Ocean, Globe, NH, SH
#EXPTS=['xqbwn']
(controlTann, controlTjja, controlTdjf)=get_globmean_temperature(CONTROL)

mean_ann_warming=[]
mean_jja_warming=[]
mean_djf_warming=[]
exptnames_used=[]

for expt in EXPTS:
    (exptTann,exptTjja,exptTdjf) = get_globmean_temperature(expt)
    mean_ann_warming.append(exptTann - controlTann)
    mean_jja_warming.append(exptTjja - controlTjja)
    mean_djf_warming.append(exptTdjf - controlTdjf)
    exptnames_used.append(EXPTNAMES.get(expt))
    print(expt,exptTann-controlTann)
      
sys.exit(0)
# plot as a bar chart with ann / djf / jja next to each other

x=np.arange(len(EXPTS))
width=0.25 # the width of the bars
bar_colors=['tab:green','tab:orange','tab:blue']
fig,ax=plt.subplots(layout='constrained')
#plot mean / djf and jja
rects=ax.bar(x-0.25,mean_ann_warming,width,label='Mean Annual',
             color=bar_colors[0])
rects=ax.bar(x,mean_jja_warming,width,label='JJA',
             color=bar_colors[1])
rects=ax.bar(x+0.25,mean_djf_warming,width,label='DJF',
             color=bar_colors[2])


# labels
ax.set_title(LAND_OCN_IND + ' Warming relative to PI experiment')
ax.set_ylabel('degC')
ax.set_xticks(x,exptnames_used,rotation=90)
ax.legend(loc='lower left',ncol=3,bbox_to_anchor=(-0.05, -0.3))
#ax.set_ylim(1.0,4.0)
ax.set_xlim(-0.5,11.5)

# put some colors dividing the experiments
plt.axvspan(-0.5,2.5,facecolor='pink',alpha=0.2)
plt.axvspan(2.5,7.5,facecolor='yellow',alpha=0.2)
plt.axvspan(7.5,9.5,facecolor='purple',alpha=0.2)
plt.axvspan(9.5,11.5,facecolor='green',alpha=0.2)
plt.axvline(x=2.5,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=7.5,ymin=0.0,ymax=1.0,color='black')
plt.axvline(x=9.5,ymin=0.0,ymax=1.0,color='black')
if LAND_OCN_IND == 'Globe':
    textpos = 4.75
if LAND_OCN_IND == 'Land':
    textpos = 8.00
if LAND_OCN_IND == 'Ocean':
    textpos = 4.00

plt.text(-0.25,textpos,'core / extension')
plt.text(4.0, textpos,'CO2 sensitivity')
plt.text(8.0,textpos,'orbital')
plt.text(10.0,textpos,'vegetaton')

#plot the data again because of the background
rects=ax.bar(x-0.25,mean_ann_warming,width,label='Mean Annual warming',
             color=bar_colors[0])
rects=ax.bar(x,mean_jja_warming,width,label='Mean JJA warming',
             color=bar_colors[1])
rects=ax.bar(x+0.25,mean_djf_warming,width,label='Mean DJF warming',
             color=bar_colors[2])

if LAND_OCN_IND == 'Globe':
    plt.savefig('plots/warming_all_models.eps')
else:
    plt.savefig('plots/warming_all_models_' + LAND_OCN_IND + '.eps')
plt.show()
           
           
