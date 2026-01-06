#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.09.2019 by Julia

As a first test for looking at the climate extremes we will plot 
 standard deviation EOI400
 =========================
 standard deviation E280

firstly for temperature. 

If the ratio is approximately 1.0 then the standard deviation is the same.  If it is >1 then there is more interannual variability in the Pliocene (which could be worrying).  If it is less than 1 there is less interannual variability in the Pliocene.


"""

import os
import sys
import numpy as np
#import matplotlib as mp
import matplotlib.pyplot as plt
#from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#import netCDF4
#from netCDF4 import Dataset, MFDataset
import iris
import iris.quickplot as qplt
import iris.analysis.cartography
import iris.coord_categorisation

#os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
from mpl_toolkits.basemap import Basemap, shiftgrid




def getmodelfield(modelname, period):
    """
    get the mean values from the model data
    inputs: modelname (ie HadCM3)
            period (likely EOI400 or E280)
    returns:  a cube contatining the mean data from the model

    """

    modfile = (FILESTART + 'regridded100/' + modelname + '/' +
               period + '.' + FIELDNAME + '.sd_month.nc')

    tempcube = iris.load(modfile)
    cube2 = tempcube[0]
    cube2.units = UNITS
    print(cube2)
    print(cube2.coord('time'))
    
    cube = cube2[MONTH_REQ, :, :]
   
    #this will make all the dimensions of all the cubes match.


    for coord in cube.coords():
        name = coord.standard_name
        if name not in ['latitude', 'longitude']:
            if name is None:
                if coord.long_name is None:
                    cube.remove_coord(coord.var_name)
                else:
                    cube.remove_coord(coord.long_name)
            else:
                cube.remove_coord(name)

    for coord in cube.coords():   # now this will be longitude or latitude
        coord.points = coord.points.astype('float32')
        coord.var_name = coord.standard_name
        coord.long_name = coord.standard_name

    return cube


class Plotalldata:
    """
    This will plot the data from the timeperiod (ie mpwp or pi)
    """
    def __init__(self, timeperiod, anom_cubes):
        self.nmodels = len(MODELNAMES)
        self.filestart = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/'+ 
                          'extremes/sd_ratio/' + timeperiod + '_' + 
                          FIELDNAME + '_' + MONTHNAMES.get(MONTH_REQ))
        self.timeperiod = timeperiod
        self.anom_cubes = anom_cubes

        if timeperiod == 'Ratio':
            self.valmin = 0.5
            self.valmax = 1.6
            self.diff = 0.1
            self.colormap = 'RdBu_r'
            self.cbarlabel = 'ratio'
        else:            
            self.valmin = 0.0
            self.valmax = 2.0
            self.diff = 0.2
            self.colormap = 'RdBu_r'
            self.cbarlabel = 'sdev'


    def plotdata(self):
        """
        this will plot all the cubes to a .eps or .png file
        input anom_cubes : a list of cubes containing the anomalies from the mean
        """


        fig = plt.figure(figsize=(5.5, 4.25))
        for i in range(0, self.nmodels):

            cubedata = self.anom_cubes[i].data
            latitudes = self.anom_cubes[i].coord('latitude').points
            lon = self.anom_cubes[i].coord('longitude').points
            datatoplot, longitudes = (shiftgrid(180., cubedata,
                                                lon, start=False))
            #if (np.mod(i, 8) + 1) == 1:
            #    title_ = (MODELNAMES[i] + ':' +
            #              self.timeperiod + ' (model - MMM)')
            #else:
            #    title_ = (MODELNAMES[i])

            title_ = (MODELNAMES[i])
            self.plotmap(i, title_,
                         datatoplot, longitudes, latitudes, fig)


        return

    def plotmap(self, i, titlename, datatoplot, longitudes, latitudes, fig):
        """
        will plot the data in a map format

        """

        xplot = 4
        yplot = 4


        plotpos = np.mod(i, xplot * yplot) + 1
        plt.subplot(xplot, yplot, plotpos)
        lons, lats = np.meshgrid(longitudes, latitudes)

        map = Basemap(llcrnrlon=-180.0, urcrnrlon=180.0,
                      llcrnrlat=-90.0, urcrnrlat=90.0,
                      projection='cyl', resolution='l')

        #map.drawmapboundary
        x, y = map(lons, lats)
        map.drawcoastlines(linewidth=0.5)

        V = np.arange(self.valmin, self.valmax, self.diff)
        cs = map.contourf(x, y, datatoplot, V, cmap=self.colormap,
                          extend='both')
        plt.title(titlename,fontsize=8)


        if plotpos == (xplot * yplot) or (i + 1) == self.nmodels:
             # Shrink current axis by 20% and put a legend to the right
            plt.subplots_adjust(left=0.05, bottom=0.1, right=0.82, top=0.9,
                                wspace=0.1, hspace=0.0)

            cb_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
           
            cbar = fig.colorbar(cs, cax=cb_ax, orientation='vertical')
            #cbar = plt.colorbar(fig, orientation='horizontal')
            #fig.colorbar(fix, ax=axs[:, col], shrink=0.6)
            cbar.ax.tick_params(labelsize=8)
            print('plotted colorbar')
            #plt.show()
            #plt.tight_layout()
            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                       + '.eps')
            plt.savefig(fileout, bbox_inches='tight')

            fileout = (self.filestart + np.str(np.int(np.ceil(i/8)))
                       + '.pdf')

            plt.savefig(fileout, bbox_inches='tight')
            plt.close()

def plot_mmm(mean_ratio, mean_pi, mean_plio, ocn_mask):
    """
    plots the multimodel mean we will mask out the ocean if necessary
    """
    land_ocn = {'y' : '_land', 'n': '_globe'}

    if ocn_mask == 'y':
       tempcube=iris.load_cube('/nfs/hera1/earjcti/regridded/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc')
       cubegrid=iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
       exptlsmcube=tempcube.regrid(cubegrid,iris.analysis.Linear())
      
       tempcube=iris.load_cube(FILESTART+'regridded/PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc')
       cntllsmcube=tempcube.regrid(cubegrid,iris.analysis.Linear())

       mean_ratio.data.mask = (exptlsmcube.data - 1.0) * (-1.0)
       mean_pi.data.mask = (cntllsmcube.data - 1.0) * (-1.0)
       mean_plio.data.mask = (exptlsmcube.data - 1.0) * (-1.0)
    

    qplt.contourf(mean_ratio, levels=np.arange(0.5, 1.6, 0.1), 
                  cmap='RdBu_r',extend='both')
    plt.title(MONTHNAMES.get(MONTH_REQ) + 'Ratio:  Plio /  PI st_dev : ' + FIELDNAME)
    plt.gca().coastlines()
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/sd_ratio/ratio_MMM_' + FIELDNAME + land_ocn.get(ocn_mask) + '_' + MONTHNAMES.get(MONTH_REQ) + '.eps', bbox_inches='tight')
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/sd_ratio/ratio_MMM_' + FIELDNAME + land_ocn.get(ocn_mask) + '_' + MONTHNAMES.get(MONTH_REQ) + '.png', bbox_inches='tight')
    plt.close()

  
    qplt.contourf(mean_pi, levels=np.arange(0.0, 2.1, 0.1), 
                  cmap='RdBu_r',extend='max')
    plt.title(MONTHNAMES.get(MONTH_REQ) + 'PI st_dev : ' + FIELDNAME)
    plt.gca().coastlines()
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/sd_ratio/pi_sdev_MMM_' + FIELDNAME + land_ocn.get(ocn_mask) + '_' + MONTHNAMES.get(MONTH_REQ) + '.eps', bbox_inches='tight')
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/sd_ratio/pi_sdev_MMM_' + FIELDNAME + land_ocn.get(ocn_mask) + '_' + MONTHNAMES.get(MONTH_REQ) + '.png', bbox_inches='tight')
    plt.close()

    qplt.contourf(mean_plio, levels=np.arange(0.0, 2.1, 0.1), 
                  cmap='RdBu_r',extend='max')
    plt.title(MONTHNAMES.get(MONTH_REQ) + 'Plio st_dev : ' + FIELDNAME)
    plt.gca().coastlines()
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/sd_ratio/plio_sdev_MMM_' + FIELDNAME + land_ocn.get(ocn_mask) + '_' + MONTHNAMES.get(MONTH_REQ) + '.eps', bbox_inches='tight')
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/sd_ratio/plio_sdev_MMM_' + FIELDNAME + land_ocn.get(ocn_mask) + '_' + MONTHNAMES.get(MONTH_REQ) + '.png', bbox_inches='tight')
    plt.close()

  
   
     

def get_mean(pi_sdev_cubes, plio_sdev_cubes):
    """
    this will get the mean standard deviation from all the models and
    find the ratio between them.
    Note to find the mean standard deviation we want to add up all the variances    and then squareroot
    """
    count=0
    for i, cube in enumerate(pi_sdev_cubes):
        if i == 0:
            pi_var_cube = cube * cube
            count=count+1
        else:
            print(pi_var_cube)
            print(cube)
            pi_var_cube = pi_var_cube + (cube * cube)
            count=count+1
    pi_var_cube = pi_var_cube / count
    mean_pi_sd = pi_var_cube.copy(data=np.sqrt(pi_var_cube.data))
   

    count=0
    for i, cube in enumerate(plio_sdev_cubes):
        if i == 0:
            plio_var_cube = cube * cube
            count=count+1
        else:
            plio_var_cube = plio_var_cube + (cube * cube)
            count=count+1
    plio_var_cube = plio_var_cube / count
    mean_plio_sd = plio_var_cube.copy(data=np.sqrt(plio_var_cube.data))
   
  
    return mean_plio_sd / mean_pi_sd, mean_pi_sd, mean_plio_sd 
##########################################################
# main program
# set up variable information
FIELDNAME = 'NearSurfaceTemperature'
UNITS = 'Celsius'
#FIELDNAME = 'SST'
#UNITS = 'Celsius'
#FIELDNAME = 'TotalPrecipitation'
#UNITS = 'mm/day'
#FIELDNAME = 'SST'
LINUX_WIN = 'l'

MONTHNAMES = {0:'Jan',1:'Feb',2:'Mar',3:'Apr',4:'May',5:'Jun',6:'Jul',7:'Aug',8:'Sep',9:'Oct',10:'Nov',11:'Dec',}
MONTH_REQ = 6

if LINUX_WIN == 'l':
    FILESTART = '/nfs/hera1/earjcti/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'

MODELNAMES = ['CESM2', 'IPSLCM6A', 'COSMOS', 
            'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
            'MIROC4m', 'IPSLCM5A2', 'HadCM3',
            'GISS2.1G', 'CCSM4', 
            'CCSM4-Utr', 'CCSM4-UoT', 
            'NorESM-L', 'MRI2.3', 'NorESM1-F'
             ]

#MODELNAMES = ['NorESM-L']



# set up cubelists to store data
mpwp_sdev_cubes = iris.cube.CubeList([])
pi_sdev_cubes = iris.cube.CubeList([])
ratio_sdev_cubes = iris.cube.CubeList([])

#################################################
# get standard deviation data

for model, modelname in enumerate(MODELNAMES):
    model_plio_cube = getmodelfield(modelname, 'EOI400')
    model_pi_cube = getmodelfield(modelname, 'E280')
  
    if modelname == 'EC-Earth3.1' and FIELDNAME == 'SST':
       model_pi_cube.coord('latitude').bounds = None
       model_pi_cube.coord('longitude').bounds = None

    model_anom_cube = model_plio_cube - model_pi_cube

    mpwp_sdev_cubes.append(model_plio_cube)
    pi_sdev_cubes.append(model_pi_cube)
    ratio_sdev_cubes.append(model_plio_cube / model_pi_cube)

##################################################
# plot the cubes for the model anomalies relative to the mean

#obj = Plotalldata('mPWP', mpwp_sdev_cubes)
#obj.plotdata()

#obj = Plotalldata('PI', pi_sdev_cubes)
#obj.plotdata()

#obj = Plotalldata('ratio', ratio_sdev_cubes)
#obj.plotdata()


################################################################
# get mean ratios 
# note that the mean standard deviation is the square root of the sum
# of the variances

mean_ratio, mean_pi_sd, mean_plio_sd = get_mean(pi_sdev_cubes, mpwp_sdev_cubes)
plot_mmm(mean_ratio, mean_pi_sd, mean_plio_sd,'y')
plot_mmm(mean_ratio, mean_pi_sd, mean_plio_sd,'n')


#
