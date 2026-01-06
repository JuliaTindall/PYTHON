#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
#Created on 05/03/2020


#@author: earjcti

This program will plot the land temperature relationships
requested by IPCC in particular

- Land vs sea: I suggest reporting both (1) all sea vs GMAT, and 
(2) sea from 60S - 60N vs GMST. 
These metrics are also needed to help inform GMAT reconstructions from the proxy data, 
which are primarily marine-based.


"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.analysis.cartography
import iris.coord_categorisation
from scipy import stats

warnings.filterwarnings("ignore")

def resort_coords(cube):
    """
    this will make all the dimensions of the cube match.
    """

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

def get_lsm(filein):
    """
    returns a  numpy array of land points and one of sea points
    """
    tempcube = iris.load_cube(filein)
    cubegrid = iris.load_cube('one_lev_one_deg.nc')
    lsmcube = tempcube.regrid(cubegrid, iris.analysis.Linear())
    landpoints = lsmcube.data
    seapoints = (lsmcube.data - 1.0) * (-1.0)

    return landpoints, seapoints

def get_mean_data(model, expt, field):
    """
    gets the cube of mean data for a single model

    Parameters
    ----------
    model : the name of the model we are interested in
    expt : whether it is the experiment or the control

    Returns
    -------
    a cube with the mean data from this file
    grid_areas = the size of the grid for averaging
    """

    filename = (FILESTART + 'regridded/' + model
                + '/' + expt + '.' + field + '.allmean.nc')

    print(field, filename)
   
    cube = iris.load_cube(filename)
    cube2 = resort_coords(cube)

    cube2.coord('latitude').guess_bounds()
    cube2.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube2)

    return cube2, grid_areas


def get_region(latmin, latmax, cube, mask, grid_areas):
    """
    Gets the average temeprature over land within the
    bounded range

    Parameters
    ----------
    latmin, latmax : the range of latitudes we are extracting data from
    cube : the cube containing average temperatures t
    mask : numpy array containing the mask we want.  Could be a land mask
           a sea mask or ones everywhere (ie all points)
    grid_areas : the size of each gridpoint for weighting

    Returns
    -------
    region_avg : scalar containing the average temperature over the
                     required region

    """

    grid_areas_mask = grid_areas * mask
    grid_areas_band = np.zeros(grid_areas.shape)
    lats = cube.coord('latitude').points

    for j, lat in enumerate(lats):
        if latmin <= lat <= latmax:
            grid_areas_band[j, :] = grid_areas_mask[j, :]

    region_avg = cube.collapsed(['longitude', 'latitude'],
                                iris.analysis.MEAN,
                                weights=grid_areas_band)

    return region_avg.data

def scatter_sea_vs_global(allanom, seaanom, plottype, txtfile):
    """
    plot the global temperature anomaly vs the 
    ocean temperature anomaly for theregion on one plot
    also outputs the regression equation

    Parameters
    ----------
    allanom : temperature anomaly from the globe
    seaanom : temperature anomaly from the sea
    plot type: '' - global or latitude range

    Returns
    -------
    None.

    """

    titlename = 'MPWP - PI: global vs ocean ' + plottype + ' temperature '
    ax = plt.subplot(1, 1, 1)

    for i, model in enumerate(MODELNAMES):
        if i % 4 == 0: # i divides 4 with no remainder
            ax.scatter(seaanom[i], allanom[i], label = model) 
        elif i % 4 == 1 :
            ax.scatter(seaanom[i], allanom[i], label = model, marker='^') 
        elif i % 4 == 2 :
            ax.scatter(seaanom[i], allanom[i], label = model, marker='<') 
        else:
            ax.scatter(seaanom[i], allanom[i], label = model, marker='v') 
    #plt.title(titlename)
    plt.xlabel('ocean temperature ' + plottype + ' anomaly (' + UNITS + ')')
    plt.ylabel('global SAT ' + ' anomaly (' + UNITS + ')')


    # Shrink current axis by 20% and put a legend to the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    
    fileout = (FILEOUTSTART + '/sea_SAT'  + plottype + '_anomaly.eps')
    plt.savefig(fileout)
    fileout = (FILEOUTSTART + '/sea_SAT'  + plottype + '_anomaly.pdf')
    #plt.show()
    plt.savefig(fileout)
    plt.close()


    # print out the regression equation and the data to a file
    print(plottype)
    print('=======')
    (slope, intercept, r_value, 
            p_value, std_err) = stats.linregress(seaanom, allanom)
    print('GMSAT = ' + np.str(np.around(slope, 3)) + 
          'x OSAT + ' + np.str(np.around(intercept, 3)))
    print('pvalue = ' + np.str(np.around(p_value, 3)) + 
          ' rvalue = ' + np.str(np.round(r_value, 3)) + 
          ' rsq = ' + np.str(np.round(r_value * r_value, 3)))
    print('   ')

    # print out the values to a file
    
    txtfile.write("modelname, global mean SAT, oceanSAT " + plottype + '\n')
    for i, model in enumerate(MODELNAMES):
        txtfile.write((model + ',' + np.str(np.around(allanom[i],3)) + 
                       ',' + np.str(np.around(seaanom[i],3)) + '\n'))
    
      

#####################################
def main():
    """
    Tha main control of the program to plot the
    polar amplification by temperature by latitude band

    """

    if LINUX_WIN == 'w': 
        exptlsm = (FILESTART + 'regridded/PlioMIP2_Boundary_conds' +
                   '/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc')
        cntllsm = (FILESTART+'regridded/PlioMIP2_Boundary_conds' +
                   '/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc')
    else:
        exptlsm = (FILESTART + 'PlioMIP2_Boundary_conds/Plio_enh' +
                   '/Plio_enh/Plio_enh_LSM_v1.0.nc')
        cntllsm = (FILESTART+'PlioMIP2_Boundary_conds/Modern_std/' +
                   'Modern_std/Modern_std_LSM_v1.0.nc')


    ########################################################
    # setup: get the lsm for the land sea contrast plot

    exptland, exptsea = get_lsm(exptlsm)
    cntlland, cntlsea = get_lsm(cntllsm)
    fullmask = np.ones(np.shape(exptsea))
    

    #########################################################
    # need to get data from annual mean plot

    sea_anomaly = np.zeros((len(MODELNAMES)))
    sea_6060_anomaly = np.zeros((len(MODELNAMES))) # from 60N-60S
    all_anomaly = np.zeros((len(MODELNAMES)))
    
    for modelno, modeluse in enumerate(MODELNAMES):

        # get mean data
        (exptcube_SAT, grid_areas_expt) = get_mean_data(modeluse, EXPTNAME, FIELD_SAT)
        (cntlcube_SAT, grid_areas_cntl) = get_mean_data(modeluse, CNTLNAME, FIELD_SAT)
        (exptcube_SST, grid_areas_expt) = get_mean_data(modeluse, EXPTNAME, FIELD_SST)
        (cntlcube_SST, grid_areas_cntl) = get_mean_data(modeluse, CNTLNAME, FIELD_SST)
        
        # get all data
        expt_data = get_region(-90.0, 90.0, exptcube_SAT, 
                               fullmask, grid_areas_expt)
        cntl_data = get_region(-90.0, 90.0, cntlcube_SAT, 
                               fullmask, grid_areas_cntl)
        print(expt_data, cntl_data)
        all_anomaly[modelno] = expt_data - cntl_data
        

        # get sea anomaly
        expt_data = get_region(-90.0, 90.0, exptcube_SST, 
                               exptsea, grid_areas_expt)
        cntl_data = get_region(-90.0, 90.0, cntlcube_SST, 
                               cntlsea, grid_areas_cntl)
        sea_anomaly[modelno] = expt_data - cntl_data

        # get sea anomaly from 60N-60S
        expt_data = get_region(-60.0, 60.0, exptcube_SST, 
                               exptsea, grid_areas_expt)
        cntl_data = get_region(-60.0, 60.0, cntlcube_SST, 
                               cntlsea, grid_areas_cntl)
        sea_6060_anomaly[modelno] = expt_data - cntl_data

   # plot everything on one plot and write results to a text file
    filetext = open((DATAOUTSTART + '/data_for_sea_SAT_anomaly.txt'), "w+")
   
    scatter_sea_vs_global(all_anomaly, sea_anomaly, '', filetext)
    scatter_sea_vs_global(all_anomaly, sea_6060_anomaly, '_60N-60S_', filetext)
    
    filetext.close  

##########################################################
# DEFINITIONS

LINUX_WIN = 'l'
FIELD_SAT = 'NearSurfaceTemperature'
FIELD_SST = 'SST'
UNITS = 'degC'
EXPTNAME = 'EOI400'
CNTLNAME = 'E280'

MODELNAMES=['CESM2', 'IPSLCM6A', 'COSMOS', 
            'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
            'MIROC4m', 'IPSLCM5A2', 'HadCM3',
            'GISS2.1G', 'CCSM4', 
            'CCSM4-Utr', 'CCSM4-UoT', 
            'NorESM-L', 'MRI2.3', 'NorESM1-F'
            ]


if LINUX_WIN == 'w':
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
    FILEOUTSTART = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\' +
                        'regridded\\allplots\\' + FIELD_SAT + '\\')
    DATAOUTSTART = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\' +
                        'regridded\\allplots\\' + FIELD_SAT + '\\')
else:
    FILESTART = '/nfs/hera1/earjcti/'
    FILEOUTSTART = '/nfs/hera1/earjcti/regridded/allplots/' + FIELD_SAT + '/'
    DATAOUTSTART = '/nfs/hera1/earjcti/regridded/alldata/'

FILEOUT = FILESTART + 'regridded/alldata/data_for_sea_SAT_relationships.txt'


main()
