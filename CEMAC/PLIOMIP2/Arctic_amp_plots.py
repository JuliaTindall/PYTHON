#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
#Created on 05/03/2020


#@author: earjcti

This program will plot the polar amplification diagrams
requested by IPCC in particular

Arctic amplification: I suggest reporting the
average mean annual land + sea temperature for six,
30 deg latitude bands,
 or at least the temperatures for 60-90N vs 0-90N.

So we will do:

For land/sea/total the temperature anomaly for the 6 bands
"""
import sys
import warnings
import numpy as np
import matplotlib.pyplot as plt
import iris
import iris.analysis.cartography
import iris.coord_categorisation
import re

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

def get_mean_data(model, expt):
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
                + '/' + expt + '.' + FIELD + '.allmean.nc')

    cube = iris.load_cube(filename)
    cube2 = resort_coords(cube)

    cube2.coord('latitude').guess_bounds()
    cube2.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube2)

    return cube2, grid_areas


def get_region(band_upper, band_lower, cube, mask, grid_areas):
    """
    Gets the average temeprature over land within the
    bounded range

    Parameters
    ----------
    band_upper : scalar value denoting the upper latitude of
                 the band
    band_lower : scalar value denoting the lower latitude of
                 the band
    cube : the cube containing average temperatures that
           we want to get the land temperature over
    mask : numpy array containing the mask we want.  Could be a land mask
           a sea mask or ones everywhere (ie all points)
    grid_areas : the size of each gridpoint for weighting

    Returns
    -------
    bound_land_avg : scalar containing the average land
                    temperature over the bounded region

    """

    grid_areas_mask = grid_areas * mask
    grid_areas_band = np.zeros(grid_areas.shape)
    lats = cube.coord('latitude').points

    for j, lat in enumerate(lats):
        if band_lower <= lat <= band_upper:
            grid_areas_band[j, :] = grid_areas_mask[j, :]

    bound_mask_avg = cube.collapsed(['longitude', 'latitude'],
                                    iris.analysis.MEAN,
                                    weights=grid_areas_band)

    return bound_mask_avg.data

def get_pliomip1_data():
    """
    Gets the pliomip1 data for each of the bands

    Returns
    -------
    pliomip1 data

    """
    if LINUX_WIN == 'l':
        PLIOMIP1_FILE = (FILESTART + '/PLIOMIP/means_for_' + FIELD + '.txt')
    else:
        PLIOMIP1_FILE = FILESTART + 'PLIOMIP1/means_for_' + FIELD + '.txt'
    f1 = open(PLIOMIP1_FILE)
    lines = f1.readlines()
    lines[:] = [line.rstrip('\n') for line in lines]
    # find line index which has the title 'modelname, latband mean '
    string = 'modelname, latband mean'
    for i, line in enumerate(lines):
        if line[0:23] == string:
            index = i
        
    # get bands by splitting the line
    bands_line = lines[index + 1]
    bands_str_array = bands_line[10:]
    print(bands_str_array)
    res = bands_str_array.split("], [")
    res = [x.strip('[') for x in res]
    res = [x.strip(']') for x in res]
    nbands = len(res)
    bands_array = np.zeros((nbands, 2))
    for i, x in enumerate(res):
        x1, x2 = x.split(',')
        bands_array[i, 0] = x1
        bands_array[i, 1] = x2
     
    # get the data from the next lines find the mean, min and max for each band
    # for the anomaly only
    minval = np.zeros(nbands)
    minval[:] = 100.
    maxval = np.zeros(nbands)
    meanval = np.zeros(nbands)
    
    for i in range(index + 2, len(lines)):
        line = lines[i]
        if line[0:9] == 'modelname':
            break      
        modname, eoi400, e280, anom = line.split(',')

        eoi400_val = np.array(eoi400.strip('[]').split(), dtype=float)
        e280_val = np.array(e280.strip('[]').split(), dtype=float)
        
        anom_val = np.array(anom.strip('[]').split(), dtype=float)
        
        for i, anom in enumerate(anom_val):
            if anom < minval[i]:
                minval[i] = anom
            if anom > maxval[i]:
                maxval[i] = anom
        if modname == 'MEAN':
            meanval = anom_val
        
    return meanval, minval, maxval
  
   
   
   
def plot_temp_by_lat(anomaly, uppervals, lowervals, plottype, fileoutstart,
                     mean_p1, min_p1, max_p1):
    """
    plot the temperature anomaly vs the region on one plot

    Parameters
    ----------
    anomaly : temperature anomaly to plot np.shape= nmodels, nbounds
    uppervals : the upper limit of the boundary range (np.shape = nbounds)
    lowervals : the lower limit of the boundary range (np.shape = nbounds)
    plottype : 'Land' 'Ocean' or '' if empty it is all surface types

    Returns
    -------
    None.

    """

    titlename = 'MPWP - PI ' + plottype + ' SAT anomaly'
    latns = {-90.0 : '90S', -60.0: '60S', -30.0 : '30S', 0.0 : '0N',
             30.0 : '30N', 60.0 : '60N', 90.0 : '90N'}


    labels = []
    for i, upval in enumerate(uppervals):
        labels.append(latns.get(upval) + '-' + latns.get(lowervals[i]))

    ax = plt.subplot(1, 1, 1)
    for i, model in enumerate(MODELNAMES):
        if i < len(MODELNAMES) / 2.0:
            ax.plot(anomaly[i, :], labels, label=model)
        else:
            ax.plot(anomaly[i, :], labels, label=model, linestyle='dashed')


    ax.plot(np.mean(anomaly, axis=0), labels, color='black',
            linestyle='dashed',
            linewidth=2, label='avg')
    plt.title(titlename)
    plt.xlabel(UNITS)

    if PLIOMIP1 == 'y' and plottype == '':
        ax.plot(mean_p1, labels, label = 'PlioMIP1', color='black',
                linestyle = 'dotted', linewidth=2)
        ax.fill_betweenx(labels, min_p1, max_p1, alpha=0.2, 
                        color="grey")
       
    FILETXT = (FILESTART + '/regridded/alldata/data_for_supp_fig2_' 
               + plottype + '.txt')

    txtout = open(FILETXT, "w+")
    
    writedata = 'modelname'
    for j in range(0, len(labels)):
        writedata = writedata + ',' + labels[j]
    writedata = writedata + '\n'
    txtout.write(writedata)

    for i, mod in enumerate(MODELNAMES):
        writedata = mod 
        for j in range(0, len(labels)):
            writedata = writedata + ',' + (np.str(np.around(anomaly[i,j],2)))
        writedata = writedata + '\n'
        txtout.write(writedata)
    txtout.close()
   

    # Shrink current axis by 20% and put a legend to the right
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    fileout = (fileoutstart + '/polar_amplification'  + plottype + '_anomaly.eps')
    plt.savefig(fileout)
    fileout = (fileoutstart + '/polar_amplification'  + plottype + '_anomaly.pdf')
    plt.savefig(fileout)
    plt.show()
    plt.close()


#####################################
def main():
    """
    Tha main control of the program to plot the
    polar amplification by temperature by latitude band

    """

    if LINUX_WIN == 'w':
        fileoutstart = ('C:\\Users\\julia\\OneDrive\\WORK\\DATA\\' +
                        'regridded\\allplots\\' + FIELD + '\\')
        exptlsm = (FILESTART + 'regridded/PlioMIP2_Boundary_conds' +
                   '/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc')
        cntllsm = (FILESTART+'regridded/PlioMIP2_Boundary_conds' +
                   '/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc')
    else:
        fileoutstart = '/nfs/hera1/earjcti/regridded/allplots/' + FIELD + '/'
        dataoutstart = '/nfs/hera1/earjcti/regridded/alldata/'
        exptlsm = (FILESTART + 'PlioMIP2_Boundary_conds/Plio_enh' +
                   '/Plio_enh/Plio_enh_LSM_v1.0.nc')
        cntllsm = (FILESTART+'PlioMIP2_Boundary_conds/Modern_std/' +
                   'Modern_std/Modern_std_LSM_v1.0.nc')


    ########################################################
    # setup: get the lsm for the land sea contrast plot

    exptland, exptsea = get_lsm(exptlsm)
    cntlland, cntlsea = get_lsm(cntllsm)
    maskall = np.ones(np.shape(exptland))


    #########################################################
    # need to get data from annual mean plot

    bandmax = [-60., -30., 0., 30., 60., 90.]
    bandmin = [-90., -60., -30., 0., 30., 60.]


    land_anomaly = np.zeros((len(MODELNAMES), len(bandmax)))
    sea_anomaly = np.zeros((len(MODELNAMES), len(bandmax)))
    bands_anomaly = np.zeros((len(MODELNAMES), len(bandmax)))

    land_expt = np.zeros((len(MODELNAMES), len(bandmax)))
    sea_expt = np.zeros((len(MODELNAMES), len(bandmax)))
    bands_expt = np.zeros((len(MODELNAMES), len(bandmax)))

    land_cntl = np.zeros((len(MODELNAMES), len(bandmax)))
    sea_cntl = np.zeros((len(MODELNAMES), len(bandmax)))
    bands_cntl = np.zeros((len(MODELNAMES), len(bandmax)))
    
    land_expt_avg = np.zeros((len(MODELNAMES)))
    sea_expt_avg = np.zeros((len(MODELNAMES)))
    expt_avg = np.zeros((len(MODELNAMES)))

    land_cntl_avg = np.zeros((len(MODELNAMES)))
    sea_cntl_avg = np.zeros((len(MODELNAMES)))
    cntl_avg = np.zeros((len(MODELNAMES)))

    

    for modelno, modeluse in enumerate(MODELNAMES):

        # get mean data
        (exptcube, grid_areas_expt) = get_mean_data(modeluse, EXPTNAME)
        (cntlcube, grid_areas_cntl) = get_mean_data(modeluse, CNTLNAME)

        # get data within bounds
        for i, band_upper in enumerate(bandmax):
            band_lower = bandmin[i]
            land_expt[modelno, i] = get_region(band_upper, band_lower,
                                               exptcube, exptland,
                                               grid_areas_expt)

            land_cntl[modelno, i] = get_region(band_upper, band_lower,
                                               cntlcube, cntlland,
                                               grid_areas_cntl)

            sea_expt[modelno, i] = get_region(band_upper, band_lower,
                                              exptcube, exptsea,
                                              grid_areas_expt)

            sea_cntl[modelno, i] = get_region(band_upper, band_lower,
                                              cntlcube, cntlsea,
                                              grid_areas_cntl)

            bands_expt[modelno, i] = get_region(band_upper, band_lower,
                                                exptcube, maskall,
                                                grid_areas_expt)

            bands_cntl[modelno, i] = get_region(band_upper, band_lower,
                                                cntlcube, maskall,
                                                grid_areas_cntl)
            
        #get average data for calculating polar amplification
        land_expt_avg[modelno] = get_region(90.0, -90.0, exptcube, exptland,
                                            grid_areas_expt)
        land_cntl_avg[modelno] = get_region(90.0, -90.0, cntlcube, cntlland,
                                            grid_areas_cntl)
        sea_expt_avg[modelno] = get_region(90.0, -90.0, exptcube, exptsea,
                                            grid_areas_expt)
        sea_cntl_avg[modelno] = get_region(90.0, -90.0, cntlcube, cntlsea,
                                            grid_areas_cntl)
        expt_avg[modelno] = get_region(90.0, -90.0, exptcube, maskall,
                                            grid_areas_expt)
        cntl_avg[modelno] = get_region(90.0, -90.0, cntlcube, maskall,
                                            grid_areas_cntl)

    land_anomaly = land_expt - land_cntl
    sea_anomaly = sea_expt - sea_cntl
    bands_anomaly = bands_expt - bands_cntl
    
    ##############################################################
    #  get pliomip1 data if required.  Note we will only get annual mean
    
    if PLIOMIP1 == 'y':
        pliomip1_mean, pliomip1_min, pliomip1_max = get_pliomip1_data()
    else:
        pliomip1_mean = 0
        pliomip1_min = 0
        pliomip1_max = 0

    #print polar amplification
    print('SH polar amplification')
    print('model, land amplification, sea amp, all amp')
    #for i, model in enumerate(MODELNAMES):
    #    print('bands',bandmax[0],bandmin[0])
    #    print(model, land_anomaly[i, 0] / (land_expt_avg[i] - land_cntl_avg[i]), 
    #          sea_anomaly[i, 0] / (sea_expt_avg[i] - sea_cntl_avg[i]),
    #          bands_anomaly[i,0] / (expt_avg[i] - cntl_avg[i]))
    all_sh_land_amp = [x1 / x2 for (x1, x2) in zip(land_anomaly[:,0], 
                                                   land_expt_avg - land_cntl_avg)]
    all_sh_sea_amp = [x1 / x2 for (x1, x2) in zip(sea_anomaly[:,0], 
                                                  sea_expt_avg - sea_cntl_avg)]
    all_sh_amp = [x1 / x2 for (x1, x2) in zip(bands_anomaly[:,0], 
                                              expt_avg - cntl_avg)]
    print('mean', 
          np.mean(land_anomaly[:, 0]) / np.mean(land_expt_avg[:] - land_cntl_avg[:]), 
          np.mean(sea_anomaly[:, 0]) / np.mean(sea_expt_avg[:] - sea_cntl_avg[:]),
          np.mean(bands_anomaly[:,0] / np.mean(expt_avg[:] - cntl_avg[:])))
    print('median', np.median(all_sh_land_amp), np.median(all_sh_sea_amp), 
          np.median(all_sh_amp))
    
    print('checking land',  np.mean(land_anomaly[:, 0]), np.mean(land_expt[:,0]),
          np.mean(land_cntl[:,0]))
          #,np.mean(land_expt_avg[:])-
          #np.mean(land_cntl_avg[:]))
    print('checking sea',  np.mean(sea_anomaly[:, 0]), np.mean(sea_expt[:,0]),
          np.mean(sea_cntl[:,0]))#,np.mean(sea_expt_avg[:])-
        #  np.mean(sea_cntl_avg[:]))
    print('checking avg',  np.mean(bands_anomaly[:, 0]),np.mean(bands_expt[:,0]),
          np.mean(bands_cntl[:,0]))#,np.mean(expt_avg[:])-
        #  np.mean(cntl_avg[:]))
      
    print(' ')    
    print('NH polar amplification')
    all_nh_land_amp = [x1 / x2 for (x1, x2) in zip(land_anomaly[:,5], 
                                                   land_expt_avg - land_cntl_avg)]
    all_nh_sea_amp = [x1 / x2 for (x1, x2) in zip(sea_anomaly[:,5], 
                                                  sea_expt_avg - sea_cntl_avg)]
    all_nh_amp = [x1 / x2 for (x1, x2) in zip(bands_anomaly[:,5], 
                                              expt_avg - cntl_avg)]
    print('model, land amp, sea amp, all amp')
    #for i, model in enumerate(MODELNAMES):
    #    print('bands',bandmax[5],bandmin[5])
    #    print(model, land_anomaly[i, 5] / (land_expt_avg[i] - land_cntl_avg[i]), 
    #          sea_anomaly[i, 5] / (sea_expt_avg[i] - sea_cntl_avg[i]),
    #          bands_anomaly[i,5] / (expt_avg[i] - cntl_avg[i]))
    print('mean', 
          np.mean(land_anomaly[:, 5]) / np.mean(land_expt_avg[:] - land_cntl_avg[:]), 
          np.mean(sea_anomaly[:, 5]) / np.mean(sea_expt_avg[:] - sea_cntl_avg[:]),
          np.mean(bands_anomaly[:,5] / np.mean(expt_avg[:] - cntl_avg[:])))
    print('median', np.median(all_nh_land_amp), np.median(all_nh_sea_amp), 
          np.median(all_nh_amp))

    # plot everything on one plot.

    plot_temp_by_lat(land_anomaly, bandmax, bandmin, 'Land', 
                     fileoutstart, 0, 0, 0)
    plot_temp_by_lat(sea_anomaly, bandmax, bandmin, 'Ocean', 
                     fileoutstart, 0, 0, 0)
    plot_temp_by_lat(bands_anomaly, bandmax, bandmin, '', 
                     fileoutstart, pliomip1_mean, pliomip1_min, pliomip1_max)



##########################################################
# DEFINITIONS

LINUX_WIN = 'l'
if LINUX_WIN == 'w':
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
else:
    FILESTART = '/nfs/hera1/earjcti/'


MODELNAMES = ['CESM2', 'IPSLCM6A', 'COSMOS', 
              'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
              'MIROC4m', 'IPSLCM5A2', 'HadCM3',
              'GISS2.1G', 'CCSM4', 
              'CCSM4-Utr', 'CCSM4-UoT', 
              'NorESM-L', 'MRI2.3', 'NorESM1-F'
              ]

PLIOMIP1 = 'y' # overplot PlioMIP1 models.

#MODELNAMES=['HadCM3','NorESM-L']
FIELD = 'NearSurfaceTemperature'
UNITS = 'degC'

#fieldnames=['TotalPrecipitation']
#units=['mm/day']
EXPTNAME = 'EOI400'
CNTLNAME = 'E280'

if FIELD == 'NearSurfaceTemperature':
    FILEOUT = FILESTART + 'regridded/alldata/data_for_arctic_amplification.txt'


main()
