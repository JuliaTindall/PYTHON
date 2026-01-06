#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created on August 2020


#@author: earjcti
#
# This program plot a figure for IPCC.  This includes
# a) MPWP - PI SAT anomaly over land (MMM)
# b) MPWP - PI SST anomaly over ocean (MMM)
# c) data overplotted
# d) Pliocene LSM


#import os
import numpy as np
import pandas as pd
#import scipy as sp
#import cf
import iris
#import iris.util
#import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#import netCDF4
#from mpl_toolkits.basemap import Basemap, shiftgrid
#from netCDF4 import Dataset, MFDataset
#import iris.analysis.cartography
#import iris.coord_categorisation
import iris.quickplot as qplt
import iris.plot as iplt
#import cf_units as unit
#from iris.experimental.equalise_cubes import equalise_attributes
import cartopy.crs as ccrs
#import matplotlib.ticker as mticker
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
#from mpl_toolkits.basemap import Basemap

import sys


def make_cmap(colors, position=None, bit=False):
    '''
    I didn't write this I found it on the web.
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap

def customise_cmap():
    """
    customises colormap
    """
    colors = [(5, 48, 97),(6, 49, 98),(7, 51, 100),(8, 53, 102),
               (9, 55, 104),(11, 57, 106),(12, 59, 108),(13, 61, 110),
               (14, 63, 112),(15, 65, 114),(17, 67, 116),
               (18, 69, 118),(19, 71, 120),(20, 73, 121),(22, 75, 123),
               (23, 77, 125),(24, 79, 127),(25, 81, 129),(26, 82, 131),
               (28, 84, 133),(29, 86, 135),(30, 88, 137),(31, 90, 139),
               (32, 92, 141),(34, 94, 143),(35, 96, 145),
               (36, 98, 146),(37, 100, 148),(39, 102, 150),(40, 104, 152),
               (41, 106, 154),(42, 108, 156),(43, 110, 158),(45, 112, 160),
               (46, 113, 162),(47, 115, 164),(48, 117, 166),(49, 119, 168),
               (51, 121, 170),(52, 123, 171),(53, 125, 173),
               (54, 127, 175),(56, 129, 177),(57, 131, 179),(58, 133, 181),
               (59, 135, 183),(60, 137, 185),(62, 139, 187),(63, 141, 189),
               (64, 143, 191),(65, 145, 193),(67, 147, 195),(69, 148, 195),
               (71, 149, 196),(74, 150, 197),(76, 152, 197),
               (78, 153, 198),(81, 155, 199),(83, 156, 199),(86, 157, 200),
               (88, 159, 201),(90, 160, 202),(93, 161, 202),(95, 163, 203),
               (97, 164, 204),(100, 165, 204),(102, 166, 205),(105, 168, 206),
               (107, 169, 207),(109, 171, 207),(112, 172, 208),(114, 173, 209),
               (116, 175, 209),(119, 176, 210),(121, 177, 211),(124, 179, 211),
               (126, 180, 212),(128, 181, 213),(131, 183, 214),(133, 184, 214),
               (135, 185, 215),(138, 187, 216),(140, 188, 216),(143, 189, 217),
               (145, 191, 218),(147, 192, 219),(150, 193, 219),(152, 195, 220),
               (155, 196, 221),(157, 197, 221),(159, 198, 222),(162, 200, 223),
               (164, 201, 223),(166, 203, 224),(169, 204, 225),(171, 205, 226),
               (174, 207, 226),(176, 208, 227),(178, 209, 228),(181, 211, 228),
               (183, 212, 229),(185, 213, 230),(188, 214, 230),(190, 216, 231),
               (193, 217, 232),(195, 219, 233),(197, 220, 233),(200, 221, 234),
               (202, 223, 235),(204, 224, 235),(207, 225, 236),(209, 227, 237),
               (212, 228, 238),(214, 229, 238),(216, 230, 239),(219, 232, 240),
               (221, 233, 240),(224, 235, 241),(226, 236, 242),(228, 237, 243),
               (231, 239, 243),(233, 240, 244),
               (235, 241, 245),(238, 243, 245),(240, 244, 246),(243, 245, 247),
               (245, 246, 247),(247, 248, 248),(248, 248, 247),(248, 246, 245),
               (247, 243, 243),(247, 242, 241),(246, 240, 238),(246, 238, 236),
               (246, 235, 234),(245, 234, 232),(245, 232, 229),(244, 230, 227),
               (244, 227, 225),(243, 226, 223),(243, 224, 220),(242, 222, 218),
               (242, 220, 216),(241, 218, 214),(241, 216, 211),(240, 214, 209),
               (240, 211, 207),(240, 210, 205),(239, 208, 202),(239, 206, 200),
               (238, 203, 198),(238, 202, 196),(237, 200, 193),(237, 198, 191),
               (236, 195, 189),(236, 194, 187),(235, 192, 184),(235, 190, 182),
               (235, 187, 108),(234, 186, 178),(234, 184, 175),(233, 181, 173),
               (233, 179, 171),(232, 178, 169),(232, 176, 166),(231, 174, 164),
               (231, 172, 162),(230, 170, 160),(230, 168, 157),(230, 166, 155),
               (229, 163, 153),(229, 162, 151),(228, 160, 148),(228, 158, 146),
               (227, 156, 144),(227, 154, 142),(226, 152, 139),(226, 149, 137),
               (225, 147, 135),(225, 146, 133),(224, 144, 130),(224, 142, 128),
               (224, 140, 126),(223, 138, 124),(223, 135, 121),(222, 134, 119),
               (222, 132, 117),(221, 130, 115),(221, 128, 112),(220, 125, 110),
               (220, 124, 108),(219, 121, 106),(219, 120, 103),(219, 118, 101),
               (218, 115, 99),(218, 113, 97),(217, 112, 94),(217, 110, 92),
               (216, 108, 90),(216, 105, 88),(215, 104, 85),(215, 102, 83),
               (214, 100, 81),(214, 97, 79),(214, 96, 76),(211, 94, 76),
               (209, 92, 75),(207, 90, 74),(205, 88, 73),(203, 86, 72),
               (200, 84, 71),(198, 82, 70),(196, 80, 69),(194, 79, 68),
               (192, 77, 67),(190, 75, 67),(187, 73, 66),(185, 71, 65),
               (183, 69, 64),(181, 67, 63),(179, 65, 62),(177, 64, 61),
               (174, 62, 60),(172, 60, 59),(170, 58, 58),(168, 56, 58),
               (166, 54, 57),(163, 52, 56),(161, 50, 55),(159, 48, 54),
               (157, 47, 53),(155, 45, 52),(153, 43, 51),(150, 41, 50),
               (148, 39, 49),(146, 37, 49),(144, 35, 48),(142, 33, 47),
               (140, 32, 46),(137, 30, 45),(135, 28, 44),(133, 26, 43),
               (131, 24, 42),(129, 22, 41),(126, 20, 40),(124, 18, 40),
               (122, 16, 39),(120, 15, 38),(118, 13, 37),(116, 11, 36),
               (113, 9, 35),(111, 7, 34),(109, 5, 33),(107, 3, 32),
               (105, 1, 31),(103, 0, 31)]
    my_cmap = make_cmap(colors, bit=True)

    return my_cmap

def customise_cmap2():
    """
    as customise_cmap but 19 colors only + 2 white in middle added by Julia
    """

    colors = [(5,48,97),(22,75,124),(39,102,151),(56,130,178),
              (80,154,199),(114,173,209),(147,192,219),(181,211,228),
              (215,230,238),(255,255,255),(255,255,255),(255,255,255),
              (242,220,217),(236,192,185),
              (229,163,153),(223,135,121),(216,107,89),(195,80,69),
              (164,53,56),(133,26,43),(103,0,31)]
    my_cmap = make_cmap(colors, bit=True)
    return my_cmap

def get_lsm():
    """
    land sea mask is where the point is ocean in both pliocene and pi
    """
    lsm_pi_cube = iris.load_cube(LSM_PI_FILE)
    lsm_plio_cube = iris.load_cube(LSM_PLIO_FILE)
    lsm_cube_data = np.maximum(lsm_pi_cube.data, lsm_plio_cube.data)
    lsm_cube_ = lsm_pi_cube.copy(data=lsm_cube_data)
  
    return lsm_cube_, lsm_plio_cube

def get_data_cube(filename, fieldname, grid_cube):
    """
    get the data from the given file and regrid change missing data to 1000.
    """

    cube = iris.load_cube(filename, fieldname)
    regrid_cube = cube.regrid(grid_cube, iris.analysis.Linear())

    regrid_data = regrid_cube.data
    regrid_data = np.where(regrid_data < 1000., regrid_data, 1000.)
    new_cube = regrid_cube.copy(data=regrid_data)
  
    return new_cube


def get_model_data():
    """
    first read in the lsm
    read in data from the pliocene and the preindustrial and regrid
    if a point is land get the data from the NSAT file
    if a point is ocean get the data from the SST file
    """
    
    (lsm_cube, lsm_plio_cube) = get_lsm()

    nsat_cube = get_data_cube(NSAT_MMM_FILE,
                              'NearSurfaceTemperaturemean_anomaly',
                              lsm_cube)
    sst_cube = get_data_cube(SST_MMM_FILE,'SSTmean_anomaly', lsm_cube)
   
   
    # anom_cube = nsat anomaly over land and sst anomaly over ocean)
    anom_cube = ((nsat_cube * lsm_cube) - 
                 (sst_cube * (lsm_cube - 1.0)))

    # there are still a few points which are using the ocean value and 
    # should be using the land value  Change these.
    # check there are not too many

    anom_data = anom_cube.data
    nsat_data = nsat_cube.data
    anom_new_data = np.where(anom_data > 900, nsat_data, anom_data)
    anom_new_cube = anom_cube.copy(anom_new_data)

    return anom_new_cube, lsm_plio_cube

def get_data():
    """
    this function willl open the file containing the data and will return 
    arrays containing:
        1. site_longitude
        2. site_latitudes
        3. T anomaly from (NOAA-ERSSTv5)
        4. standard deviation of the data
        5. number of sites
    """
    dfs = pd.read_excel(DATAFILE)
    dfs_subset = dfs[["Latitude (¡N)", "Longitude (¡E)", "NOAA_anom", "Standard dev.", "N"]]
    
    lats = dfs_subset.iloc[:,0]
    lons = dfs_subset.iloc[:,1]
    data_tanom = dfs_subset.iloc[:,2]
    data_stdev = dfs_subset.iloc[:,3]
    Npoints = dfs_subset.iloc[:,4]
    
    
    return lats, lons, data_tanom, data_stdev, Npoints

class GetPliovar:
    """
    this class is to do with getting everything from Heathers excel files
    """
    def __init__(self, interval, datatype):
        """
        the interval is esentially which excel sheet we are getting data from
        t1 t2 or t3
        datatype = UK37 or MGCA
        """
        
        if datatype == 'UK37':
            self.filename = DATASTART + 'pliovar_uk37_ori_vs_bayspline.xlsx'
            self.bsloc = 8
        if datatype == 'MGCA':
            self.filename = DATASTART +  'pliovar_mgca_OrivsBaymag.xlsx'
            self.bsloc = 7
        self.metafile = DATASTART + 'pliovar_metadata_global_02102019.csv'
        self.pifile = DATASTART + 'modeloutput_pliovar.xls'
        self.interval = interval # this is the time range likely t1 t2 or t3
           
    def get_proxydata(self):
        """
        this will obtain in an array the latitude, longitude and SST of the 
        proxy data.  It will put them in an array
        
        returns for each latitude bound
        boundtemp : the average temperature in the latitude band
        boundtemp_bs : the average temperature in the latitude band using bayspline
        boundmin ; the minimum latitude of the band
        boundmax : the maximum latitude of the band
        nval: the number of points in the band (for weighting)
        """
        
        # reads into a dictionary
        dfs = pd.read_excel(self.filename, sheet_name=None)
        
        t1sheet = dfs.get(self.interval)
        

        self.sitenames = t1sheet.iloc[1:,0]
        self.nsites = len(self.sitenames)
        self.lon = np.zeros(self.nsites)
        self.lat = np.zeros(self.nsites)
        self.temppi = np.zeros(self.nsites)
        
        
        # get the temperatures
        self.sitetemp = t1sheet.iloc[1:,1]
        self.sitetemp_bs = t1sheet.iloc[1:,self.bsloc]
        
        
        # get the latitudes and longitudes
        self.get_lonlat() 
        
        # get the preindustrial temperatures
        self.get_piT() 
        
        if HARRY_ERIN == 'Eb':
            data_tanom = self.sitetemp_bs - self.temppi
        if HARRY_ERIN == 'En':
            data_tanom = self.sitetemp - self.temppi

        latuse = []
        lonuse = []
        tanom_use = []
        nsites_use = 0
        for i, tanom in enumerate(data_tanom):
            if np.isfinite(tanom):
                latuse.append(self.lat[i])
                lonuse.append(self.lon[i])
                tanom_use.append(tanom)
                nsites_use = nsites_use + 1

        
        return latuse, lonuse, tanom_use, nsites_use
       
    def get_lonlat(self):
        """
        will get the longitude and laitude from each site
        and add them to the self.lon and self.lat array
        """
        
        # gets the dictionary of longitudes and latitudes
        # from the metadatafile
        df = pd.read_csv(self.metafile, encoding='latin-1')
        metadf = df[["name", "lon", "lat"]]
        lonlatdict = metadf.set_index('name').T.to_dict()
        
        #print(lonlatdict)
        #sys.exit(0)
        
        for i in range(0, self.nsites):
            sitedata = lonlatdict.get(self.sitenames.iloc[i],'lat')
            self.lat[i] = sitedata.get('lat')
            self.lon[i] = sitedata.get('lon')
            
        return
    
 
    def get_piT(self):
        """
        will get the pi temperature from each site from NOAASST
        and add to self.pitemp array
        """
        
        dfs = pd.read_excel(self.pifile, sheet_name='E280near')
        # gets the dictionary of longitudes and latitudes
        # from the metadatafile
        metadf = dfs[["site", "NOAAERSST5"]]
       
        pitempdict = metadf.set_index(['site']).T.to_dict()
        
        
        for i in range(0, self.nsites):
            noaadata = pitempdict.get((self.sitenames.iloc[i]))
            self.temppi[i] = noaadata.get('NOAAERSST5')
           
        return
   
def shift_lons(lons,lats,temp):
    """ 
    if two points are in the same location then shift longitude slightly so that both are 
    visible
    """

    new_lons =  np.zeros(np.shape(lons))
    new_lons[:] = lons[:]

    for i, lon in enumerate(lons):
        subscript_same = []
        for j in range(i+1, len(lons)):
            if (np.abs(lon - lons[j]) < 1.0 and np.abs(lats[i] - lats[j]) < 1.0):
                subscript_same.append(j)
                print(i,j)
        for s, subscript in enumerate(subscript_same):
            if lons[subscript] == new_lons[subscript]:
                new_lons[i] = lons[i] - 2.0
                print('here',i,new_lons[i],lons[i],lons[i]-0.5,s)
                new_lons[subscript] = lons[subscript] + 2.0 + (4.0 *s)
            
    print(lons[0],new_lons[0], new_lons[14], new_lons[15])
    return new_lons

def plot(model_cube, mask_cube, lats, lons, data):

    """
    plots the model anomaly with the data anomaly on top
    """

   
    # plot model
    vmin = -10.0
    vmax = 10.0
    incr = 1.0
    V = np.arange(vmin, vmax + incr, incr)
    mycmap = customise_cmap2()

    #brewer_cmap = cm.get_cmap('brewer_RdBu_11')
    ax = plt.axes(projection=ccrs.Robinson(central_longitude=0))
    cs = iplt.contourf(model_cube, levels=V,  extend='both',
                       cmap=mycmap)
    cbar = plt.colorbar(cs,  orientation= 'horizontal',
                        ticks=[-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10])
    cbar.set_label('deg C')
    iplt.contour(mask_cube, colors='black', linewidths=0.1)
    plt.title('3.205Ma - PI temperature anomaly')
    

    # overplot data
  
    #norm = colors.Normalize(vmin = vmin, vmax = vmax)
    norm = colors.BoundaryNorm(boundaries=V, ncolors=mycmap.N)
    print(norm)

    plt.scatter(lons, lats, c='black',  marker='o', s=60, transform=ccrs.Geodetic())

    plt.scatter(lons, lats, c=data,  marker='o', s=45,
                norm = norm , cmap=mycmap, transform=ccrs.Geodetic())

  
    plt.savefig('/nfs/hera1/earjcti/regridded/IPCC_Tanom_' + OUTSS + '.png')
    plt.savefig('/nfs/hera1/earjcti/regridded/IPCC_Tanom_' + OUTSS + '.eps')


  
def main():
    """
    calling structure
    a) get's model data
    b) get's proxy data
    c) plots model data with proxy data on top
    """

    model_anom_cube, lsmplio_cube = get_model_data()

    if HARRY_ERIN == 'H':
        lats, lons, data_tanom, data_stdev, Npoints = get_data()
    if HARRY_ERIN == 'En' or HARRY_ERIN == 'Eb':
        obj = GetPliovar('t1', 'MGCA') # get data for t1 timeslice
        lats, lons, data_tanom, Npoints = obj.get_proxydata() 
        obj = GetPliovar('t1', 'UK37') # get data for t1 timeslice
        lats_UK37, lons_UK37, data_tanom_UK37, Npoints_UK37 = obj.get_proxydata() 
        for i in range(0, Npoints_UK37):
            lats.append(lats_UK37[i])
            lons.append(lons_UK37[i])
            data_tanom.append(data_tanom_UK37[i])

    # if two points are same shift them so they are both visible
    lons_shift = shift_lons(lons, lats, data_tanom)
        
    plot(model_anom_cube, lsmplio_cube, lats, lons_shift, data_tanom)

##########################################################
# main program


LINUX_WIN = 'l'
FILESTART = '/nfs/hera1/earjcti/'
HARRY_ERIN = 'H' # H=Harry, En Erin Normal, Eb Erin Bayspline



LSM_PLIO_FILE = (FILESTART + 'regridded/PlioMIP2_Boundary_conds/Plio_enh' 
            + '/Plio_enh/Plio_enh_LSM_v1.0.nc')
LSM_PI_FILE = (FILESTART + 'regridded/PlioMIP2_Boundary_conds/Modern_std' 
            + '/Modern_std/Modern_std_LSM_v1.0.nc')
NSAT_MMM_FILE = (FILESTART + 
                 'regridded/NearSurfaceTemperature_multimodelmean.nc')
SST_MMM_FILE = (FILESTART + 'regridded/SST_multimodelmean.nc')


if HARRY_ERIN == 'H':
    DATAFILE = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/cs_mp_sst_data_30k_plusNOAA.xlsx'
    OUTSS = 'FD30'

if HARRY_ERIN == 'Eb':
    DATASTART = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/'
    OUTSS = 'McClymont_Bayspline'

if HARRY_ERIN == 'En':
    DATASTART = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/'
    OUTSS = 'McClymont_Standard'
    
main()

#sys.exit(0)
