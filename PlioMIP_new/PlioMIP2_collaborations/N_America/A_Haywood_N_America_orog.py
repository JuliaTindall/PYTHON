#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#Created October 2021

#@author: earjcti
#
# This program plot a figure because of the menemenlis paper
# about hydroclimate in the southwest US

# a) MPWP - PI orography  PlioMIP2 (MMM)
# b) MPWP - PI SAT orography  PlioMIP1 (MMM)
# c) b-a
# d) MPWP - PI precip anomaly  PlioMIP2 (MMM)
# e) MPWP - PI precip anomaly  PlioMIP1 (MMM)
# f) e-d


#import os
import numpy as np
import pandas as pd
#import scipy as sp
#import cf
import iris
#import iris.util
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
import iris.plot as iplt
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

def customise_cmap3():
    """
    as customise_cmap but 19 colors only + 2 white in middle added by Julia
    """
    colors = [(84, 48, 5), (113, 70, 16), (143, 93, 27), (173, 115, 38),
              (195, 137, 60), (206, 160, 97), (216, 182, 135),
              (227, 204, 173), (238, 226, 211), (248, 248, 247),
              (212, 230, 229), (176, 212, 209), (140, 194, 190),
              (103, 176, 170), (67, 158, 150), (44, 135, 127),
              (29, 110, 100), (14, 85, 74), (0, 60, 48)]
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


   

def average_cubes(cubelist, name, units):
    """
    finds the average of the cube
    """
    newcubelist = iris.cube.CubeList([])

    for i, cube in enumerate(cubelist):
        # makes cube ready for concatenation
        print(cube.var_name)
  
        cube.long_name = name
        cube.short_name = name
        cube.var_name = name
        cube.standard_name = None
        cube.attributes = None
        cube.units = units
        cube = iris.util.squeeze(cube)
        cube.cell_methods=None
        cube.data = cube.data.astype('float32')
        for coord in cube.coords():
            coodname = coord.standard_name
            coord.points = coord.points.astype('float32')
            print(coord)
            if coodname !='latitude' and coodname !='longitude':
                if coodname==None:
                    if coord.long_name==None:
                        cube.remove_coord(coord.var_name)
                    else:
                        cube.remove_coord(coord.long_name)
                else:
                    cube.remove_coord(coodname)
        # add new axis for merging
        tempcube=iris.util.new_axis(cube)
        tempcube.add_dim_coord(iris.coords.DimCoord(i, 
            standard_name='model_level_number', long_name='model', 
            var_name='model', 
            units=None,
            bounds=None,
            coord_system=None, circular=False),0) 
    
        newcubelist.append(tempcube)




    iris.experimental.equalise_cubes.equalise_attributes(newcubelist)
    allcubes = newcubelist.concatenate_cube()
    cubeavg = allcubes.collapsed(['model_level_number'], iris.analysis.MEAN)
   
    return cubeavg

def get_model_p1():
    """
    first read in the lsm
    read in data from the pliocene and the preindustrial and regrid
    if a point is land get the data from the NSAT file
    if a point is ocean get the data from the SST file
    """
    
 
    modelnames = ['COSMOS', 'Had', 'CCSM', 'IPSL', 'MIROC', 'MRI', 'NOR']
    filename = '/nfs/hera1/earjcti/PLIOMIP/PlioMIP1_regridded.nc'
    cubeall = iris.load(filename)
    plio_sat_cubes = iris.cube.CubeList([])
    pi_sat_cubes = iris.cube.CubeList([])
    plio_sst_cubes = iris.cube.CubeList([])
    pi_sst_cubes = iris.cube.CubeList([])
    plio_precip_cubes = iris.cube.CubeList([])
    pi_precip_cubes = iris.cube.CubeList([])
       
    for model in modelnames:
        
        sat_pi = model + '_ctrl_sat'
        sat_plio = model + '_plio_sat'
        sst_pi = model + '_ctrl_sst'
        sst_plio = model + '_plio_sst'
        precip_pi = model + '_ctrl_precip'
        precip_plio = model + '_plio_precip'
    
        for cubetemp in cubeall:
            var = cubetemp.var_name
            if sat_pi.lower() in var.lower():
                pi_sat_cubes.append(cubetemp)
            if sat_plio.lower() in var.lower():
                plio_sat_cubes.append(cubetemp)
            if sst_pi.lower() in var.lower():
                pi_sst_cubes.append(cubetemp)
            if sst_plio.lower() in var.lower():
                plio_sst_cubes.append(cubetemp)
            if precip_pi.lower() in var.lower():
                pi_precip_cubes.append(cubetemp)
            if precip_plio.lower() in var.lower():
                plio_precip_cubes.append(cubetemp)



    plio_avgsatcube = average_cubes(plio_sat_cubes, 'sat', 'degC')
    pi_avgsatcube = average_cubes(pi_sat_cubes, 'sat', 'degC')   
    plio_avgsstcube = average_cubes(plio_sst_cubes, 'sst', 'degC')
    pi_avgsstcube = average_cubes(pi_sst_cubes, 'sat', 'degC')   
    plio_avgprecipcube = average_cubes(plio_precip_cubes, 'precip', 'mm')
    pi_avgprecipcube = average_cubes(pi_precip_cubes, 'precip', 'mm')   
    anom_sat_cube = plio_avgsatcube - pi_avgsatcube
    anom_sst_cube = plio_avgsstcube - pi_avgsstcube
    anom_precip_cube = plio_avgprecipcube - pi_avgprecipcube
   
    # get land sea mask

    lsminit = iris.load_cube('/nfs/hera1/earjcti/PRISM/PLIOMIP/exp2_preferred/land_fraction_v1.1.nc')
    lsm_newdata = np.where(lsminit.data > 0.0, 1.0, 0.0)
    lsm_cube = lsminit.copy(data=lsm_newdata)

    Tanom_data =  ((anom_sat_cube.data * lsm_cube.data) - 
                  (anom_sst_cube.data * (lsm_cube.data - 1.0)))
    Tanom_data.mask = False
    Tanom_cube = anom_sat_cube.copy(data=Tanom_data)

        
    return Tanom_cube, anom_precip_cube, lsm_cube

 
def get_land_obs():
    """
    reads in the spredsheet from ulrich and returns temperatures
    """

    dfs = pd.read_excel(LAND_DATAFILE)
    sites = []
    lats = []
    lons = []
    temps = []
    temp_modern = []

    row_locs = [2, 3, 4, 5, 6, 7, 8, 9, 11, 12]
    for rl in row_locs:
        # if temp ne nan then move to array
        temp = dfs.iloc[rl, 9]
        print(temp,'julia')
        if np.isfinite(temp):
            sites.append(dfs.iloc[rl, 0])
            lats.append(dfs.iloc[rl, 2])
            lons.append(dfs.iloc[rl, 3])
            temp_modern.append(dfs.iloc[rl, 4])
            temps.append(temp)

    return lats, lons, temps, temp_modern



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
            self.filename = P2_DATASTART + 'pliovar_uk37_ori_vs_bayspline.xlsx'
            self.bsloc = 8
        if datatype == 'MGCA':
            self.filename = P2_DATASTART +  'pliovar_mgca_OrivsBaymag.xlsx'
            self.bsloc = 7
        self.metafile = P2_DATASTART + 'pliovar_metadata_global_02102019.csv'
        self.pifile = P2_DATASTART + 'modeloutput_pliovar.xls'
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
        
        data_tanom = self.sitetemp_bs - self.temppi
        
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


  
def main_orog():
    """
    calling structure
    a) get's model data
    """

    filename = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_topo_v1.0.nc'
    p2_orog_cube_temp = iris.load_cube(filename,'p4_topo')
    # if data less than 0 set to 0, we are not interested in bathymetry
    data = p2_orog_cube_temp.data
    newdata = np.where(data > -0, data, 0)
    p2_orog_cube = p2_orog_cube_temp.copy(data=newdata)

    filename = '/nfs/hera1/earjcti/PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_topo_v1.0.nc'
    pi_orog_cube_temp = iris.load_cube(filename,'etopo1_topo')
    # if data less than 0 set to 0, we are not interested in bathymetry
    data = pi_orog_cube_temp.data
    newdata = np.where(data > 0, data, 0)
    pi_orog_cube = pi_orog_cube_temp.copy(data=newdata)


    # p1 put on same grid
    filename = '/nfs/hera1/earjcti/PRISM/PLIOMIP/exp2_preferred/topo_v1.1.nc'
    #p1_orog = iris.load(filename)
    p1_orog_cube_temp = iris.load_cube(filename,'height')
    p1_orog_cube = p1_orog_cube_temp.regrid(p2_orog_cube, iris.analysis.Linear())



    return [p2_orog_cube, p1_orog_cube, pi_orog_cube]


def plot_orog(p2_orog_cube, p1_orog_cube, pi_orog_cube):
    """
    plots the orography
    """

    # P2
    V = np.arange(0,2500,250)
    ax1=plt.subplot(2,3,1,projection=ccrs.PlateCarree())
    ax1.set_extent([-150, -60, 10, 70])
    cs = iplt.contourf(p2_orog_cube, V, extend='both')
    cbar = plt.colorbar(cs, orientation='horizontal')
    cbar.ax.tick_params(labelsize=7, rotation=90)
    plt.gca().coastlines()
    plt.title('PlioMIP2 orography (m)',fontsize=10)
    
    #P1
    ax2=plt.subplot(2,3,2,projection=ccrs.PlateCarree())
    ax2.set_extent([-150, -60, 10, 70])
    cs = iplt.contourf(p1_orog_cube, V, extend='both')
    cbar = plt.colorbar(cs, orientation='horizontal')
    cbar.ax.tick_params(labelsize=7, rotation=90)
    plt.gca().coastlines()
    plt.title('PlioMIP1 orography (m)',fontsize=10)
  
    #PI
    ax3=plt.subplot(2,3,3,projection=ccrs.PlateCarree())
    ax3.set_extent([-150, -60, 10, 70])
    cs = iplt.contourf(pi_orog_cube, V, extend='both')
    cbar = plt.colorbar(cs, orientation='horizontal')
    cbar.ax.tick_params(labelsize=7, rotation=90)
    plt.gca().coastlines()
    plt.title('PI orography (m)',fontsize=10)

    # P2 - PI
    V = np.arange(-400,440,40)
    ax1=plt.subplot(2,3,4,projection=ccrs.PlateCarree())
    ax1.set_extent([-150, -60, 10, 70])
    cs = iplt.contourf(p2_orog_cube - pi_orog_cube, V, extend='both')
    cbar = plt.colorbar(cs, orientation='horizontal')
    cbar.ax.tick_params(labelsize=7, rotation=90)
    plt.gca().coastlines()
    plt.title('Plio2 - PreInd orog (m)',fontsize=10)
    
    #P1 - PI
    ax2=plt.subplot(2,3,5,projection=ccrs.PlateCarree())
    ax2.set_extent([-150, -60, 10, 70])
    cs = iplt.contourf(p1_orog_cube - pi_orog_cube, V, extend='both')
    cbar = plt.colorbar(cs, orientation='horizontal')
    cbar.ax.tick_params(labelsize=7, rotation=90)
    plt.gca().coastlines()
    plt.title('Plio1 - PreInd orog (m)',fontsize=10)
  
    #P2 -P1
    ax3=plt.subplot(2,3,6,projection=ccrs.PlateCarree())
    ax3.set_extent([-150, -60, 10, 70])
    cs = iplt.contourf(p2_orog_cube - p1_orog_cube, V, extend='both')
    cbar = plt.colorbar(cs, orientation='horizontal')
    cbar.ax.tick_params(labelsize=7, rotation=90)
    plt.gca().coastlines()
    plt.title('Plio2 - Plio1 orog (m)',fontsize=10)
  
    plt.savefig('orogdiff.png')
    
    


def get_precip():

    precfile = '/nfs/hera1/earjcti/regridded/TotalPrecipitation_multimodelmean.nc'
    precip_cube = iris.load_cube(precfile, 'TotalPrecipitationmean_anomaly')
   
   
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

    return anom_new_cube, lsm_plio_cube, precip_cube

   
        

    return [model_anom_cube, lsmplio_cube, lats, lons_shift, data_tanom, 
            land_lats, land_lons, land_tanom, model_anom_precip_cube]


def get_p1_land():
    """
    gets the land data from p1
    """

    lons = []
    lats = []
    temps = []

    filename = '/nfs/hera1/earjcti/PRISM/prism3_pliocene/SAT_2012_07.txt'
    f1 = open(filename)
    lines = f1.readlines()
    for i in range(1, len(lines)):
        line = lines[i]
        lon, lat, temp, dummy = line.split()
        lons.append(np.float(lon))
        lats.append(np.float(lat))
        temps.append(np.float(temp))

    return lons, lats, temps

def get_p1_ocn():
    """
    gets the ocean data from p1
    """

    lons = []
    lats = []
    temps = []

    filename = '/nfs/hera1/earjcti/PRISM/prism3_pliocene/SST_2012_07_reformat.txt'
    f1 = open(filename)
    lines = f1.readlines()
    for i in range(1, len(lines)):
        line = lines[i]
        print(line)
        print(line.split())
        name, lat, lon, modern, prism, anomaly = line.split()
        lons.append(np.float(lon))
        lats.append(np.float(lat))
        temps.append(np.float(anomaly))

    return lons, lats, temps

def main_p1_sat():
    """
    gets the pliomip1 data
    calling structure
    a) get's model data
    b) get's proxy data
    """

    # get model data
    (model_satanom_cube, model_precipanom_cube, lsm_cube) = get_model_p1()
   
    # get data from p1
    (land_lon, land_lat, land_data) = get_p1_land()
    (ocn_lon, ocn_lat, ocn_data) = get_p1_ocn()
  
   
    
    return [model_satanom_cube, lsm_cube, land_lon, land_lat, land_data,
            ocn_lon, ocn_lat, ocn_data, model_precipanom_cube]

########################################################

def plot(ax, i, model_cube, mask_cube, lats, lons, data, land_lats, land_lons,
         land_data, plot_ocean, plot_land):

    """
    plots the model anomaly with the data anomaly on top
    """
   
    plt.subplot(2,3,i)
    # plot model
    if i ==1 or i==2:  # SAT plot
        vmin = -10.0
        vmax = 10.0
        incr = 1.0
        ticks = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10]
        mycmap = customise_cmap2()
        if i ==1:
            title = 'a) PlioMIP2 '
        if i==2:
            title = 'b) PlioMIP1'
        units = 'deg C'


    if i ==3: #SAT anom
        vmin = -3.0
        vmax = 3.0
        incr = 0.5
       # ticks = [-2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5]
        ticks = [-3, -2, -1, 0,  1,  2, 3]
        mycmap = customise_cmap2()
        title = 'c) PlioMIP2 - PlioMIP1'
        units = 'deg C'



    if i ==4 or i==5:  # precip plot
        vmin = -2.0
        vmax = 2.0
        incr = 0.2
        ticks = [-1.8, -1.2, -0.6, 0, 0.6,  1.2, 1.8]
        mycmap = customise_cmap3()
        if i ==4:
            title = 'd) PlioMIP2 '
        if i==5:
            title = 'e) PlioMIP1'
        units = 'mm / day'



    if i ==6:  # precip anom
        vmin = -2.0
        vmax = 2.0
        incr = 0.2
        ticks = [-1.8, -1.2, -0.6, 0, 0.6,  1.2, 1.8]
        mycmap = customise_cmap3()
        title = 'f) PlioMIP2 - PlioMIP1'
        units = 'mm / day'




    V = np.arange(vmin, vmax + incr, incr)
   
    # turn the iris cube data structure into numpy arrays
    gridlons_d = model_cube.coord('longitude').points
    gridlats_d = model_cube.coord('latitude').points
    moddata = model_cube.data


    print(model_cube)
    print(moddata.shape, gridlons_d.shape,  gridlats_d.shape)
    #cs1 = ax.contourf(gridlons_d, gridlats_d, moddata, levels=V,  extend='both',
    #                   cmap=mycmap)
   # 
    #plt.title('3.205Ma - PI temperature anomaly')
 
    cs = iplt.contourf(model_cube, V, extend='both', cmap=mycmap)
    print('j1',gridlons_d.shape, gridlats_d.shape, mask_cube.data.shape)
    plt.contour(mask_cube.coord('longitude').points, 
                mask_cube.coord('latitude').points,mask_cube.data, colors='black', linewidths=0.1)
  
    plt.title(title, fontsize=9)
    cbar =plt.colorbar(cs,  orientation= 'horizontal', ticks = ticks)
    cbar.set_label(units, fontsize=7)
    cbar.ax.tick_params(labelsize=7) 
 

    # overplot data ocean
  
    #norm = colors.Normalize(vmin = vmin, vmax = vmax)
    norm = colors.BoundaryNorm(boundaries=V, ncolors=mycmap.N)
    print(norm)

    if plot_ocean:
        plt.scatter(lons, lats, c='black',  marker='o', s=30, transform=ccrs.Geodetic())

        plt.scatter(lons, lats, c=data,  marker='o', s=15,
                norm = norm , cmap=mycmap, transform=ccrs.Geodetic())


    # overplot data lane
 
    if plot_land:
        plt.scatter(land_lons, land_lats, c='black',  
                    marker='o', s=30, transform=ccrs.Geodetic())#

        plt.scatter(land_lons, land_lats, c=land_data,  marker='o', s=15,
                    norm = norm , cmap=mycmap, transform=ccrs.Geodetic())
  
  
###################################################################

def main():
    """
    calling structure
    a) get's orog (model and data) from pliomip2
    b) get's orog (model and data) from pliomip1
    c) get's total precip (model and data) from pliomip2
    d) get's total precip (model and data) from pliomip1
    d) plot
    """
    
  #  fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(7.09, 6.0),
#                      subplot_kw={'projection': ccrs.Robinson(central_longitude=0)})
  #  fig.set_dpi(300.0)

  #  fig.text(0.02, 0.86, 'Temperature anomaly:', fontsize = 12)
   
  #  fig.text(0.02, 0.41, 'Precipitation anomaly:', fontsize = 12)
   

    # process pliomip2 temperature
    (p2_orog_cube, p1_orog_cube, pi_orog_cube) = main_orog()
   

    # plot_orography
    plot_orog(p2_orog_cube, p1_orog_cube, pi_orog_cube)

    plot(ax[0,0], 1, model_satanom_p2_cube, lsmplio_p2_cube, lats_p2, lons_shift_p2,
         p2_data_tanom, p2_land_lats, p2_land_lons, p2_land_tanom,True,True)


    # process pliomip1 temperature
    (model_satanom_p1_cube, lsmplio_p1_cube, land_lat_p1, 
     land_lon_p1, land_data_p1, ocean_lon_p1, ocean_lat_p1, sst_p1,
     model_precipanom_p1_cube) = main_p1_sat()
    dummy=-999.
    plot(ax[0,1], 2, model_satanom_p1_cube, lsmplio_p1_cube, ocean_lat_p1, ocean_lon_p1,
         sst_p1, land_lat_p1, land_lon_p1, land_data_p1, True, True)

  
    # plot p2 - p1
    p1_regrid = model_satanom_p1_cube.regrid(model_satanom_p2_cube, iris.analysis.Linear())
    p2_p1_anom = model_satanom_p2_cube.data - p1_regrid.data
    p2_p1_cube = model_satanom_p2_cube.copy(data=p2_p1_anom)

    plot(ax[0,2],3, p2_p1_cube, lsmplio_p2_cube, dummy, dummy, dummy, dummy, dummy, dummy,
         False, False)



    # plot precipitation
    plot(ax[1,0],4, model_precipanom_p2_cube, lsmplio_p2_cube,
         dummy, dummy, dummy, dummy, dummy, dummy, False, False)
  
    plot(ax[1,1],5, model_precipanom_p1_cube, lsmplio_p1_cube,
         dummy, dummy, dummy, dummy, dummy, dummy, False, False)

    # plot p2 - p1
    p1_regrid = model_precipanom_p1_cube.regrid(model_precipanom_p2_cube, iris.analysis.Linear())
    p2_p1_anom = model_precipanom_p2_cube.data - p1_regrid.data
    p2_p1_cube = model_precipanom_p2_cube.copy(data=p2_p1_anom)

    plot(ax[1,2],6, p2_p1_cube, lsmplio_p2_cube, dummy, dummy, dummy, dummy, dummy, dummy,
         False, False)

    #plt.tight_layout(w_pad=3)
    plt.tight_layout()
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/A_Haywood_PAGES.png')
    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/A_Haywood_PAGES.eps')
    plt.close()
   

##########################################################
# main program


LINUX_WIN = 'l'
FILESTART = '/nfs/hera1/earjcti/'

OROG_P2_FILE = (FILESTART + 'regridded/PlioMIP2_Boundary_conds/Plio_enh' 
            + '/Plio_enh/Plio_enh_topo_v1.0.nc')
OROG_PI_FILE = (FILESTART + 'regridded/PlioMIP2_Boundary_conds/Modern_std' 
            + '/Modern_std/Modern_std_LSM_v1.0.nc')
NSAT_MMM_FILE = (FILESTART + 
                 'regridded/NearSurfaceTemperature_multimodelmean.nc')
SST_MMM_FILE = (FILESTART + 'regridded/SST_multimodelmean.nc')


P2_DATASTART = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/'
P2_OUTSS = 'McClymont_Bayspline'
    
LAND_DATAFILE = ('/nfs/hera1/earjcti/PLIOMIP2/proxydata/' + 
                 'PlioceneTerrestrial_IPCCAR6.xlsx')

main()

#sys.exit(0)
