#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created February 2021 by Julia

This program will plot the biomes in a nice way

; tropical evergreen broadleaf forest                = 1
; tropical semi-evergreen broadleaf forest           = 2
; tropical deciduous broadleaf forest & woodland     = 3
; temperate deciduous broadleaf forest               = 4
; temperate evergreen needleleaf forest              = 5
; warm-temperate evergreen broadleaf & mixed forest  = 6
; cool mixed forest                                  = 7
; cool evergreen needleleaf forest                   = 8
; cool-temperate evergreen needleleaf & mixed forest = 9
; cold evergreen needleleaf forest                   = 10
; cold deciduous forest                              = 11
; tropical savanna                                   = 12
; tropical xerophytic shrubland                      = 13
; temperate xerophytic shrubland                     = 14
; temperate sclerophyll woodland and shrubland       = 15
; temperate deciduous broadleaf savanna              = 16
; temperate evergreen needleleaf open woodland       = 17
; cold parkland                                      = 18
; tropical grassland                                 = 19
; temperate grassland                                = 20
; desert                                             = 21
; graminoid and forb tundra                          = 22
; low and high shrub tundra                          = 23
; erect dwarf-shrub tundra                           = 24
; prostrate dwarf-shrub tundra                       = 25
; cushion-forb tundra                                = 26
; barren                                             = 27
; ice                                                = 28


"""
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
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

def get_biome_details():
    """
    gets the biome names and colormap.  These are hard coded
    """
    
    names = ["Tropical evergreen forest",    
             "Tropical semi-deciduous forest",
             "Tropical decididous forest",
             "Temperate deciduous forest",
             "Temperate conifer forest",
             "warm mixed forest",
             "cool mixed forest",
             "cool conifer forest",  
             "cool-temperate mixed forest",
             "Evergreen taiga/montane forest",                  
             "Deciduous taiga/montane forest",                           
             "Tropical savanna",                            
             "Tropical xerophytic shrubland",                  
             "Temperate xerophytic shrubland",                   
             "Temperate sclerophyll woodland",     
             "Temperate broadleaf savanna",              
             "Open conifer woodland",     
             "Boreal parkland",                                  
             "Tropical grassland",                              
             "Temperate grassland",                              
             "Desert",                                           
             "Steppe tundra",                         
             "Shrub tundra",                         
             "Dwarf-shrub tundra",                          
             "Prostrate shrub tundra",                      
             "Cushion-forb tundra",                               
             "Barren",                                           
             "Land ice"] 

    colors = [ #( 1.000, 1.000, 1.000), 
            # ( 0.000, 0.000, 0.000 ),
             ( 0.110, 0.333, 0.063 ),
             ( 0.396, 0.573, 0.031 ),
             ( 0.682, 0.490, 0.125 ),
             ( 0.333, 0.922, 0.286 ),
             ( 0.094, 0.510, 0.443 ),
             ( 0.000, 0.000, 0.396 ),
             ( 0.792, 1.000, 0.792 ),
             ( 0.000, 0.604, 0.094 ),
             ( 0.443, 0.141, 0.208 ),
             ( 0.000, 0.125, 0.792 ),
             ( 0.396, 0.698, 1.000 ),
             ( 0.729, 1.000, 0.208 ),
             ( 1.000, 0.729, 0.604 ),
             ( 1.000, 0.875, 0.792 ),
             ( 0.557, 0.635, 0.157 ),
             ( 0.459, 1.000, 0.208 ),
             ( 1.000, 0.604, 0.875 ),
             ( 0.396, 0.490, 1.000 ),
             ( 1.000, 0.729, 0.208 ),
             ( 1.000, 0.875, 0.604 ),
             ( 0.969, 1.000, 0.792 ),
             ( 0.906, 0.906, 0.094 ),
             ( 0.396, 1.000, 0.604 ),
             ( 0.475, 0.525, 0.286 ),
             ( 0.824, 0.620, 0.588 ),
             ( 0.604, 0.396, 1.000 ),
             ( 0.729, 0.714, 0.667 ),
             ( 0.714, 0.824, 0.875 ) 
             # this last one was commented out
             #,( 0.700, 0.700, 0.700 )
             ]

    #cmap = mpl.colors.LinearSegmentedColormap.from_list(
    #    'Custom cmap', cmaplist, cmap.N)
    
   # cmap = mpl.colors.ListedColormap(['red',    'green',  'blue', 
   #                                   'cyan', 'red',   'green',  'blue', 
   #                                   'cyan',   'red',    'green',    'blue', 
   #                                   'cyan', 'red',  'green',  'blue', 
   #                                   'cyan',    'red',  'green', 'blue', 
   #                                   'cyan',   'red',  'green',  'blue', 
   #                                   'cyan',  'red',  'green',  'blue', 
   #                                   'cyan'    ])
    cmap = mpl.colors.ListedColormap(colors)
    bounds = np.linspace(0.5, 27.5, 28)
    print(bounds)
    #sys.exit(0)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    return names, colors, cmap, norm, bounds

def get_data(colors):
    """
    reads in the data and puts it in a array (M, N, 3) with the
    third dimensio being the colors
    """
    cube = iris.load_cube(FILENAME, 'biome')
    cubedata = cube.data
    ny, nx = np.shape(cubedata)
    colordata = np.zeros((ny, nx, 3))

    print(np.shape(cubedata))
    print(np.shape(colordata))

    for index, point in np.ndenumerate(cubedata):
        if point >= 0:
            colordata[index[0], index[1], :] = colors[point-1]
        else:
            colordata[index[0], index[1], :] = (1, 1, 1)

    return cube, colordata
      
def plot_biomes(names, rgbdata, biomecube, cmap, norm, V,):
    """
    plots the biomes
    """
    fig = plt.figure(figsize=(11.0, 11.0))
    ax = plt.axes(projection = ccrs.PlateCarree())
    #ax.set_extent([-180, 180, -65, 90])
    ax.coastlines()

    # turn the iris Cube data structure into numpy arrays
    gridlons = biomecube.coord('longitude').contiguous_bounds()
    gridlats = biomecube.coord('latitude').contiguous_bounds()
    biomedata = biomecube.data

    print(gridlons)
    print(gridlats)
    print(np.shape(biomedata))
    cs = plt.pcolormesh(gridlons, gridlats, biomedata,
                   cmap=cmap, norm=norm)
    plt.title(MODELNAME + ': ' + ABS_ANOM)

    cbar = plt.colorbar(cs, ticks=V+0.5, )
    cbar.ax.set_yticklabels(names)
    cbar.ax.tick_params (labelsize=6)
    cbar.ax.invert_yaxis()
    #ax.imshow(rgbdata, transform=ccrs.PlateCarree())
    if MODELNAME == 'BIOME4':
        midfname = MODELNAME + '/'
    elif MODELNAME == 'BIOME4_NEWPARAMETERS':
        midfname = 'BIOME4/' + MODELNAME + '/'
    else:
        midfname =  MODELNAME + '/biome4/'
    fileout = (FILESTART + midfname + '_biome4out_' + ABS_ANOM + '.eps')
    print(fileout)
    plt.savefig(fileout)
    print(fileout)
    fileout = (FILESTART + midfname + '_biome4out_' + ABS_ANOM + '.png')
    plt.savefig(fileout)
    sys.exit(0)

def get_prism3_biomes():
    """
    get the biomes from prism3 as a cube
    """
    fname = ('/nfs/hera1/earjcti/PRISM/prism3_pliocene/exp2_alternate/' + 
             'biome_veg_v1.2.nc')
    cube = iris.load_cube(fname)

    # change zeros to nan
    cubedata = cube.data
    cubedata2 = np.where(cubedata != 0, cubedata, np.nan)
    cube_prism3 = cube.copy(data = cubedata2)
   
    fname = ('/nfs/hera1/earjcti/PRISM/prism3_pliocene/BAS_Observ_BIOME.nc')
    cube = iris.load_cube(fname)

    # change zeros to nan
    cubedata = cube.data
    cubedata2 = np.where(cubedata != 0, cubedata, np.nan)
    cube_pi = cube.copy(data = cubedata2)
    

    return cube_prism3, cube_pi

def plot_biomes_with_obs(names, rgbdata, modelcube, datacube, cmap, norm, V,
                         model_ind):
    """
    plots the biomes and the observations on the same figure
    """
    fig,ax = plt.subplots(ncols=1,nrows=2,figsize=(11,8),
                      subplot_kw={'projection': ccrs.PlateCarree()})
    ax1 = ax[0]
    ax2 = ax[1]

    ax1.set_extent([-180, 180, -65, 90])
    ax1.coastlines()

    # do observations
    # turn the iris cube data structure into numpy arrays
    gridlons_d = datacube.coord('longitude').contiguous_bounds()
    gridlats_d = datacube.coord('latitude').contiguous_bounds()
    data = datacube.data

    cs = ax1.pcolormesh(gridlons_d, gridlats_d, data,
                   cmap=cmap, norm=norm)
    if model_ind == 'y':
        ax1.title.set_text('reconstructed biomes for mPWP')
    else:
        ax1.title.set_text('biomes for preindustrial')
 

    # do biome cube
    ax2.set_extent([-180, 180, -65, 90])
    ax2.coastlines()

    # turn the iris Cube data structure into numpy arrays
    gridlons = modelcube.coord('longitude').contiguous_bounds()
    gridlats = modelcube.coord('latitude').contiguous_bounds()
    modeldata = modelcube.data

    cs = ax2.pcolormesh(gridlons, gridlats, modeldata,
                   cmap=cmap, norm=norm)
    if model_ind == 'y':
        if MODELNAME == 'BIOME4':
            ax2.title.set_text('simulated biomes from PlioMIP2 multimodel mean')
        else:
            ax2.title.set_text(MODELNAME + ': ' + ABS_ANOM)
    else:
        ax2.title.set_text('reconstructed biomes for mPWP')


    # adjust plots so we have room for colorbar
    fig.subplots_adjust(bottom=0.05, top = 0.95, left = 0.05, right = 0.70)

    # do colorbar
    cbar_ax = fig.add_axes([0.75, 0.05, 0.05, 0.9])
    cbar = fig.colorbar(cs, cax= cbar_ax, ticks=V+0.5, )
    cbar.ax.set_yticklabels(names)
    cbar.ax.tick_params (labelsize=9)
    cbar.ax.invert_yaxis()
 
    if model_ind == 'y':
        if MODELNAME == 'BIOME4':
            midfname = MODELNAME + '/'
        else:
            midfname =  MODELNAME + '/biome4/'
        fileout = (FILESTART + midfname + '_biome4out_' + ABS_ANOM + '_data.eps')
        plt.savefig(fileout)
        print(fileout)
        fileout = (FILESTART + midfname + '_biome4out_' + ABS_ANOM + '_data.png')
        plt.savefig(fileout)
    else:
        fileout = (FILESTART +  'BIOME4/' + 'pi_vs_plio_obs.eps')
        plt.savefig(fileout)
        print(fileout)
        fileout = (FILESTART +  'BIOME4/' + 'pi_vs_plio_obs.png')
        plt.savefig(fileout)

def get_temp_ranges(biome_prism3_cube, plio_mmmT_cube, i):
        """
        get's the minimum and maximum temperature in this cube for Jan/Apr/Jul/Oct
        """

        biome_mean_temps = np.ma.masked_where(np.repeat(biome_prism3_cube.data[np.newaxis,: :],12, axis=0) != i, plio_mmmT_cube.data)
        biome_mean_temp_cube = plio_mmmT_cube.copy(data=biome_mean_temps)
        biome_NHtemp_cube = biome_mean_temp_cube.extract(iris.Constraint(latitude = lambda cell: cell >= 0))

        plio_temp = []
           
        # min and maximum temperatures are
        month1= biome_NHtemp_cube.extract(iris.Constraint(time=1))
        plio_temp.append([np.float(month1.collapsed(['latitude','longitude'], iris.analysis.MEAN).data), np.float(month1.collapsed(['latitude','longitude'], iris.analysis.MAX).data)])
            
        month4 = biome_NHtemp_cube.extract(iris.Constraint(time=4))
        plio_temp.append([np.float(month4.collapsed(['latitude','longitude'], iris.analysis.MEAN).data), np.float(month4.collapsed(['latitude','longitude'], iris.analysis.MAX).data)])

        month7= biome_NHtemp_cube.extract(iris.Constraint(time=7))
        plio_temp.append([np.float(month7.collapsed(['latitude','longitude'], iris.analysis.MEAN).data), np.float(month7.collapsed(['latitude','longitude'], iris.analysis.MAX).data)])
            
        month10 = biome_NHtemp_cube.extract(iris.Constraint(time=10))
        plio_temp.append([np.float(month10.collapsed(['latitude','longitude'], iris.analysis.MEAN).data), np.float(month10.collapsed(['latitude','longitude'], iris.analysis.MAX).data)])

        return np.asarray(plio_temp)
           
def check_seasonal_for_each_biome(biome_pi_cube_lowres, biome_prism3_cube_lowres,
                                  biome_names):
    """
    1. for pliocene find out which biomes are present polewards of 50N
    2. for each of these biomes find the temperature of these biomes in the
       NH for the Pliocene and the PI throughout the annual cycle
    3. Extract Jan/Apr/July/Oct temperature
    4.
    """
  
    cubegrid = iris.load_cube('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/CEMAC/PLIOMIP2/one_lev_one_deg.nc')
      
    biome_pi_cube = biome_pi_cube_lowres.regrid(cubegrid,iris.analysis.Nearest())
    biome_prism3_cube = biome_prism3_cube_lowres.regrid(cubegrid,iris.analysis.Nearest())

   
    # 1. find out which biomes are present in pliocene polewards of 50N
    present=np.zeros(28,dtype='bool')
    biome_data = biome_prism3_cube.data
    for j,lat in enumerate(biome_prism3_cube.coord('latitude').points):
        if lat >= 50.0:
            for i, lon in enumerate(biome_prism3_cube.coord('latitude').points):
                if np.isfinite(biome_data[j,i]):
                    print(biome_data[j,i])
                    present[np.int(biome_data[j,i])] = True

    # find mean temperature of each biome
    plio_mmmT_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/NearSurfaceTemperature_multimodelmean_month.nc','NearSurfaceTemperaturemean_plio') # load cube
    pi_mmmT_cube = iris.load_cube('/nfs/hera1/earjcti/regridded/NearSurfaceTemperature_multimodelmean_month.nc','NearSurfaceTemperaturemean_pi') # load cube

    allnames = []
    plio_all_range = []
    pi_all_range = []

    # get the ranges of each temperature for months Jan/Apr/July/Oct 
    # where biomes exist
    for i, reqd in enumerate(present):
        if reqd:
           
            allnames.append(biome_names[i-1])
            plio_temps_range = get_temp_ranges(biome_prism3_cube, plio_mmmT_cube, i)
            pi_temps_range = get_temp_ranges(biome_pi_cube, pi_mmmT_cube,  i)

            # note this is a list of arrays
            # each biome corresponding to names, within this each 'month', within this min and max temp
            plio_all_range.append(plio_temps_range)
            pi_all_range.append(pi_temps_range)
            
 
    plio_range = np.asarray(plio_all_range)
    pi_range = np.asarray(pi_all_range)
 

    # plot
    
    fig1 = plt.figure(figsize=[16.0, 20.0])
  

    nbiomes =  len(allnames)
    yarray = np.arange(0,nbiomes,1)
    print(np.shape(plio_range))
#    plt.hlines(y=yarray[0:8], xmin=plio_range[0:8,0,0], xmax=plio_range[0:8,0,1], color='tab:blue', label='January Plio')
#    plt.hlines(y=yarray[0:8]+0.2, xmin=plio_range[0:8,1,0], xmax=plio_range[0:8,1,1], color='tab:green', label='April Plio')
#    plt.hlines(y=yarray[0:8]+0.4, xmin=plio_range[0:8,2,0], xmax=plio_range[0:8,2,1], color='tab:red', label='July Plio')
   
#    plt.hlines(y=yarray[0:8]+0.1, xmin=pi_range[0:8,0,0], xmax=pi_range[0:8,0,1], color='tab:blue', label='January Pi', linestyle='dashed')
#    plt.hlines(y=yarray[0:8]+0.3, xmin=pi_range[0:8,1,0], xmax=pi_range[0:8,1,1], color='tab:green', label='April Pi', linestyle='dashed')
#    plt.hlines(y=yarray[0:8]+0.5, xmin=pi_range[0:8,2,0], xmax=pi_range[0:8,2,1], color='tab:red', label='July Pi', linestyle='dashed')
 
    plt.scatter(plio_range[0:8, 0, 0],yarray[0:8], color='tab:blue', label='January Plio')
    plt.errorbar(plio_range[0:8, 0, 0],yarray[0:8],xerr = plio_range[0:8, 0, 1], color='tab:blue')
    plt.scatter(plio_range[0:8, 1, 0],yarray[0:8], color='tab:green', label='apr Plio')
    plt.scatter(plio_range[0:8, 2, 0],yarray[0:8], color='tab:red', label='Jul Plio')
    plt.scatter(pi_range[0:8, 0, 0],yarray[0:8], color='tab:blue', marker='^',label='Jan Pi')
    plt.scatter(pi_range[0:8, 1, 0],yarray[0:8], color='tab:green', marker='^', label='apr Pi')
    plt.scatter(pi_range[0:8, 2, 0],yarray[0:8], color='tab:red', marker='^',label='jul pi')
  
    for j in range(0, 8):  # just do first 8 biomes, only interested in forest tyeps
        plt.text(-40, yarray[j], allnames[j])
    plt.legend()
    plt.show()

    
            
           
        
    
    

def main():
    """
    driver to get biomes
    """

    print('j1')
    biome_names, colors, cmap, norm, bounds= get_biome_details()
    print('j2')

    cube_model, rgbdata = get_data(colors)
    print('j3')
 
    plot_biomes(biome_names, rgbdata, cube_model, cmap, norm, bounds)
    #print('j4')

    biome_prism3_cube, biome4_pi_cube = get_prism3_biomes()
    plot_biomes_with_obs(biome_names, rgbdata, cube_model, biome_prism3_cube, cmap, norm, bounds,'y')
    #print('j5')
    plot_biomes_with_obs(biome_names, rgbdata, biome_prism3_cube, biome4_pi_cube, cmap, norm, bounds,'n')

    print('j6')
    # look at each biome in turn and find out the Jan / Apr / July/ Oct 
    # temperatures for the PI and the Pliocene
    check_seasonal_for_each_biome(biome4_pi_cube, biome_prism3_cube,
                                  biome_names)

    print('prog finished')
    
 

FILESTART = '/nfs/hera1/earjcti/regridded/'
MODELNAME = 'IPSLCM6A_origgrid' # BIOME4 -  MMM
#MODELNAME = 'xozzb'
ABS_ANOM = 'Anom'
#FILESTART = '/nfs/see-fs-02_users/earjcti/BIOME4/biome4_pliomip2/'
#MODELNAME = ''
#ABS_ANOM = ''
if MODELNAME == 'BIOME4':
    midfname = MODELNAME
elif MODELNAME == 'BIOME4_NEWPARAMETERS':
    midfname = 'BIOME4/BIOME4_NEWPARAMETERS/'
else:
    midfname =  MODELNAME + '/biome4/' 
   
if ABS_ANOM == 'Anom':
    FILENAME = FILESTART + midfname + '/biome4out_anomaly.nc'
if ABS_ANOM == 'Abs':
    FILENAME = FILESTART + midfname + '/biome4out_absolute.nc'
if ABS_ANOM == '':
    FILENAME = FILESTART + midfname + 'biome4out.nc'

main()
