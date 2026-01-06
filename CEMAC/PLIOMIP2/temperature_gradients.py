#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:43:50 2019

@author: earjcti

This will plot the SST zonal and meridional gradients acropss the atlantic and 
the pacific

"""

import sys
import iris
import iris.quickplot as qplt
#from mpl_toolkits.basemap import Basemap, shiftgrid
import numpy as np
import matplotlib.pyplot as plt


def get_data(filereq, field, modeluse):
    """
    gets the field (field) from the file (filereq) and loads it
    into an iris cube (the model name is in modeluse)
    outputs a cube of the data that is as simple as possible
    """

    print(modeluse)
    if modeluse == 'MMM':
        print(filereq,field)
        cube = iris.load_cube(filereq, field)
    else:
        cubes = iris.load(filereq)
        cube = cubes[0]
    cube.data = cube.data.astype('float32')

    if (modeluse == 'MIROC4m' or modeluse == 'COSMOS'):
        cube.units = 'Celsius'
    else:
        cube.convert_units('Celsius')

    for coord in cube.coords():
        name = coord.standard_name
        if name != 'latitude' and name != 'longitude':
            if name is None:
                if coord.long_name is None:
                    cube.remove_coord(coord.var_name)
                else:
                    cube.remove_coord(coord.long_name)
            else:
                cube.remove_coord(name)

    cube.cell_methods = None

    return cube

def get_within_region(cube, region):
    """
    This routine will mask out all the regions that are not 
    input: a temperature cube, a region name
    output: a new cube with all data not within the region masked out.
    """


    if region == 'TROPATL':
        latmin = -20.
        latmax = 20.
        lonmin = 300.
        lonmax = 360.
        
    if region == 'TROPPAC':
        latmin = -20.
        latmax = 20.
        lonmin = 140.
        lonmax = 260.
        
    if region == 'MERIDATL':
        latmin = -70.
        latmax = 70.
        lonmin = 290.
        lonmax = 360.
        
    if region == 'MERIDPAC':
        latmin = -70.
        latmax = 60.
        lonmin = 150.
        lonmax = 260.
        
        
    cubedata = cube.data
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points
    

    for j in range(0, len(lats)):
        cubedata[j, :] = np.ma.masked_where(lons < lonmin, 
                                            cubedata[j, :])
        cubedata[j, :] = np.ma.masked_where(lons > lonmax, 
                                            cubedata[j, :])
  
        
    for i in range(0, len(lons)):
        cubedata[:, i] = np.ma.masked_where(lats < latmin, 
                                            cubedata[:, i])
        cubedata[:, i] = np.ma.masked_where(lats > latmax, 
                                            cubedata[:, i])
        
    newcube = cube.copy(data=cubedata)
    
    #if region == 'MERIDPAC':
    #plot
    #    contour = qplt.contourf(newcube)
    #    plt.gca().coastlines()
    #    plt.show()
    #    sys.exit(0)

    return [newcube, lonmin, lonmax, latmin, latmax]



class Plotalldata:
    """
    This will plot the data from the timeperiod (ie mpwp or pi)
    it can either plot the zonal mean or the meridional mean
    """
    def __init__(self, region, timeperiod, lonmin, lonmax, latmin, latmax):

        """
        inputs are:
                    
        """

        fullregname  =  {
                         "TROPATL" : "Atlantic ",
                         "TROPPAC" : "Pacific ",
                         "MERIDATL" : "Atlantic ",
                         "MERIDPAC": "Pacific "}
    
        colorreq = {
                    "mPWP" : "red",
                    "pi" : "blue",
                    "mPWP-pi" : "green"}
        
        hatchreq = {
                    "mPWP" : "",
                    "pi" : "",
                    "mPWP-pi" : ""}
        
        self.region = region
        self.timeperiod = timeperiod
        self.regiontitle = fullregname.get(region)
        self.colorreq = colorreq.get(timeperiod)
        self.hatchreq = hatchreq.get(timeperiod)


        if region == 'TROPATL':
            #self.regiontitle = (self.regiontitle + 
            #                    np.str(np.int(latmin)) + 
            #                    '-' + np.str(np.int(latmax)))
            if timeperiod == 'mPWP-pi':
                self.valmin = 0.
                self.valmax = 10.
            else:
                self.valmin = 20.
                self.valmax = 30.
                
        if region == 'TROPPAC':
            #self.regiontitle = (self.regiontitle + 
            #                    np.str(np.int(latmin)) + '-' 
            #                    + np.str(np.int(latmax)))
            if timeperiod == 'mPWP-pi':
                self.valmin = 0.
                self.valmax = 10.
            else:
                self.valmin = 20.
                self.valmax = 32.
                
        if region == 'MERIDATL':
            #self.regiontitle = (self.regiontitle + 
            #                    np.str(np.int(lonmin)) + 
            #                    '-' + np.str(np.int(lonmax)))
            if timeperiod == 'mPWP-pi':
                self.valmin = -0.
                self.valmax = 20.
            else:
                self.valmin = -10.
                self.valmax = 30.
                
        if region == 'MERIDPAC':
            #self.regiontitle = (self.regiontitle + 
            #                    np.str(np.int(lonmin)) + '-' 
            #                    + np.str(np.int(lonmax)))
            if timeperiod == 'mPWP-pi':
                self.valmin = 0.
                self.valmax = 20.
            else:
                self.valmin = -5.
                self.valmax = 35.
            



    def plotzm(self, cube, cubemin, cubemax, ax, colorname):
        """
        plot the zonal mean
        """
    

        cube_zm = cube.collapsed('latitude', iris.analysis.MEAN)
        cube_zm_min = cubemin.collapsed('latitude', iris.analysis.MEAN)
        cube_zm_max = cubemax.collapsed('latitude', iris.analysis.MEAN)
        lons = cube_zm.coord('longitude').points
        
        #cube_zm_var = cubevar.collapsed('latitude', iris.analysis.MEAN)
        #zm_std = np.sqrt(cube_zm_var.data)
        #cube_zm_2sigma = cube_zm_var.copy(data=zm_std * 2)
        ax.plot(lons, cube_zm.data, label=self.timeperiod, color=self.colorreq)
        ax.fill_between(lons, cube_zm_min.data, cube_zm_max.data, alpha=0.2, 
                        color=self.colorreq,
                        hatch=self.hatchreq)
        plt.title('20N-20S Mean SST:' +  self.regiontitle)
        ax.set_xlabel('longitude')
        ax.set_ylabel('deg C', color=colorname)
        ax.set_ylim(self.valmin, self.valmax)
        ax.tick_params(axis='y', labelcolor=colorname)
        
    
        return

    def plotmm(self, cube, cubemin, cubemax, ax, colorname, fig):
        """
        plot the meridional mean from the cube
        """

        cube_mm = cube.collapsed('longitude', iris.analysis.MEAN)
        cube_mm_min = cubemin.collapsed('longitude', iris.analysis.MEAN)
        cube_mm_max = cubemax.collapsed('longitude', iris.analysis.MEAN)
        lats = cube_mm.coord('latitude').points
        
        for i, lat in enumerate(lats):
            print(lat, cube_mm.data[i], cube_mm_min.data[i], cube_mm_max.data[i])
        ax.plot(cube_mm.data, lats, label=self.timeperiod, color=self.colorreq)
        ax.fill_betweenx(lats, cube_mm_min.data, cube_mm_max.data, alpha=0.2, 
                        color=self.colorreq,
                        hatch=self.hatchreq)
        fig.suptitle('Zonal Mean SST: ' +  self.regiontitle)
        ax.set_xlabel('deg C', color=colorname)
        ax.set_ylabel('latitude')
        ax.set_xlim(self.valmin, self.valmax)
        ax.tick_params(axis='x', labelcolor=colorname)
   
    
        return


# emd of class

def main_time(timeperiod, region):
    """
    the main routine for a single timeperiod (likely mPWP or PI)
    """

    exptname = {"pi" : "E280",
                "mPWP" : "EOI400"}

    filename = (FILESTART + FIELDNAME + '_multimodelmean.nc')

    mean_cube = get_data(filename, FIELDNAME + 'mean_' + timeperiod, 
                         'MMM')
    
    max_cube = get_data(filename, FIELDNAME + 'max_' + timeperiod, 
                         'MMM')
    
    min_cube = get_data(filename, FIELDNAME + 'min_' + timeperiod, 
                         'MMM')
    
    if timeperiod == 'anomaly':
        sd_cube = get_data(filename, 'SSTanomaly_multimodel_stddev', 
                         'MMM')
    else:
        sd_cube = get_data(filename, FIELDNAME + 'std_' + timeperiod, 
                         'MMM')
    variance_cube = np.square(sd_cube)

    # get the mean value within the region and the standard deviation within the 
    # region
    regioncube, lonmin, lonmax, latmin, latmax = get_within_region(mean_cube, region)
    regioncubemax, lonmin, lonmax, latmin, latmax = get_within_region(max_cube, region)
    regioncubemin, lonmin, lonmax, latmin, latmax = get_within_region(min_cube, region)
    regionvariance, lonmin, lonmax, latmin, latmax = get_within_region(variance_cube, region)


    return regioncube, regionvariance, regioncubemin, regioncubemax, lonmin, lonmax, latmin, latmax

def main():
    """
    the main routine that will split the timeperiod up if appropriate
    if there are two timeperiods it will plot them on the same figure
    """
    
    cubelist = iris.cube.CubeList([])
    
    figno =  {"MERIDATL": "a)",
            "MERIDPAC": "b)",
            "TROPATL": "c)",
            "TROPPAC": "d)"
    }
    
    for j, region in enumerate(REGIONNAMES):
        
        
        
        # get cubes for each timeperiod and plot
        fig, ax1 = plt.subplots() 
        for i in range(0, len(TIMEPERIODS)):
            (cube_time_region, cube_region_variance, 
             cube_region_min, cube_region_max,
             lonmin, lonmax, latmin, latmax)  = main_time(TIMEPERIODS[i], region)
            
            cubelist.append(cube_time_region)
            cubelist.append(cube_region_min)
            cubelist.append(cube_region_max)
       
            plobj = Plotalldata(region, TIMEPERIODS[i], lonmin, lonmax, latmin, latmax)
            if region[0:3] == 'TRO': 
                plobj.plotzm(cube_time_region, cube_region_min,
                             cube_region_max,
                             ax1, 'black')
                outname = OUTSTART + 'zonalmean' + region
                
    
            
                
            if region[0:3] == 'MER':             
                plobj.plotmm(cube_time_region, cube_region_min,
                             cube_region_max,
                             ax1,'black', fig)
                outname = OUTSTART + 'mm_' + region
        
       
        box = ax1.get_position()
        #print(box)
        #sys.exit(0)
        #ax1.set_position([box.x0+(box.height*0.2), box.y0, 1.0, 1.0])
        # plot zonal or meridional mean for the difference
        
        (cubediff_region, cubediff_region_variance, 
         cubediff_region_min, cubediff_region_max,
             lonmin, lonmax, latmin, latmax)  = main_time('anomaly', region)
        print(np.mean(cubediff_region_max.data))
       
        
        plobjdiff = Plotalldata(region, TIMEPERIODS[1] + '-' + TIMEPERIODS[0],
                                lonmin, lonmax, latmin, latmax)
       
        if region[0:3] == 'TRO':             
           ax2 = ax1.twinx() 
           
           plobjdiff.plotzm(cubediff_region, cubediff_region_min,
                            cubediff_region_max, ax2, 'green') 
           
           print('julia',np.mean(cubediff_region.data), 
                 np.mean(cubediff_region_min.data),
                 np.mean(cubediff_region_max.data)
                 )
           #sys.exit(0)
           
           
        if region[0:3] == 'MER':             
            ax2 = ax1.twiny() 
            ax2.xaxis.set_ticks_position("bottom")
            ax2.xaxis.set_label_position("bottom")
            ax2.spines["bottom"].set_position(("axes", -0.2))
            plobjdiff.plotmm(cubediff_region, cubediff_region_min,
                             cubediff_region_max, ax2, 'green', fig) 
           
         # plot the legend and close the plot 
        
        #fig.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        #fig.legend()
        plt.figtext(0.02, 0.97,figno.get(region),
                    horizontalalignment='left',
                    verticalalignment='top',
                    fontsize=20)
        #plt.subplots_adjust(top=0.85)
        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        #plt.show()
        #sys.exit(0)
        plt.savefig(outname + '.eps')
        plt.savefig(outname + '.pdf')
       
        iris.save(cubelist, OUTWRITE, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)        
    
       
    return    
        
        

# variable definition
LINUX_WIN = 'l'
FIELDNAME = 'SST'

if LINUX_WIN == 'l':
    FILESTART = ('/nfs/hera1/earjcti/regridded/')
    OUTSTART = FILESTART + 'allplots/' + FIELDNAME + '/'
    OUTWRITE = FILESTART + 'alldata/data_for_fig4.nc'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\'
    OUTSTART = FILESTART + 'allplots\\' + FIELDNAME + '\\'
    OUTWRITE = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\regridded\\alldata\\data_for_fig4.nc'
    

REGIONNAMES = ['TROPATL', 'TROPPAC', 'MERIDATL', 'MERIDPAC']
#REGIONNAMES = ['TROPATL']

UNITS = 'deg C'
#TIMEPERIODS = ['pi']
TIMEPERIODS = ['pi', 'mPWP']


main()
