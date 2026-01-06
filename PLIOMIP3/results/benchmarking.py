"""
#NAME
#    benchmarking.py
#PURPOSE 
#
#  This program will create a t-test to see whether the means from
# our climate runs are significantly different between two experiments
#
#NOTES
#  what is a t-test
#  1.  Our sample has 100 years so our data has 99degrees of freedom
#  2.  If our |t-value| >2  then our Null hypothesis (that the distributions are the same)
#                           is rejected at the 95% confidence level/
#  3.  So we are looking for |t-value| < 2
#  4.  t-value = (mean_expt - mean_climate) / s
#      I can't figure out whether s=standard dev or s=standard-dev / sqrt(n)

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
import sys
#from netCDF4 import Dataset, MFDataset

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

def round_to_n(x, n):
    """
    rounds x to n significant figures
    """
    if not x: return 0
    power = -int(math.floor(math.log10(abs(x)))) + (n - 1)
    factor = (10 ** power)
    return round(x * factor) / factor

def get_plot_limits(cube, ndiv, nsig):
    """
    gets the limit for the plot
    imput: a cube and the number of divisions on the plot
    nsig is the number of significant figures to round the difference to
    """

   
    valmin = round_to_n(np.min(cube.data), 2)
    valmax = round_to_n(np.max(cube.data), 2)
        
    valmax = np.min([valmin * -1.0, valmax])
    valmin = valmax * -1.0
    valdiff = round_to_n((valmax - valmin) / ndiv, nsig)
    
    # adjust valmax and valmin so that the difference is an exact multiple
    # of valdiff
    valmax = valdiff * (ndiv / 2.0)
    valmin = -1.0 * valmax
    
       
    return valmin, valmax, valdiff

def fix_mask(cube):
    """
    this is for if the valid min and valid max was wrong and the data has 
    ended up being masked everywhere 
    """
    #1. check if data is masked everywhere
   
    nmask = np.ma.count_masked(cube.data)
    nsize = cube.data.size
    
    if nmask > 0.9 * nsize:  # too much masking unmask everything
        cube.data.mask = False
        
    
        cube.attributes['valid_min'] = -100000.
        cube.attributes['valid_max'] = -100000.
        cube.attributes['fill_value'] = 2.0E20
        
    newcube = cube.copy()
    
    
    return newcube

class main():
    """
    this class will basically run the program
    """
    def __init__(self, field):
        """
        input: field # a short field name
        
        gets data such as the list of the filenames we need
        and the longfieldname and also which type of file it is
        """
 
        if LINUX_WIN == 'w':
            startfname = "C:\\Users\\julia\\OneDrive\\WORK\\DATA\\"
            startout = "C:\\Users\\julia\\OneDrive\\WORK\\DATA\\"
        else:
            startout = "/nfs/hera1/earjcti/um/" + EXPT + "/avgplots/"
            if EXPT_TYPE == 'S':
                startfnamee = "/nfs/hera3/palaeo_share/PlioMIP2/processed/"
            else:
                startfnamee = "/nfs/hera1/earjcti/um/"
            if CNTL_TYPE == 'S':
                startfnamec = "/nfs/hera3/palaeo_share/PlioMIP2/processed/"
            else:
                startfnamec = "/nfs/hera1/earjcti/um/"
            
        longfield = {'temp1.5' : 'TEMPERATURE AT 1.5M',
                     'precip' : 'TOTAL PRECIPITATION RATE     KG/M2/S',
                     'precipmm' : 'TOTAL PRECIPITATION RATE MM/DAY',
                     'cloud_cover' : 'TOTAL CLOUD AMOUNT - RANDOM OVERLAP',
                     'mslp' : 'PRESSURE AT MEAN SEA LEVEL',
                     'mslp_hPa' : 'PRESSURE AT MEAN SEA LEVEL hPa',
                     'evapsea' : 'EVAPORATION FROM SEA (GBM)   KG/M2/S',
                     'seaiceconc' : 'AICE : ICE CONCENTRATION',
                     'icefrac' : 'SEA ICE FRACTION AFTER TIMESTEP',
                     'oceansurftemp' : 'OCN TOP-LEVEL TEMPERATURE          K',
                     'oceansurftempK' : 'OCN TOP-LEVEL TEMPERATURE K',
                     'surfsalinity': 'SALINITY (OCEAN)       (PSU-35)/1000',
                     'surfsalinitypsu': 'SALINITY (OCEAN) (PSU)',
                     'MLD' : 'MIXED LAYER DEPTH (OCEAN)          M',
                     'MLDm' : 'MIXED LAYER DEPTH (OCEAN) M',
                     'AMOC' : 'Meridional Overturning Stream Function (Atlantic)'
                    
                     }
        atm_ocn_ind = {'temp1.5' : 'a',
                       'precip' : 'a',
                       'cloud_cover' : 'a',
                       'mslp' : 'a',
                       'evapsea' : 'a',
                       'seaiceconc' : 'o',
                       'oceansurftemp' :  'o',
                       'surfsalinity': 'o',
                       'MLD' : 'o',
                       'AMOC' : 'o'
                       }
      
        fileletter = {'temp1.5' : 'pd',
                      'precip' : 'pd',  
                      'cloud_cover' : 'pd',
                      'mslp' : 'pd',
                      'evapsea' : 'pd',
                      'seaiceconc' : 'pf',
                      'oceansurftemp' : 'pf',
                      'surfsalinity': 'pf',
                      'MLD' : 'pf'}
        
        if EXPT_TYPE == 'S':
            self.meanfile = startfnamee + EXPT + '_Annual.nc'
            self.stdevfile = startfnamee + EXPT + '_Annual_std.nc'
        else:
            self.meanfile = (startfnamee + EXPT + '/average/' + EXPT + 
                         '_' + field + '_' + SEASON + '_mean_' + 
                         np.str(EXPT_STARTYEAR) + '_' +  
                         np.str(EXPT_STARTYEAR + NYEARS) + '.nc')
            self.stdevfile = (startfnamee + EXPT + '/stdev/' + EXPT + 
                          '_' + field + '_' + SEASON + '_stdev_' +
                          np.str(EXPT_STARTYEAR) + '_' +  
                          np.str(EXPT_STARTYEAR + NYEARS) + '.nc')
        print(self.meanfile)
       
        
        # get benchmarking data
        if CNTL_TYPE == 'B':
            self.meanbench = (startfnamec +  'Benchmarking/' + 
                              CNTL + atm_ocn_ind.get(field) + 
                              '.' +  fileletter.get(field) + 'cl' + 
                              SEASON + '.nc')
            self.sdbench = (startfnamec + 'Benchmarking/' + 
                              CNTL + atm_ocn_ind.get(field) + 
                              '.' +  fileletter.get(field) + 'sd' + 
                              SEASON + '.nc')
            
        elif CNTL_TYPE == 'S':
            self.meanbench = startfnamec + CNTL + '_Annual.nc'
            self.sdbench = startfnamec + CNTL + '_Annual_std.nc'
            
        else:
            self.meanbench = (startfnamec + CNTL + '/average/' + CNTL + 
                              '_' + field + '_' + SEASON +'_mean_' +
                               np.str(CNTL_STARTYEAR) + '_' +  
                               np.str(CNTL_STARTYEAR + NYEARS) + '.nc')
            self.sdbench = (startfnamec + CNTL + '/stdev/' + CNTL + 
                            '_' + field + '_' + SEASON +'_stdev_' +
                            np.str(CNTL_STARTYEAR) + '_' +  
                            np.str(CNTL_STARTYEAR + NYEARS) + '.nc')
                           
            
        self.fieldname = longfield.get(field)
        if CNTL_TYPE=='S':
            cntl_out=CNTL
        else:
            cntl_out = (CNTL + '_' + np.str(CNTL_STARTYEAR) + 
                     '_' +  np.str(CNTL_STARTYEAR + NYEARS) + '_')
        if EXPT_TYPE=='S':
            expt_out=EXPT
        else:
            expt_out = (EXPT + '_' + np.str(EXPT_STARTYEAR) + 
                     '_' +  np.str(EXPT_STARTYEAR + NYEARS) + '_')
                       
        self.outstart = (startout + expt_out +'_' + cntl_out + 
                         '_' + field  + '_' + SEASON)
        self.field = field
        
      
           
    def main_benchmark_field(self):
        """
        this will compare the field with the benchmarked field and
        plot to a file
        """
        
        cube_exptmean = iris.load_cube(self.meanfile, self.fieldname)
        cube_exptsd = iris.load_cube(self.stdevfile, self.fieldname)
      
        if CNTL_TYPE == 'B':
            variable_constraint = iris.Constraint(cube_func=(lambda c: c.long_name == self.fieldname))
            cube_cntlmean = iris.load_cube(self.meanbench, variable_constraint)
            cube_temp = iris.load_cube(self.sdbench, variable_constraint)
            cube_cntlsd = fix_mask(cube_temp)
        else:
            cube_cntlmean = iris.load_cube(self.meanbench, self.fieldname)
            cube_cntlsd = iris.load_cube(self.sdbench, self.fieldname)
        
        # if steve#s salinity we just need two levels
        if self.field ==  'surfsalinitypsu':
            #cube_exptmean.coord('z2').var_name = 'depth'
            #cube_exptmean.coord('depth').standard_name = 'depth'
            #print(cube_exptmean.coords())
            cube_exptmean = cube_exptmean.extract(iris.Constraint(depth=0.0))
            print(cube_exptmean)
            sys.exit(0)
            cube_exptsd = cube_exptsd.extract(iris.Constraint(z2=0.0))
            cube_cntlmean = cube_cntlmean.extract(iris.Constraint(z2=0.0))
            cube_cntlsd = cube_cntlsd.extract(iris.Constraint(z2=0.0))
            
            
       
        # makes the exptmean, exptsd, cntlmean and cntlsd cubes all look the 
        # same
        cube_exptmean, cbarunits = self.mean_equalise_cubes(cube_exptmean)
        cube_cntlmean, cbarunits = self.mean_equalise_cubes(cube_cntlmean)
        cube_exptsd, cbarunits = self.mean_equalise_cubes(cube_exptsd)
        cube_cntlsd, cbarunits = self.mean_equalise_cubes(cube_cntlsd)
        
        cube_anom = cube_exptmean - cube_cntlmean
        # to calculate the divisor in the t-statistic
        # you will have to look at how to do a t-test on the web for exlanation
        data1 = ((((cube_exptsd.data ** 2.0) * (NYEARS-1.0))
                          + ((cube_cntlsd.data ** 2.0) * (NYEARS-1.0)))
                     / (2.0*(NYEARS-1.0)))
        data2 = (np.sqrt(data1)) * (np.sqrt(2.0/NYEARS))
        cube_div = cube_anom.copy(data=data2)
        
       
        cube_tstat = cube_anom / cube_div
        print(cube_tstat.data)
       
      
        if self.field == 'AMOC':
            self.main_plot_AMOC(cube_anom, cbarunits, cube_tstat)
        else:
            self.main_plot_latlon(cube_anom, cbarunits, cube_tstat)
        # plot the anomaly
           
    def main_plot_AMOC(self, cube_anom, cbarunits, cube_tstat):
        """
        do a lat-depth plot of the field
        we are plotting the anomaly and the tvalue
        """
        
        ax = plt.subplot(1,2,1)
      
        anommin, anommax, anomdiff = get_plot_limits(cube_anom, 10, 1)
        cs = iplt.contourf(cube_anom, np.arange(anommin, anommax+anomdiff, anomdiff), 
                      cmap='RdYlBu_r', extend='both')
        cbar=plt.colorbar(cs,orientation="horizontal")
        plt.title(self.field + ' anomaly', fontsize=10)
        cbar.set_label('SV')
        
        # plot the tstatistic
        ax2 = plt.subplot(1, 2, 2)
        tmin, tmax, tdiff = get_plot_limits(cube_tstat, 10, 1)
        cs = iplt.contourf(cube_tstat, np.arange(tmin, tmax+tdiff, tdiff), cmap='RdYlBu_r', extend='both')
        cbar=plt.colorbar(cs,orientation="horizontal")
        
        iplt.contour(cube_tstat, [-np.inf, -2, 2, np.inf], hatches=[3*'\\\'', None, 3*'\\\''], colors='black')
        iplt.contourf(cube_tstat, [-np.inf, -2, 2, np.inf], hatches=[3*'\\\'', None, 3*'\\\''], colors='none')
        plt.title(self.field + ' t-value', fontsize=10)
        
        
        
        plt.tight_layout()
        print(self.outstart)
        plt.savefig(self.outstart + '.png')
        plt.savefig(self.outstart + '.pdf')
        #plt.show()
        plt.close()
        
       
      
        
    def main_plot_latlon(self, cube_anom, cbarunits, cube_tstat):
        """
        do a latlon plot of the field
        we are plotting the anomaly and the tvalue
        """

        ax = plt.subplot(1, 2, 1, projection = ccrs.PlateCarree())
        anommin, anommax, anomdiff = get_plot_limits(cube_anom, 10, 1)
        print(anommin, anommax, anomdiff)
        print(np.arange(anommin, anommax+anomdiff, anomdiff))
        
        if cube_anom.data.ndim==3:
            cube_anom = cube_anom[0,:,:]
            cube_tstat=cube_tstat[0,:,:]
       
        #levels=np.arange(anommin, anommax+anomdiff, anomdiff) 
        levels=np.arange(-2.0, 2.4, 0.4)
        cs = iplt.contourf(cube_anom, levels, 
                      cmap='RdYlBu_r', extend='both')
        cbar=plt.colorbar(cs,orientation="horizontal")
        if cube_anom.units == '':
            cbar.set_label(cbarunits)
        else:
            cbar.set_label(cube_anom.units)
        cbar.ax.tick_params(labelsize=8, labelrotation=60) 
        plt.title(self.fieldname + ' anomaly', fontsize=10)
        plt.gca().coastlines()
        #gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-', draw_labels='true')
        #gl.xlabels_top = False
        #gl.ylabels_left = False
        #gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
        #gl.xformatter = LONGITUDE_FORMATTER
        #gl.yformatter = LATITUDE_FORMATTER
        
        
        # plot the tstatistic
        ax2 = plt.subplot(1, 2, 2, projection = ccrs.PlateCarree())
        tmin, tmax, tdiff = get_plot_limits(cube_tstat, 10, 1)
        cs = iplt.contourf(cube_tstat, np.arange(tmin, tmax+tdiff, tdiff), cmap='RdYlBu_r', extend='both')
        cbar=plt.colorbar(cs,orientation="horizontal")
        cbar.ax.tick_params(labelsize=8, labelrotation=60) 
        iplt.contour(cube_tstat, [-np.inf, -2, 2, np.inf], hatches=[3*'\\\'', None, 3*'\\\''], colors='black')
        iplt.contourf(cube_tstat, [-np.inf, -2, 2, np.inf], hatches=[3*'\\\'', None, 3*'\\\''], colors='none')
        plt.title(self.fieldname + ' t-value', fontsize=10)
        plt.gca().coastlines()
        #gl = ax.gridlines(crs=ccrs.PlateCarree(), linestyle='-', draw_labels='true')
        #gl.xlabels_top = False
        #gl.ylabels_left = False
        #gl.ylocator = mticker.FixedLocator([-90, -45, 0, 45, 90])
        #gl.xformatter = LONGITUDE_FORMATTER
        #gl.yformatter = LATITUDE_FORMATTER
        
        plt.tight_layout()
    
        plt.savefig(self.outstart + '.png')
        plt.savefig(self.outstart + '.eps')
        #plt.show()
        plt.close()
        
       
          
    def mean_equalise_cubes(self, cube):
        """
        makes the 2 cubes in the main program look the same for analysis
        they will look the same as cube_expt
        """    
        
        # check dimensions
        newcube = iris.util.squeeze(cube)
        
        
        # if cube is precipitation: change to mm/month
        if self.field == 'precip':
           
           newcube.long_name = 'Precipitation'
           newcube.name('precipitation_flux')
           newcube.title = 'Precipitation'
           self.fieldname = 'Precipitation'
           #newcube.units = 'kg m-2 s-1'
           #newcube.convert_units('mm day-1')
           newcube.data = newcube.data * 60. * 60. * 24.
           newcube.units = 'mm day-1'
           
        # check names etc
        cbarunits = ''
        if self.field == 'cloud_cover':
            self.fieldname = 'Cloud amount'
            newcube.units = ''
            cbarunits = 'fraction'
            
        if self.field == 'seaiceconc':
            newcube.units = ''
            cbarunits = 'fraction'
            
        if self.field == 'oceansurftemp':
            self.fieldname = 'Ocean Surface Temp'
        
        if self.field == 'surfsalinity':
            self.fieldname = 'Salinity'
            newcube.units = ''
            newcube.data = newcube.data * 1000.
            cbarunits = 'psu'
            
        if self.field == 'MLD':
            self.fieldname = 'Mixed Layer Depth'
            #cbarunits = 'm'
            
        if self.field == 'evapsea':
            self.fieldname = 'Evaporation from sea'
            newcube.data = newcube.data * 60. * 60. * 24.
            newcube.units = 'mm day-1'
           
           
        
 
               
        return newcube, cbarunits
       
    
#=================================================================
# MAIN PROGRAM STARTS HERE



LINUX_WIN='l'
NYEARS = 100
SEASON = 'ann'

# data from new experiemnt
MODELTYPE = 'y' # n=HadGEM, y=HadCM3, F=Famous

EXPT = 'xpsir'  # xpsic PI,  xpsij-lp490  xpsik - lp560
EXPT_STARTYEAR = 500
#EXPT = 'Eoi400_ARC4_2450-2499'
EXPT_TYPE = 'A'#  A = average file like I made using 
               #      program COMMON/basic_plots/avg_fields.py
               # B = Bristol file like Paul made
               # S = average file like Steve made

# data from good experiment
CNTL_TYPE = 'A' # A = average file like I made
               # B = Bristol file like Paul made
                 # S = like file steve made
CNTL = 'xpsid'  # xpsic pi, xpsid lp400
CNTL_STARTYEAR = 500
#CNTL = 'Eoi400_2450-2499'
#CNTL = 'tcfze'
#CNTL_STARTYEAR = 0

FIELDS  = ['temp1.5','precip','cloud_cover','mslp',
          'seaiceconc','oceansurftemp',
          'surfsalinity','MLD', 'evapsea']

#FIELDS  = ['surfsalinity','MLD', 'evapsea']

#FIELDS_STEVE = ['temp1.5','precipmm','cloud_cover','mslp_hPa',
#          'icefrac','oceansurftempK',
#          'surfsalinity','MLD', 'evapsea']

#FIELDS = ['AMOC']
#FIELDS = ['temp1.5']

for i, field in enumerate(FIELDS):
    print(field)
    obj = main(field)
    obj.main_benchmark_field() # benchmarks this field

    
