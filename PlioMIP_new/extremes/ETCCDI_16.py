
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia

We are looking at ETCCDI Climate change indicies.  This program will deal with indices 16 this is monthly average daily temperature range

DTR, Daily temperature range: Monthly mean difference between TX and TN

Let TXij and TNij be the daily maximum and minimum temperature respectively on day i in month j. If I represents the number of days in j, then:

DTR j = SIGMA(TXij - TN) / 30.  

"""
import numpy as np
import iris
from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
from iris.coords import DimCoord
import matplotlib.pyplot as plt
import sys
from iris.experimental.equalise_cubes import equalise_attributes
from os.path import exists


def get_HCM3_year_data(filestart,year):
    """
    reads in the maximum  and minimum temperature for the year and puts it in 
    a single cube
    """
    DTRcubelist = iris.cube.CubeList([])
    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    for monthno, month in enumerate(months):
        filename = filestart + np.str(year).zfill(2) + month + '.nc'
        # load in data
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp'))
        cube = iris.load(filename, constraints=variable_constraint)
        maxTcube = cube[0]

        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_1'))
        cube = iris.load(filename, constraints=variable_constraint)
        minTcube = cube[0]

        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == 'temp_2'))
        cube = iris.load(filename, constraints=variable_constraint)
        meanTcube = cube[0]

        # check you have got maxT, meanT and minT in correct order        
        if np.max(maxTcube.data) < np.max(meanTcube.data):
            print('cubes not in right order')
            sys.exit(0)

        if np.max(meanTcube.data) < np.max(minTcube.data):
            print('cubes not in right order2')
            sys.exit(0)

        dtrdailycube = maxTcube - minTcube
        dtrmonthcube = dtrdailycube.collapsed('t', iris.analysis.MEAN)
        dtrmonthcube.coord('t').points = monthno + 1
        DTRcubelist.append(dtrmonthcube)
       
    equalise_attributes(DTRcubelist)
    iris.util.unify_time_units(DTRcubelist)
    DTRcube = DTRcubelist.merge_cube()
    DTRcube.coord('ht').rename('year')
    DTRcube.coord('year').points = year
    DTRcube.coord('t').rename('month')
    DTRcube.coord('month').attributes=None
    DTRcube.coord('month').units=None
    DTRcube.coord('month').bounds=None
    DTRcube.cell_methods = None
   
    return DTRcube
   
def calculate_daily_temperature_range():
    """
    for each year and for each month calculates the average daily
    temperature range Tmax - Tmin
    """
    filestart = ('/nfs/hera1/earjcti/um/' + EXPTNAME + '/pb/' + 
                 EXPTNAME + 'a@pb' + EXTRA)
  
    DTR_year_cubes = iris.cube.CubeList([])
    for year in range(0, 100):
        yearuse = np.str(year).zfill(2)
        DTRcube  = get_HCM3_year_data(filestart, year)
        DTR_year_cubes.append(DTRcube)

    equalise_attributes(DTR_year_cubes)
    iris.util.unify_time_units(DTR_year_cubes)
   

    allDTR_cube = DTR_year_cubes.concatenate_cube()
    allDTR_cube.long_name = 'Diurnal Temperature Range'
    allDTR_cube.units='Celsius'

    fileout = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
               EXPTNAME + '/diag16/DTR.nc')
    iris.save(allDTR_cube,fileout)

#########################################################
def plot_daily_temperature_range(ocn_mask):
    """
    plots the daily temperature range in EXPTNAME and also the 
    difference between the experiment and the control
    """ 
    TIME = {'tenvj' : 'mPWP',
            'xozza' : 'PI'}
    filestart = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/'

    # read in control data (normally pi)
    cube_cntl = iris.load_cube(filestart + CNTLNAME + '/diag16/DTR.nc')
    cube_cntl_avg = cube_cntl.collapsed('year',iris.analysis.MEAN)
    # read in experiment data (normally pliocene)
    cube_expt = iris.load_cube(filestart + EXPTNAME + '/diag16/DTR.nc')
    cube_expt_avg = cube_expt.collapsed('year',iris.analysis.MEAN)
    # read in masks in case they are needed
    maskplio=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
    maskpi=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/e280/qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
       
      
   
    for mon in range(1,13):
        cube_cntl_mon = cube_cntl_avg.extract(iris.Constraint(month=mon))
        cube_expt_mon = cube_expt_avg.extract(iris.Constraint(month=mon))
       
        if ocn_mask == 'y':
            if TIME.get(EXPTNAME) == 'mPWP':
                cube_expt_mon.data.mask = (maskplio.data - 1.0) * (-1.0)
            if TIME.get(CNTLNAME) == 'PI':
                cube_cntl_mon.data.mask = (maskpi.data - 1.0) * (-1.0)
      
        plt.subplot(2,2,1)
        qplt.contourf(cube_expt_mon, levels=np.arange(0,22,2), extend='both')
        plt.gca().coastlines()
        plt.title(MONTHNAMES.get(mon-1) + 
                  ': Tmax - Tmin:' +  TIME.get(EXPTNAME))

        plt.subplot(2,2,2)
        qplt.contourf(cube_cntl_mon, levels=np.arange(0,22,2), extend='both')
        plt.gca().coastlines()
        plt.title('Tmax - Tmin:' +  TIME.get(CNTLNAME))

        plt.subplot(2,2,3)
        print(cube_expt_mon)
        print(cube_cntl_mon)
        diffcube = cube_expt_mon - cube_cntl_mon
        qplt.contourf(diffcube,levels=np.arange(-3,3.5,0.5), extend='both',
                      cmap='RdBu_r')
        plt.gca().coastlines()
        plt.title('Tmax - Tmin:' + TIME.get(EXPTNAME) + '-'+ TIME.get(CNTLNAME))
        

        plt.tight_layout()
        fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/' + 
                    'diag16/' + EXPTNAME + '_' + CNTLNAME + '_' + 
                   MONTHNAMES.get(mon-1))
        plt.savefig(fileout + '.eps')
        plt.savefig(fileout + '.png')
        plt.close()
    
   
   
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'
NYEARS = 100
MONTHNAMES = {0:'January',1:'February',2:'March',3:'April',4:'May',5:'June',
              6:'July', 7:'August',8:'September',9:'October',10:'November',
              11:'December'}
  
######################################
#  step 1.  Calculate daily temperature range
EXPTNAME = 'xozzb'
EXTRA='o'

#calculate_daily_temperature_range() # and write to a file

#####################################
# step 2.  Plot daily temperature range
# can only do this after we have done step 1

EXPTNAME = 'xozzb'
CNTLNAME = 'xozza'
EXTRA = 'o'
plot_daily_temperature_range('y') # pass 'y' for mask ocean, 'n' for no ocn mask
