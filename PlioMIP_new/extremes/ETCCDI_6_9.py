#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 07.02.2022 by Julia

We are looking at ETCCDI Climate change indicies.  This program will write
indices 6-9 to a file.  These are:

6  TXx, Monthly maximum value of daily maximum temperature:

Let TXx be the daily maximum temperatures in month k, period j. The maximum daily maximum temperature each month is then:

TXxkj=max(TXxkj)

7.  TNx, Monthly maximum value of daily minimum temperature:

Let TNx be the daily minimum temperatures in month k, period j. The maximum daily minimum temperature each month is then:

TNxkj=max(TNxkj)

8.TXn, Monthly minimum value of daily maximum temperature:

Let TXn be the daily maximum temperatures in month k, period j. The minimum daily maximum temperature each month is then:

9. TXnkj=min(TXnkj)

TNn, Monthly minimum value of daily minimum temperature:

Let TNn be the daily minimum temperatures in month k, period j. The minimum daily minimum temperature each month is then:

TNnkj=min(TNnkj) 

"""
import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import matplotlib.pyplot as plt
import sys


def get_HCM3_year_data(filestart,year):
    """
    reads in the maximum  and minimum temperature for the year and puts it in 
    a single cube
    """
    maxmon_maxT_cubelist = iris.cube.CubeList([])
    minmon_maxT_cubelist = iris.cube.CubeList([])
    maxmon_minT_cubelist = iris.cube.CubeList([])
    minmon_minT_cubelist = iris.cube.CubeList([])

    months = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    for month in months:
        filename = filestart + np.str(year).zfill(2) + month + '.nc'
        # load in data
        cubes=iris.load(filename)
        temperature_cubes = iris.cube.CubeList([])
        max_temp = []
        for cube in cubes:
            if (cube.var_name == 'temp' or cube.var_name == 'temp_1' 
                or cube.var_name == 'temp_2'):
               try:
                   cube.coord('t_1').rename('t')
               except:
                   pass
               temperature_cubes.append(cube)
               max_temp.append(np.max(cube.data))
     
        maxindex = np.argmax(max_temp)
        minindex = np.argmin(max_temp)
        indexes = np.argsort(max_temp)
        minTcube = temperature_cubes[indexes[0]]
        meanTcube = temperature_cubes[indexes[1]]
        maxTcube = temperature_cubes[indexes[2]]
    
       
        # check you have got maxT, meanT and minT in correct order        
        if np.max(maxTcube.data) < np.max(meanTcube.data):
            print('cubes not in right order')
            sys.exit(0)

        if np.max(meanTcube.data) < np.max(minTcube.data):
            print('cubes not in right order2')
            sys.exit(0)

        maxmon_of_maxT_cube = maxTcube.collapsed('t', iris.analysis.MAX)
        minmon_of_maxT_cube = maxTcube.collapsed('t', iris.analysis.MIN)
        maxmon_of_minT_cube = minTcube.collapsed('t', iris.analysis.MAX)
        minmon_of_minT_cube = minTcube.collapsed('t', iris.analysis.MIN)

        maxmon_maxT_cubelist.append(iris.util.squeeze(maxmon_of_maxT_cube))
        minmon_maxT_cubelist.append(iris.util.squeeze(minmon_of_maxT_cube))
        maxmon_minT_cubelist.append(iris.util.squeeze(maxmon_of_minT_cube))
        minmon_minT_cubelist.append(iris.util.squeeze(minmon_of_minT_cube))
        
    equalise_attributes(maxmon_maxT_cubelist)
    equalise_attributes(minmon_maxT_cubelist)
    equalise_attributes(maxmon_minT_cubelist)
    equalise_attributes(minmon_minT_cubelist)

    #for cube in maxmon_maxT_cubelist:
    #    print(cube.coord('t'))
    max_maxdayTcube = maxmon_maxT_cubelist.merge_cube()
    min_maxdayTcube = minmon_maxT_cubelist.merge_cube()
    max_mindayTcube = maxmon_minT_cubelist.merge_cube()
    min_mindayTcube = minmon_minT_cubelist.merge_cube()
  
    return max_maxdayTcube, min_maxdayTcube, max_mindayTcube, min_mindayTcube 
   

def get_HadCM3_diagnostics(expt, extra):
    """
    gets the diagnostics (frost days, summer days, icing days tropical nights
    from HadCM3)
    """
    filestart = '/nfs/hera1/earjcti/um/' + expt + '/pb/' + expt + 'a@pb' + extra
  
    for year in range(0, 100):
        yearuse = np.str(year).zfill(2)
        (max_Tmax_cube, min_Tmax_cube, 
         max_Tmin_cube, min_Tmin_cube)= get_HCM3_year_data(filestart, year)

        # rename
        max_Tmax_cube.long_name = 'Monthly maximum value of daily maximum temperature'
        max_Tmin_cube.long_name = 'Monthly maximum value of daily minimum temperature'
        min_Tmax_cube.long_name = 'Monthly minimum value of daily maximum temperature'
        min_Tmin_cube.long_name = 'Monthly minimum value of daily minimum temperature'

        cubelist = [max_Tmax_cube, max_Tmin_cube,min_Tmax_cube,
                    min_Tmin_cube]

        outfile = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
                   expt + '/' + extra + '_' + 
                   'diag6-9_' + np.str(year) + '.nc')
        iris.save(cubelist, outfile, netcdf_format="NETCDF3_CLASSIC")


#########################################################################
def read_data(expt,extra,startyear,endyear):
    """
    reads in the data for each year, finds the sum and returns
    """
    filestart = '/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/'

    maxTmaxcubelist = iris.cube.CubeList([])
    maxTmincubelist = iris.cube.CubeList([])
    minTmaxcubelist = iris.cube.CubeList([])
    minTmincubelist = iris.cube.CubeList([])
 
    for year in range(startyear,endyear):
        filename = (filestart + expt + '/diag6-9/' + extra + 
                    '_diag6-9_'+ np.str(year) + '.nc')
        cube = iris.load_cube(filename,
                        'Monthly maximum value of daily maximum temperature')
        maxTmaxcubelist.append(cube[MONTH_REQ, :,:])
            
        cube = iris.load_cube(filename,
                        'Monthly maximum value of daily minimum temperature')
        maxTmincubelist.append(cube[MONTH_REQ,:,:])

        cube = iris.load_cube(filename,
                        'Monthly minimum value of daily maximum temperature')
        minTmaxcubelist.append(cube[MONTH_REQ, :, :])

        cube = iris.load_cube(filename,
                        'Monthly minimum value of daily minimum temperature')
        minTmincubelist.append(cube[MONTH_REQ, :, :])

    equalise_attributes(maxTmaxcubelist)
    all_maxTmax_cube = maxTmaxcubelist.merge_cube()
    max_maxTmax_cube = all_maxTmax_cube.collapsed('t',iris.analysis.MAX)
    mean_maxTmax_cube = all_maxTmax_cube.collapsed('t',iris.analysis.MEAN)
  
    equalise_attributes(minTmaxcubelist)
    all_minTmax_cube = minTmaxcubelist.merge_cube()
    mean_minTmax_cube = all_minTmax_cube.collapsed('t',iris.analysis.MEAN)
  
    equalise_attributes(maxTmincubelist)
    all_maxTmin_cube = maxTmincubelist.merge_cube()
    max_maxTmin_cube = all_maxTmin_cube.collapsed('t',iris.analysis.MAX)
    mean_maxTmin_cube = all_maxTmin_cube.collapsed('t',iris.analysis.MEAN)
  
    equalise_attributes(minTmincubelist)
    all_minTmin_cube = minTmincubelist.merge_cube()
    min_minTmin_cube = all_minTmin_cube.collapsed('t',iris.analysis.MIN)
    mean_minTmin_cube = all_minTmin_cube.collapsed('t',iris.analysis.MEAN)
  
  
    return (mean_maxTmax_cube, mean_minTmax_cube, mean_maxTmin_cube,
            mean_minTmin_cube, max_maxTmax_cube, min_minTmin_cube)
  
##########################################################   
def plot_extreme_extremes(meanpliocube, meanpicube,meananomcube, extrpliocube,
                          extrpicube,
                          extranomcube, plottype, ocn_mask, expt,cntl):
    """
    this will do a four panel plot
    a) pliocene_annmean extreme temperature  b) pliocene-pi anomaly of a
    c) plioceene most extreme temp in 100 yrs d) plio - pi anomaly of c
    """
    months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    exptnames = {'xozzm' : 'E560', 'xozza': 'PI'}

    if ocn_mask == 'y':
       maskcube=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
      
       meanpliocube.data.mask = (maskcube.data - 1.0) * (-1.0)
       meanpicube.data.mask = (maskcube.data - 1.0) * (-1.0)
       meananomcube.data.mask = (maskcube.data - 1.0) * (-1.0)
       extrpliocube.data.mask = (maskcube.data - 1.0) * (-1.0)
       extranomcube.data.mask = (maskcube.data - 1.0) * (-1.0)
       extrpicube.data.mask = (maskcube.data - 1.0) * (-1.0)
      
    plt.subplot(2,2,1)
    if plottype == 'max':
        vals = np.arange(30,55,5)
        pivals = [40]
        picolor='white'
        cmapname = 'gist_stern_r'
    if plottype == 'min':
        vals = np.arange(-50,0,10)
        pivals = [-30]
        picolor='black'
        cmapname='gist_ncar'
    meanpliocube.convert_units('celsius')
    meanpicube.convert_units('celsius')
    meanpliocube.long_name = exptnames.get(expt) + ': Mean ' + plottype + ' of T' + plottype
   
    iplt.contourf(meanpliocube, extend='both',levels=vals,cmap=cmapname)
    iplt.contour(meanpicube, levels=pivals, linestyles='solid',
                 colors=picolor,linewidths=1)
    plt.gca().coastlines()
    print(np.str(pivals[0]))
    plt.title(exptnames.get(expt) + ' mean ' + plottype + ' ' + months[MONTH_REQ] + ' T'+
              plottype +  '(contour: pi='+ np.str(pivals[0])+'degC)',fontsize=8)

    plt.subplot(2,2,2)
    vals_a = np.arange(-12,14,2)
    meananomcube.long_name = exptnames.get(expt) + ' - ' + exptnames.get(cntl) + ': Mean ' + plottype + ' of T' + plottype  
    iplt.contourf(meananomcube, extend='both',levels=vals_a,cmap='RdBu_r')
    plt.gca().coastlines()
    plt.title(exptnames.get(expt) + '-' + exptnames.get(cntl) + 'ann mean ' + plottype + ' ' + months[MONTH_REQ] + ' T',
              fontsize=8)

    plt.subplot(2,2,3)
    extrpliocube.convert_units('celsius')
    extrpicube.convert_units('celsius')
    extrpliocube.long_name = (exptnames.get(expt) + ' ' + plottype + ' ' + months[MONTH_REQ] + 
                              ' T' + plottype + ' in ' + np.str(NYEARS) + 
                              'years')
    qplt.contourf(extrpliocube, extend='both',levels=vals,cmap=cmapname)
    iplt.contour(extrpicube, levels=pivals, linestyles='solid',
                 colors=picolor,linewidths=1)
  
    plt.gca().coastlines()
    plt.title(exptnames.get(expt) + ' ' + plottype + ' ' + months[MONTH_REQ] + ' T in '
              + np.str(NYEARS) + 'years',fontsize=8)

    plt.subplot(2,2,4)
    extranomcube.long_name = (plottype + ' ' + months[MONTH_REQ] + 
                              ' T' + plottype + ' in ' + np.str(NYEARS) + 
                              'years : '+ exptnames.get(expt) + '-' + 
                              exptnames.get(cntl))
  
    qplt.contourf(extranomcube, extend='both',levels=vals_a,cmap='RdBu_r')
  
    plt.gca().coastlines()
    plt.title(plottype + ' ' + months[MONTH_REQ] + ' T in '
              +np.str(NYEARS) + 'years: ' + exptnames.get(expt) + '-' + 
              exptnames.get(cntl) + ' anom',fontsize=8)

    
    plt.tight_layout()
   
   
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag6-9/' +
               expt + '_' + cntl + '_' + 
               plottype + 'T' + plottype + '_' + months[MONTH_REQ])
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()
    # save data to a netcdffile
   
    meanpliocube.data.mask = None
    meananomcube.data.mask = None
    extrpliocube.data.mask = None
    extranomcube.data.mask = None
  
    cubelist = [meanpliocube, meananomcube, extrpliocube, extranomcube]
    iris.save(cubelist, fileout + '.nc', netcdf_format="NETCDF3_CLASSIC")

##########################################################              
def  plot_minTmax_maxTmin(mean_minTmax_pliocube, mean_minTmax_anomcube,
                  mean_maxTmin_pliocube,mean_maxTmin_anomcube,ocn_mask,
                  expt,cntl):
    """
    this will do a four panel plot
    a) mean_minTmax_plio, b) mean_minTmax_anom
    c) mean_maxTmin_plio, d) mean_maxTmin_anom
    """
    
    months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
 
    mean_minTmax_pliocube.convert_units('celsius')
    mean_maxTmin_pliocube.convert_units('celsius')
    if ocn_mask == 'y':
       maskcube=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
      
       mean_minTmax_pliocube.data.mask = (maskcube.data - 1.0) * (-1.0)
       mean_minTmax_anomcube.data.mask = (maskcube.data - 1.0) * (-1.0)
       mean_maxTmin_pliocube.data.mask = (maskcube.data - 1.0) * (-1.0)
       mean_maxTmin_anomcube.data.mask = (maskcube.data - 1.0) * (-1.0)
       
    plt.subplot(2,2,1)
    vals = np.arange(-50,45,5)
    name = 'Plio minTmax:' + months[MONTH_REQ]
    mean_minTmax_pliocube.long_name = name
    qplt.contourf(mean_minTmax_pliocube, extend='both',levels=vals)
    plt.gca().coastlines()
    plt.title(name)

    plt.subplot(2,2,2)
    vals_a = np.arange(-12,14,2)
    name = 'Plio - PI: minTmax '+ months[MONTH_REQ]
    mean_minTmax_anomcube.long_name = name
    qplt.contourf(mean_minTmax_anomcube, extend='both',
                  levels=vals_a,cmap='RdBu_r')
    plt.gca().coastlines()
    plt.title(name)
    
    plt.subplot(2,2,3)
    name = 'Plio maxTmin:' + months[MONTH_REQ]
    mean_maxTmin_pliocube.long_name = name
    qplt.contourf(mean_maxTmin_pliocube, extend='both',levels=vals)
    plt.gca().coastlines()
    plt.title(name)

    plt.subplot(2,2,4)
    name = 'Plio - PI: maxTmin '+ months[MONTH_REQ]
    mean_maxTmin_anomcube.long_name = name
    qplt.contourf(mean_maxTmin_anomcube, extend='both',
                  levels=vals_a,cmap='RdBu_r')
    plt.gca().coastlines()
    plt.title(name)

    
    plt.tight_layout()
   
   
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag6-9/'+ 
               expt + '_' + cntl + '_maxTmin_minTmax' + months[MONTH_REQ])
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()
    # save data to a netcdffile
   
    mean_maxTmin_pliocube.data.mask = None
    mean_maxTmin_anomcube.data.mask = None
    mean_minTmax_pliocube.data.mask = None
    mean_minTmax_anomcube.data.mask = None
   
    cubelist = [mean_maxTmin_pliocube, mean_maxTmin_anomcube, 
                mean_minTmax_pliocube, mean_minTmax_anomcube]
    iris.save(cubelist, fileout + '.nc', netcdf_format="NETCDF3_CLASSIC")


def plot_textremes_toplevel(expt,exptextra,cntl,cntlextra):
    """
    this will plot the extreme temperatures (maximum of Tmax, minimum of Tmax,
    maximum of Tmin, minimum of Tmin) by month and show how these 
    have changed between the PI and the Pliocene
    """

    (mean_maxTmax_pliocube, mean_minTmax_pliocube, 
     mean_maxTmin_pliocube, mean_minTmin_pliocube, max_maxTmax_pliocube, 
     min_minTmin_pliocube) = read_data(expt,exptextra,0, 0+NYEARS)

    (mean_maxTmax_picube, mean_minTmax_picube, 
     mean_maxTmin_picube, mean_minTmin_picube, 
     max_maxTmax_picube, min_minTmin_picube)= read_data(cntl,cntlextra,0, 0+ NYEARS)
              

    # plot maximum temperature of the Tmax for this month
    plot_extreme_extremes(mean_maxTmax_pliocube, mean_maxTmax_picube,
                  mean_maxTmax_pliocube - mean_maxTmax_picube,
                  max_maxTmax_pliocube, max_maxTmax_picube,
                  max_maxTmax_pliocube - max_maxTmax_picube,'max','y',expt,cntl)

    # plot minimum temperature of Tmin for this month
    plot_extreme_extremes(mean_minTmin_pliocube, mean_minTmin_picube,
                  mean_minTmin_pliocube - mean_minTmin_picube,
                  min_minTmin_pliocube, min_minTmin_picube,
                  min_minTmin_pliocube - min_minTmin_picube,'min','y',expt,cntl)

    # plot minimum temperature of Tmax and maximum temperature of Tmin
    # for this month
    plot_minTmax_maxTmin(mean_minTmax_pliocube,
                  mean_minTmax_pliocube - mean_minTmax_picube,
                  mean_maxTmin_pliocube,
                  mean_maxTmin_pliocube - mean_maxTmin_picube,
                 'y',expt,cntl)

#########################################################
def plot_all_months(field, expt, cntl, ocn_mask):
    """
    this is for plotting the anomalies in the diagnostics for all months
    on one page
    ie max_Tmax_plio - max_Tmax_pi (mean value) for all months
    field : maxTmax, maxTmin, minTmax, minTmin
    """ 

    fig = plt.figure(figsize=(11.75, 13.0))
    # get all the stuff (names/ files etc) we might need
    filestart = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/' + 
                 'diag6-9/' + expt + '_' + cntl + '_')
    months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
    fieldname = {"maxTmax": "E560 - PI: Mean max of Tmax",
                 "minTmin": "E560 - PI: Mean min of Tmin",
                 "minTmax": "E560 - PI: minTmax ",
                 "maxTmin": "E560 - PI: maxTmin "}
    filemid = {"maxTmax": "maxTmax_", "minTmin": "minTmin_",
               "minTmax": "maxTmin_minTmax",
               "maxTmin": "maxTmin_minTmax"}
    maskcube=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
    fieldstart = fieldname.get(field)
    vals = np.arange(-12,14,2)
   
    # process and plot
    for mon in range(0,12):
        filename = filestart + filemid.get(field) + months[mon] + '.nc'
        if field == "minTmax" or field == "maxTmin":
            fielduse = fieldstart + months[mon]
        else:
            fielduse = fieldstart 
        test = iris.load(filename)
        cube = iris.load_cube(filename,fielduse)
        if ocn_mask == 'y':
           cube.data.mask = (maskcube.data - 1.0) * (-1.0)
      
        plt.subplot(4,3,mon+1)
        cs = iplt.contourf(cube, levels=vals, cmap='RdBu_r', extend='both')
        plt.gca().coastlines()
        plt.title(months[mon] + ': Plio-Pi ' + field) 
    
    # tidy up plot and add colorbar
    plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.95,
                                wspace=0.1, hspace=0.0)

    cb_ax = fig.add_axes([0.35, 0.05, 0.30, 0.02])
           
    cbar = fig.colorbar(cs, cax=cb_ax, orientation='horizontal')
    cbar.set_label('degC', fontsize=10)
         
    cbar.ax.tick_params(labelsize=10)
  

    # save to file
    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag6-9/allmonths_'+  expt + '_' + cntl + '_' + field)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()
  
def extract_avg(latreq,lonreq,expt,extra):
    """
    extract average at a given location
    """
    filestart = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + expt + 
          '/diag6-9/')
    max_Tmax = []
    max_Tmin = []
    min_Tmax = []
    min_Tmin = []
    constraint = iris.Constraint(longitude=lonreq,latitude=latreq)
    for year in range(0,100):
        filename = filestart + extra + '_diag6-9_' + str(year) + '.nc'
        cubes = iris.load(filename)
        for cube in cubes:
            smallcube = cube.extract(constraint)
            if cube.long_name == 'Monthly minimum value of daily maximum temperature':
                min_Tmax.append(np.min(smallcube.data)-273.15)
            if cube.long_name == 'Monthly maximum value of daily maximum temperature':
                max_Tmax.append(np.max(smallcube.data)-273.15)
            if cube.long_name == 'Monthly minimum value of daily minimum temperature':
                min_Tmin.append(np.min(smallcube.data)-273.15)
            if cube.long_name == 'Monthly maximum value of daily minimum temperature':
                max_Tmin.append(np.max(smallcube.data)-273.15)
      
    print('mean maxTmax',np.mean(max_Tmax))
    print('mean maxTmin',np.mean(max_Tmin))
    print('mean minTmax',np.mean(min_Tmax))
    print('mean minTmin',np.mean(min_Tmin))
##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'

###################################################
# this is for if you actually want to get the diagnostics.
# ie write them to the file /nfs/hera1/earjcti/PLIOMIP2/..ETCCDI....diags 1-4  
#if MODELNAME == 'HadCM3':
#    get_HadCM3_diagnostics('xozzm','w')

##################################################################
###  this is for plotting the extremes on a month by month basis
#NYEARS = 100
#for MONTH_REQ in range(0,12):
#    plot_textremes_toplevel('xozzm','w','xozza','o')
 

###############################################################
### this is for plotting the anomalies in the diagnostics for all months
### on one page
### ie max_Tmax_plio - max_Tmax_pi (mean value) for all months
#field='minTmax' # values maxTmax,  maxTmin, maxTmin, minTmax
#plot_all_months(field, 'xozzb', 'xozza','y')  


############################################################
# extract average at a given location

extract_avg(52.5,0,'xozzm','w')
