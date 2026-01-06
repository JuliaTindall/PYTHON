#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:43:50 2019

@author: earjcti

This will plot the SST meridional gradient temperature change
vs global SSTtemperature change for each of the PlioMIP models

changed on October 17th to add the data to the figure

"""

import sys
import iris
import iris.quickplot as qplt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


class GetProxy:
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
            self.filename = FILESTARTP + 'pliovar_uk37_ori_vs_bayspline.xlsx'
            self.bsloc = 8
        if datatype == 'MGCA':
            self.filename = FILESTARTP +  'pliovar_mgca_OrivsBaymag.xlsx'
            self.bsloc = 7
        self.metafile = FILESTARTP + 'pliovar_metadata_global_02102019.csv'
        self.pifile = FILESTARTP + 'modeloutput_pliovar.xls'
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
        
        # put the temperature anomalies into latitude bounds
        self.boundmin = -90. + (np.arange(0, 12) * 15.)
        self.boundmax = -75. + (np.arange(0, 12) * 15)
       
        boundtemp, nval = self.put_data_to_bounds(self.sitetemp.values - self.temppi)
        boundtemp_bs, nval_bs = self.put_data_to_bounds(self.sitetemp_bs.values - self.temppi)
       

        
        return boundtemp, boundtemp_bs, self.boundmin, self.boundmax, nval, nval_bs
    
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
    
    
    def put_data_to_bounds(self, tanom):
       """
       we now have the longitude, latitude and temperature of each datapoint
       we now put them into 15deg latitude bounded regions (defined by self.bounds)
       and find the average temperature in each region
       also return the number of points in each region
       """  
      
       boundtemp = np.zeros(12)
       count_boundtemp = np.zeros(12)
      
       for i in range(0, self.nsites):    
           # if temperature is a number add the temperature to the 
           # bound region
         
           if np.isfinite(tanom[i]):
               for bound in range(0,12):
              
                   if ((self.boundmin[bound] < self.lat[i]) &
                       (self.lat[i] <= self.boundmax[bound])):
            
                           boundtemp[bound] = (boundtemp[bound] + tanom[i])
                           count_boundtemp[bound] = count_boundtemp[bound] + 1
                           
                   
       # get average
      
       boundtemp = boundtemp / count_boundtemp
       
       print('bound temp',boundtemp)
       print('nbound',count_boundtemp)
       print(self.boundmin)
      
      
       return boundtemp, count_boundtemp
       
        
    
# end of class   
 
######################################################################
def combine_mgca_uk37(temp_uk37, temp_mgca, n_uk37, n_mgca,
                     boundmin_uk37, boundmax_uk37, boundmin_mgca, boundmax_mgca):
    """
    we have lots of arrays which show average temperature (temp_????) in a bounded
    box (boundmin_???? - boundmax_????).  We also have the number of points that make
    up the average temperature.
    
    We would like to average these depending on how many points are in the box
    
    returns: an array which contains the average combined temperature of uk37 and mg/ca
    """
    
    nvals = len(boundmin_uk37)
    temp_combined = np.zeros(nvals)
        
    for i in range(0, nvals):
        j = np.where(boundmin_mgca == boundmin_uk37[i])
        
        print(i, n_uk37[i], n_mgca[i])
        temp_combined[i] = (((np.nan_to_num(temp_uk37[i]) * n_uk37[i]) +
                            (np.nan_to_num(temp_mgca[j]) * n_mgca[j])) / 
                            (n_uk37[i] + n_mgca[j]))
    
    return temp_combined
 
#####################################################
def weightdata(tanom, lowerbounds, upperbounds):
        """
        we have some proxy data.  We would like to weight it to find the
        average value in each region.  
        Currently the proxy data is in latitudinal bands
        we want to find the global average and the (30S-30N) - (60N-75N)
        gradients
        
        input: temperature anomaly data, 
        lowervalue of the latituinal bound
        uppervalue of the latitudinal bound
        
        returns: merid_dy (30S-30N) - (60N-75N) - weighted by ocean area
                 glob_avg (average of (60S-75N) - weighted by ocean area)
        """
       
        # get land sea mask from pi model
        lsmcube = iris.load_cube(FILESTART + 'HadCM3/E280.SST.allmean.nc')
        lsmcube.data = lsmcube.data / lsmcube.data
        # get weights and multiply by lsm
        lsmcube.coord('latitude').guess_bounds()
        lsmcube.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(lsmcube)
        lsmcube.data = lsmcube.data * grid_areas
       
        
        # find weights area between each lowerbound and upper bound
        weights = np.zeros(len(tanom))
        for i, lat in enumerate(lsmcube.coord('latitude').points):
            latcube = lsmcube.extract(iris.Constraint(latitude=lat))
            latdata = np.sum(latcube.data)
            for bound in range(0, len(lowerbounds)):
                if lowerbounds[bound] < lat <=upperbounds[bound]:
                   weights[bound] = weights[bound] + latdata
        
        # global average = temperautre in each bound * weights / total of weights

        
        globavg = np.nansum(weights * tanom) / np.nansum(weights)
        print('weights and tanom')
        print('weights',weights)
        print('tanom',tanom)
        print('upperbound',upperbounds)
        
        
        # get average in each bound
        
        tropavg_deep = bound_avg(-15., 15., lowerbounds, upperbounds, weights, tanom)
        tropavg = bound_avg(-30., 30., lowerbounds, upperbounds, weights, tanom)
        tropavg_nh = bound_avg(0., 30.,lowerbounds, upperbounds, weights, tanom)
        highlatavg = bound_avg(60., 75.,lowerbounds, upperbounds, weights, tanom)
        
        # get average polewards of 45degrees
        weightcount_midhighlat = 0.
        midhighlatavg = 0.
        for bound in range(0, len(lowerbounds)):
            # do poleward of 45 deg
            if ((lowerbounds[bound] >=45. or upperbounds[bound] <=-45.)
                 and (np.isfinite(tanom[bound]))):
                midhighlatavg = midhighlatavg + (weights[bound] * tanom[bound])
                weightcount_midhighlat = weightcount_midhighlat + weights[bound]
        
        midhighlatavg = midhighlatavg / weightcount_midhighlat
        
        # get gradients
        merid_dy = tropavg - highlatavg
        
        merid_dy_nh = tropavg_nh - highlatavg
        
        merid_dy_45 = tropavg_deep - midhighlatavg
        
    
        return merid_dy, merid_dy_nh, merid_dy_45, globavg
    
def bound_avg(minval, maxval, lowerbounds, upperbounds, weights, tanom):
    """
    averages the temperature over the range minval, maxval
    returns:  weighted averaged temperautre
    """
    
    avgval = 0.
    weightcount = 0.
      
    for bound in range(0, len(lowerbounds)):
        if lowerbounds[bound] >= minval and upperbounds[bound] <= maxval:
            avgval = avgval + (weights[bound] * tanom[bound]) 
            weightcount = weightcount + weights[bound]
    avgval = avgval / weightcount
    
    return avgval
            
############################################################
def get_data(filereq, field, modeluse):
    """
    gets the field (field) from the file (filereq) and loads it
    into an iris cube (the model name is in modeluse)
    outputs a cube of the data that is as simple as possible
    """

    if modeluse == 'MMM':
        cube = iris.load_cube(filereq, field)
    else:
        cubes = iris.load(filereq)
        cube = cubes[0]
    cube.data = cube.data.astype('float32')

    if field == 'SST' or field == 'NearSurfaceTemperature':
        if (modeluse == 'MIROC4m' or modeluse == 'COSMOS'
                or modeluse == 'CESM1.0.5'):
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

    if modeluse == 'EC-Earth3.1' and field == 'SST':
        cube.coord('latitude').bounds = None
        cube.coord('longitude').bounds = None

    cube.cell_methods = None

    return cube

def globmean(cube):
    """
    returns the global mean value of a SST cube (single value)
    """

    cube.coord('latitude').guess_bounds()
    cube.coord('longitude').guess_bounds()
    grid_areas = iris.analysis.cartography.area_weights(cube)
    tempcube = cube.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas)

    return tempcube.data

def getgradient(cube):
    """
    gets the gradient in the cube.  This is the average SST 60-75N
    minus the average SST equatorward of 30N
    input : cube
    output : gradient_et (30N-30S) - (60N-75N)
             gradient_na (30N-30S) - (60N-75N)(290E-5E)
             gradient_na_nh (30N-0S) - (60N-75N)(290E-5E)
    """

    grid_areas = iris.analysis.cartography.area_weights(cube)
    grid_areas_tropics = np.zeros(np.shape(grid_areas))
    grid_areas_deeptropics = np.zeros(np.shape(grid_areas))
    grid_areas_tropics_nh = np.zeros(np.shape(grid_areas))
    grid_areas_et = np.zeros(np.shape(grid_areas)) # 60N-75N
    grid_areas_na = np.zeros(np.shape(grid_areas)) # 60N - 75N 290E-360E
    grid_areas_45 = np.zeros(np.shape(grid_areas))

    nlat = len(cube.coord('latitude').points)
    nlon = len(cube.coord('longitude').points)
    lats = cube.coord('latitude').points
    lons = cube.coord('longitude').points

    for j in range(0, nlat):
        if ((lats[j] >= 60.) and (lats[j] <= 75)):
            grid_areas_et[j, :] = grid_areas[j, :]
            for i in range(0, nlon):
                if (lons[i] >= 290 or lons[i]<=5.):
                    grid_areas_na[j, i] = grid_areas[j, i]
                else:
                    grid_areas_na[j, i] = 0.0
        else:
            grid_areas_et[j, :] = 0.0
            grid_areas_na[j, :] = 0.0

        if ((lats[j] <= 30.) and (lats[j] >= -30.)):
            grid_areas_tropics[j, :] = grid_areas[j, :]
        else:
            grid_areas_tropics[j, :] = 0.0
            
        if ((lats[j] <= 15.) and (lats[j] >= -15.)):
            grid_areas_deeptropics[j, :] = grid_areas[j, :]
        else:
            grid_areas_deeptropics[j, :] = 0.0
            
        if ((lats[j] >= 45.) or (lats[j] <= -45.)):
            grid_areas_45[j, :] = grid_areas[j, :]
        else:
            grid_areas_45[j, :] = 0.0
            
        if ((lats[j] <= 30.) and (lats[j] >= 0.)):
            grid_areas_tropics_nh[j, :] = grid_areas[j, :]
        else:
            grid_areas_tropics_nh[j, :] = 0.0

    temptrop = cube.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_tropics)
    temptropdeep = cube.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_deeptropics)
    temptrop_nh = cube.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas_tropics_nh)
    tempet = cube.collapsed(['latitude', 'longitude'],
                            iris.analysis.MEAN, weights=grid_areas_et)
    temp45 = cube.collapsed(['latitude', 'longitude'],
                            iris.analysis.MEAN, weights=grid_areas_45)
    tempna = cube.collapsed(['latitude', 'longitude'],
                            iris.analysis.MEAN, weights=grid_areas_na)

    gradient_et = temptrop.data - tempet.data 
    gradient_na = temptrop.data - tempna.data 
    gradient_nh = temptrop_nh.data - tempna.data
    gradient_deep_extra = temptropdeep.data - temp45.data
    
    print(np.sum(grid_areas))  
    print(np.sum(grid_areas_tropics)/np.sum(grid_areas))  
    print(np.sum(grid_areas_deeptropics)/np.sum(grid_areas))  
    print(np.sum(grid_areas_et)/np.sum(grid_areas))  
    print(np.sum(grid_areas_na)/np.sum(grid_areas)) 
    print(np.sum(grid_areas_45)/np.sum(grid_areas)) 
   
    
    return gradient_et, gradient_na, gradient_nh, gradient_deep_extra

def get_model_data(modelname):
    """
    1. gets the pliocene and the preindustrial SST data for each file
    2. calculates the global SSTA
    3. calculates the gradient as (SST polewards of 60deg) - (SST equatorward of 30deg)
    4. calculates the mPWP - PI gradient (gradient produced in 3.)

    input modelname
    output modTanom = the global mPWP-PI SSTA
           mod_gradanom = the mPWP (meridional SST gradient) minus the PI (meridional SST gradient)

    """

    #1. get data
    cubepi = get_data(FILESTART + modelname + '/E280.' + FIELDNAME + '.allmean.nc',
                      FIELDNAME, modelname)
    cubeplio = get_data(FILESTART + modelname + '/EOI400.' + FIELDNAME + '.allmean.nc',
                        FIELDNAME, modelname)
    #2 global mean anomaly
    meanpi = globmean(cubepi)
    meanplio = globmean(cubeplio)
    tempSSTA = meanplio - meanpi

    # calculate gradient
    gradpi_et, gradpi_na, gradpi_nh, gradpi_15NS_45NS = getgradient(cubepi)
    gradplio_et, gradplio_na, gradplio_nh, gradplio_15NS_45NS = getgradient(cubeplio)
    gradSSTA_et = gradplio_et - gradpi_et
    gradSSTA_na = gradplio_na - gradpi_na
    gradSSTA_nh = gradplio_nh - gradpi_nh
    gradSSTA_15NS_45NS = gradplio_15NS_45NS - gradpi_15NS_45NS

    return tempSSTA, gradSSTA_et, gradSSTA_na, gradSSTA_nh, gradSSTA_15NS_45NS

def plotdata(global_data, merid_grad_data, Tanom_model, gradanom_model, fileout):
    """
    this will plot the data and the model output
    we are plotting the global temerature anomaly vs the meridional gradient
    temperautre anomaly for both model and data
    """
    
    plt.scatter(global_data, merid_grad_data,
                label='proxy data', color='black',
                marker = 's')
    
    for i, model in enumerate(MODELNAMES):
         # add this to scatterplot
        if i < 6:
            plt.scatter(Tanom_model[i], gradanom_model[i], label=model)
        else:
            plt.scatter(Tanom_model[i], gradanom_model[i], label=model, marker='^')

    plt.xlabel('global mean SSTA')
    plt.ylabel('change in meridional gradient')
    plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
    plt.tight_layout()

   
    plt.savefig(fileout + 'pdf')
    plt.savefig(fileout + 'png')
    plt.close()
    
    # write data to textfile
    txtout = open(fileout + 'txt', 'w+')

    txtout.write('model, global_mean, merid_gradient \n')
    txtout.write('proxy_data,' + 
                 np.str(np.around(global_data, 2)) + ',' + 
                 np.str(np.around(merid_grad_data,2)) + '\n')
    for i, model in enumerate(MODELNAMES):
        txtout.write(model + ',' + 
                     np.str(np.around(Tanom_model[i],2)) + ',' + 
                     np.str(np.around(gradanom_model[i],2)) + '\n')
    txtout.write('multimodelmean,' + 
                     np.str(np.around(np.mean(Tanom_model),2)) + ',' + 
                     np.str(np.around(np.mean(gradanom_model),2)) + '\n')
    txtout.close


def main():
    """
    1. Call a program that will get the data to plot
    2. Weight the proxy data by area of each latitude banc
    3. Get the output (global SST, meridional gradient SST anomaly)
       for each of the models
    4. Plot the data on a symplot
    """

    # 1. get data
    obj = GetProxy('t1', 'UK37') # get data for t1 timeslice
    (t1_temp, t1_temp_bs, 
     boundmin, boundmax,
     nval, nval_bs) = obj.get_proxydata()
    
    
    if MG_CA == 'y': # also get Mg/Ca data
        obj = GetProxy('t1', 'MGCA') # get data for t1 timeslice
        (t1_temp_mgca, t1_temp_bs_mgca, 
         boundmin_mgca, boundmax_mgca,
         nval_mgca, nval_bs_mgca)  = obj.get_proxydata()
       
        print('doing t1_comb')
        t1_comb = (combine_mgca_uk37(t1_temp, t1_temp_mgca, nval, nval_mgca,
                                     boundmin, boundmax, boundmin_mgca, boundmax_mgca))
        print('doing t1_comb_bs')
        t1_comb_bs = (combine_mgca_uk37(t1_temp_bs, t1_temp_bs_mgca, nval_bs, nval_bs_mgca,
                                     boundmin, boundmax, boundmin_mgca, boundmax_mgca))
        
        t1_temp = t1_comb
        t1_temp_bs = t1_comb_bs
        print('about to exit')
        
     
    
    # 2. weight proxydata (standard and bayspline)
    
    (merid_grad, merid_grad_nh, 
     merid_grad_15_45, global_ssta) = weightdata(t1_temp, boundmin, boundmax)
    
    (merid_grad_bs, merid_grad_nh_bs, 
     merid_grad_15_45_bs, global_ssta_bs) = weightdata(t1_temp_bs, boundmin, boundmax)
   
    
    # 3. get model results 
    Tanom = np.zeros(len(MODELNAMES))
    gradanom_et = np.zeros(len(MODELNAMES))
    gradanom_na = np.zeros(len(MODELNAMES))
    gradanom_na_nh = np.zeros(len(MODELNAMES))
    gradanom_15_45 = np.zeros(len(MODELNAMES))
    
    for i, model in enumerate(MODELNAMES):
        (Tanom[i], gradanom_et[i], 
         gradanom_na[i], 
         gradanom_na_nh[i],
         gradanom_15_45[i]) = get_model_data(model)
   
    #4. plot standard (gradient 30N-30S - NH north atlantic)
    fileout = OUTSTART+'PlioVAR_gradients_uk37' + MGCASS + '.'
    plotdata(global_ssta, merid_grad, Tanom, gradanom_na, fileout)
    
    #plot standard (gradient 30N-0N - NH north atlantic)
    fileout = OUTSTART+'PlioVAR_gradients_uk37' + MGCASS + '_nh.'
    plotdata(global_ssta, merid_grad_nh, Tanom, gradanom_na_nh, fileout)
    
    #plot Bayspline (gradient 30N-30S - NH north atlantic)
    fileout = OUTSTART+'PlioVAR_gradients_uk37' + MGCASS + '_bayspline.'
    plotdata(global_ssta_bs, merid_grad_bs, Tanom, gradanom_na, fileout)
    
    
    #plot Bayspline (gradient 30N-0N - NH north atlantic)
    fileout = OUTSTART+'PlioVAR_gradients_uk37' + MGCASS + '_bayspline_nh.'
    plotdata(global_ssta_bs, merid_grad_nh_bs, Tanom, gradanom_na_nh, fileout)
    
    # plot original (gradient 15N-15S, 45N/S-90N/S)
    fileout = OUTSTART+'PlioVAR_gradients_uk37' + MGCASS + '_15NS_45NS-90NS.'
    plotdata(global_ssta, merid_grad_15_45, Tanom, gradanom_15_45, fileout)
    
    # plot Bayspline (gradient 15N-15S, 45N/S-90N/S)
    fileout = OUTSTART+'PlioVAR_gradients_uk37' + MGCASS + '_15NS_45NS-90NS_bayspline.'
    plotdata(global_ssta_bs, merid_grad_15_45_bs, Tanom, gradanom_15_45, fileout)



    return

    

# variable definition
LINUX_WIN = 'l'
FIELDNAME = 'SST'

if LINUX_WIN == 'l':
    FILESTART = ('/nfs/hera1/earjcti/regridded/')
    OUTSTART = FILESTART + 'allplots/' + FIELDNAME + '/'
    FILESTARTP = '/nfs/hera1/earjcti/PLIOMIP2/proxydata/' # where proxy data is
else:
    FILESTART = 'C:/Users/julia/OneDrive/WORK/DATA/regridded/'
    OUTSTART = FILESTART + 'allplots/' + FIELDNAME + '/'
    FILESTARTP = 'C:/Users/julia/OneDrive/WORK/DATA/proxydata/'

UNITS = 'deg C'
TIMEPERIODS = ['pi', 'mPWP']
MODELNAMES = ['CCSM4', 'CCSM4-UoT',
              'CCSM4-Utr',  
              'CESM1.2','CESM2',
              'COSMOS', 'EC-Earth3.3', 
              'GISS2.1G', 'HadCM3',
              'IPSLCM6A', 'IPSLCM5A2', 'IPSLCM5A',
              'MIROC4m', 'MRI2.3',
              'NorESM-L', 'NorESM1-F'
             ]


MG_CA = 'y'
if MG_CA == 'y':
    MGCASS = '_mgca'
else:
    MGCASS = ''

main()
