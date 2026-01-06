#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:11:26 2019
@author: earjcti

 This program will get the means from the PlioMIP1 models that
 can be added to the pliomip2 figures

#
 """

#import os
import warnings
import sys
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import matplotlib.cm as cm
#from matplotlib.colors import Normalize
import numpy as np
import iris
#import xlwt
#from xlwt import Workbook

warnings.filterwarnings("ignore")


###############################################################################

class Getmodeldata:
    """
    get all of the data from the model with the required averaing
    ie e280 or eoi400
    """
    def __init__(self, field, period):

        fieldunits = {
            "SST" : "degC",
            "TotalPrecipitation" : "mm/day",
            "NearSurfaceTemperature" : "degC"
                        }


        self.fieldname = field
        self.period = period
        if period == 'E280':
            if LINUX_WIN == 'l':
                self.lsm_file = ('/nfs/a103/eardjh/Datasets/PlioMIP/' +
                                 'data/Exp2_files/lsm/PlioMIP_Ctrl_landmasks_p3grid.nc')
            else:
                self.lsm_file = FILESTART + '/landmasks/PlioMIP_Ctrl_landmasks_p3grid.nc'
        if period == 'EOI400':
            if LINUX_WIN == 'l':
                self.lsm_file = ('/nfs/a103/eardjh/Datasets/PlioMIP/' +
                                 'data/Exp2_files/lsm/PlioMIP_Plio_landmasks_p3grid.nc')
            else:
                self.lsm_file = FILESTART + '/landmasks/PlioMIP_Plio_landmasks_p3grid.nc'


        self.units = fieldunits.get(field)



    def get_fieldreq(self, model_):
        """
        All fields are in a single file.  This will get the fieldname
        required via a number of dictionaries
        """


        PeriodE280Use = {
            "CCSM" : "ctrl", "COSMOS" : "Ctrl", "GISS" : "Ctrl",
            "HAD" : "ctrl", "IPSL" : "ctrl", "MIROC" : "ctrl",
            "MRI" : "ctrl", "NOR" : "ctrl"
                        }

        PeriodEoi400Use = {
            "CCSM" : "plio", "COSMOS" : "plio", "GISS" : "Plio",
            "HAD" : "plio", "IPSL" : "plio", "MIROC" : "plio",
            "MRI" : "plio", "NOR" : "plio"
                          }

        fieldname_sst = {
            "CCSM" : "sst", "COSMOS" : "SST", "GISS" : "SST",
            "HAD" : "sst", "IPSL" : "sst", "MIROC" : "sst",
            "MRI" : "sst", "NOR" : "sst"
                        }

        fieldname_sat = {
            "CCSM" : "sat", "COSMOS" : "sat", "GISS" : "SAT",
            "HAD" : "SAT", "IPSL" : "sat", "MIROC" : "sat",
            "MRI" : "sat", "NOR" : "sat"
                        }


        if FIELDNAME == 'NearSurfaceTemperature':
            fielduse = fieldname_sat.get(model_)
        if FIELDNAME == 'TotalPrecipitation':
            fielduse = 'precip'

        if self.period == 'E280':
            fieldreq = (model_ + '_' +
                             PeriodE280Use.get(model_) + '_' + fielduse)
        if self.period == 'EOI400':
            fieldreq = (model_ + '_' +
                             PeriodEoi400Use.get(model_) + '_' + fielduse)


        return fieldreq

    def extract_cube(self, allcube_, fieldreq_):
        """
        will extract the cube from the list of cubes
        this is needed because the cube comes from varname not long name
        """
        ncubes = len(allcube_)

        # look for the cube
        cube_found = 'n'
        for i in range(0, ncubes):
            if allcube_[i].var_name == fieldreq_:
                cube = allcube_[i]
                cube_found = 'y'

        # if cube not found then change the capitalisation of the first
        # letter of the time period (ie plio ==> Plio,  Plio ==> plio)
        if cube_found == 'n':
            newstring = ' '
            if fieldreq_.find("Plio") > 0:
                newstring = fieldreq_.replace("Plio", "plio")
            if fieldreq_.find("plio") > 0:
                newstring = fieldreq_.replace("plio", "Plio")
            if fieldreq_.find("ctrl") > 0:
                newstring = fieldreq_.replace("ctrl", "Ctrl")
            if fieldreq_.find("Ctrl") > 0:
                newstring = fieldreq_.replace("Ctrl", "ctrl")
            for i in range(0, ncubes):
                if allcube_[i].var_name == newstring:
                    cube = allcube_[i]
                    cube_found = 'y'

        cube = iris.util.squeeze(cube)

        return cube

    def get_lsm(self, modname, datacube):
        """
        gets the lsm from the lsm file
        """

        alllsm = iris.load(self.lsm_file)
        ncubes = len(alllsm)
        fieldreq = modname + '_landmask'
        cube_found = 'n'

        time = {'E280': 'Ctrl',
                'EOI400' : 'Plio'}
        for i in range(0, ncubes):
            if alllsm[i].var_name == fieldreq:
                cube = alllsm[i]
                cube = iris.util.squeeze(cube)
                cube_found = 'y'

        # if cube not found then add the time period before the landseamask
        if cube_found == 'n':
            newstring = fieldreq.replace("landmask",
                                         time.get(self.period) + "_landmask")
            for i in range(0, ncubes):
                if alllsm[i].var_name == newstring:
                    cube = alllsm[i]
                    cube_found = 'y'

        if cube_found == 'n': #cube still not found then error
            print('cannot find lsm ' + fieldreq + ' or ' + newstring)
            sys.exit(0)

        # check coords of lsm cube are the same as the
        # coordinates of the data cube
        if np.array_equal(datacube.coord('longitude').points,
                          cube.coord('longitude').points):
            pass

        else:
            print('longitudes dont match')
            sys.exit(0)

        if np.array_equal(datacube.coord('latitude').points,
                          cube.coord('latitude').points):
            pass
        else:
            print('latitudes dont match')
            sys.exit(0)


        return np.squeeze(cube.data)

    def get_globalmean(self, cube):
        """
        calculates the area weighted global mean from the cube given
        """
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube)

        meancube = cube.collapsed(['longitude', 'latitude'],
                                  iris.analysis.MEAN,
                                  weights=grid_areas)
        if FIELDNAME == 'NearSurfaceTemperature' and meancube.data > 270.:
            meancube.data = meancube.data - 273.15

        return meancube.data, grid_areas

    def get_latmean(self, cube, grid_areas, nbands):
        """
        calculates the area weighted mean in latitude bounds
        from the cube given
        """

        model_latbands = np.zeros((nbands))
        for boundno, thisbound in enumerate(LATBANDS):
            grid_areas_thisbound = np.zeros((np.shape(grid_areas)))
            for j, lat in enumerate(cube.coord('latitude').points):
                if thisbound[0] <= lat < thisbound[1]:
                    if cube.ndim == 3:
                        grid_areas_thisbound[:, j, :] = grid_areas[:, j, :]
                    if cube.ndim == 4:
                        grid_areas_thisbound[:, :, j, :] = grid_areas[:, :, j, :]
                    if cube.ndim == 2:
                        grid_areas_thisbound[j, :] = grid_areas[j, :]

            meancube = cube.collapsed(['longitude', 'latitude'],
                                      iris.analysis.MEAN,
                                      weights=grid_areas_thisbound)

            if FIELDNAME == 'NearSurfaceTemperature' and meancube.data > 250.:
                meancube.data = meancube.data - 273.15

            model_latbands[boundno] = meancube.data


        return model_latbands

    def get_land_sea_means(self, cube, grid_areas, modname):
        """
        extract the land and the sea anomaly
        """

        # get lsm
        lsm = self.get_lsm(modname, cube)


        # get grid_areas (_20 means from 20N-20S)
        grid_areas_land = grid_areas * lsm
        grid_areas_sea = grid_areas - grid_areas_land
        grid_areas_land_20 = np.zeros(np.shape(grid_areas))
        grid_areas_sea_20 = np.zeros(np.shape(grid_areas))

        for j, lat in enumerate(cube.coord('latitude').points):
            if -20 <= lat <= 20:
                grid_areas_land_20[j, :] = grid_areas_land[j, :]
                grid_areas_sea_20[j, :] = grid_areas_sea[j, :]

        mean_land = cube.collapsed(['longitude', 'latitude'],
                                   iris.analysis.MEAN,
                                   weights=grid_areas_land)
        mean_sea = cube.collapsed(['longitude', 'latitude'],
                                  iris.analysis.MEAN,
                                  weights=grid_areas_sea)
        mean_land_20 = cube.collapsed(['longitude', 'latitude'],
                                      iris.analysis.MEAN,
                                      weights=grid_areas_land_20)
        mean_sea_20 = cube.collapsed(['longitude', 'latitude'],
                                     iris.analysis.MEAN,
                                     weights=grid_areas_sea_20)

        if FIELDNAME == 'NearSurfaceTemperature' and mean_land.data > 250.:
            mean_land.data = mean_land.data - 273.15
            mean_sea.data = mean_sea.data - 273.15
            mean_land_20.data = mean_land_20.data - 273.15
            mean_sea_20.data = mean_sea_20.data - 273.15

        return (mean_land.data, mean_sea.data, mean_land_20.data,
                mean_sea_20.data)

    def get_highlatmeans(self, cube, grid_areas, latval):
        """
        gets the means polewards of latval(north) and polewards of latval(S)
        """
        grid_areas_NH = np.zeros(np.shape(grid_areas))
        grid_areas_SH = np.zeros(np.shape(grid_areas))
        for j, lat in enumerate(cube.coord('latitude').points):
            if lat >= latval:
               grid_areas_NH[j, :] = grid_areas[j, :]
            if lat <= (-1.0) * latval:
                grid_areas_SH[j, :] = grid_areas[j, :]
                
        mean_NH = cube.collapsed(['longitude', 'latitude'],
                                   iris.analysis.MEAN,
                                   weights=grid_areas_NH)
        mean_SH = cube.collapsed(['longitude', 'latitude'],
                                   iris.analysis.MEAN,
                                   weights=grid_areas_SH)

        return (mean_NH.data, mean_SH.data)

    def extract_means(self):
        """
        the top level function for this class.  This will return all the
        means to the main program
        """
        # arrays to store data
        meanvals = np.zeros(len(MODELNAMES))
        meanvals_land = np.zeros(len(MODELNAMES))
        meanvals_sea = np.zeros(len(MODELNAMES))
        meanvals_land_20 = np.zeros(len(MODELNAMES))
        meanvals_sea_20 = np.zeros(len(MODELNAMES))
        meanvals_NH45 = np.zeros(len(MODELNAMES))
        meanvals_SH45 = np.zeros(len(MODELNAMES))
        meanvals_NH60 = np.zeros(len(MODELNAMES))
        meanvals_SH60 = np.zeros(len(MODELNAMES))

        nbands, nlims = np.shape(LATBANDS)
        latmeans = np.zeros((len(MODELNAMES), nbands))

        filename = (FILESTART + 'PlioMIP1_regridded.nc')
        allcubes = iris.load(filename)

        for i, model in enumerate(MODELNAMES):


            fieldreq = self.get_fieldreq(model)

            # extract the cube we want
            cube = self.extract_cube(allcubes, fieldreq)
            meanvals[i], grid_areas = self.get_globalmean(cube)
            latmeans[i, :] = self.get_latmean(cube, grid_areas, nbands)

            (meanvals_land[i],
             meanvals_sea[i],
             meanvals_land_20[i],
             meanvals_sea_20[i]) = self.get_land_sea_means(cube, grid_areas, model)
            
            (meanvals_NH45[i], 
             meanvals_SH45[i]) = self.get_highlatmeans(cube,grid_areas, 45.0)
            
            (meanvals_NH60[i], 
             meanvals_SH60[i]) = self.get_highlatmeans(cube,grid_areas, 60.0)


        return (meanvals, latmeans, meanvals_land,
                meanvals_sea, meanvals_land_20, meanvals_sea_20,
                meanvals_NH45, meanvals_SH45,
                meanvals_NH60, meanvals_SH60)




# end of class Getmodeldata
##################################################################################
def get_nh_seascyc():
    """
    get the NH seasonal cycle from each model
    """


    shortfield = {"NearSurfaceTemperature" : "SAT",
                  "TotalPrecipitation" : "precip"}
    longfield_sat = {"CCSM" : "Reference height temperature",
                     "COSMOS" : "2m temperature",
                     "IPSL" : "t2m",
                     "MIROC": "tas",
                     "MRI" : "near surface air temperature [degC]"
                     }

    longfield_precip = {
        "CCSM" : "Total (convective and large-scale) precipitation rate (liq + ice)",
        "COSMOS" : "total precipitation",
        "IPSL" : "IPSL_precip",
        "MIROC": "MIROC_precip",
        "MRI" : "total precipitation [mm/day]"
                       }
    nh_means = np.zeros((len(MODELNAMES), 12))

    for i, model in enumerate(MODELNAMES):
        filename = (FILESTART + '/' + FIELDNAME + '/' +
                    model + '_Exp2_anom_' +
                    shortfield.get(FIELDNAME) + '_p3grid.nc')

        if FIELDNAME == "NearSurfaceTemperature":
            fieldreq = longfield_sat.get(model, "air_temperature")
        if FIELDNAME == "TotalPrecipitation":
            fieldreq = longfield_precip.get(model, "precipitation_flux")

        cube = iris.load_cube(filename, fieldreq)
        cube = iris.util.squeeze(cube)

        # get grid areas for seasonal average
        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()
        grid_areas = iris.analysis.cartography.area_weights(cube)
        grid_areas_nh = np.zeros(np.shape(grid_areas))

        for j, lat in enumerate(cube.coord('latitude').points):
            if lat > 0:
                grid_areas_nh[:, j, :] = grid_areas[:, j, :]

        # get seasonal average
        meancube = cube.collapsed(['longitude', 'latitude'],
                                  iris.analysis.MEAN,
                                  weights=grid_areas_nh)

        nh_means[i, :] = meancube.data
        plt.plot(meancube.data, label=model)
    plt.plot(np.mean(nh_means, axis=0), label='mean')
    plt.legend()
    #plt.show()

    return nh_means


def write_global_means(filetext, modelmean_eoi400, modelmean_e280):
    """
    write global means from each model to a text file

    """

    modelmean_anomaly = (modelmean_eoi400 -modelmean_e280)

    filetext.write("modelname, global mean EOI400," +
                   "global mean E280, anomaly \n")

    for i, model in enumerate(MODELNAMES_FULL):
        eoi400 = np.str(np.around(modelmean_eoi400[i], 3))
        e280 = np.str(np.around(modelmean_e280[i], 3))
        anom = np.str(np.around(modelmean_anomaly[i], 3))

        filetext.write((model + ',' + eoi400 + ',' + e280 +
                        ',' + anom + '\n'))
    filetext.write('MEAN,' +
                   np.str(np.around(np.mean(modelmean_eoi400), 3)) +
                   ',' +
                   np.str(np.around(np.mean(modelmean_e280), 3)) +
                   ',' +
                   np.str(np.around(np.mean(modelmean_anomaly), 3)) + '\n')



def write_lat_means(filetext, modelmean_eoi400, modelmean_e280):
    """
    write the latitudinal mean from each model to a textfile
    """


    modelmean_anomaly = (modelmean_eoi400 -modelmean_e280)
    mean_eoi400 = np.mean(modelmean_eoi400, axis=0)
    mean_e280 = np.mean(modelmean_e280, axis=0)
    mean_anomaly = mean_eoi400 - mean_e280

    filetext.write("modelname, latband mean EOI400," +
                   "latband mean E280, latband anomaly \n")
    filetext.write("bands are " + np.str(LATBANDS) + '\n')

    for i, model in enumerate(MODELNAMES_FULL):
        eoi400 = np.str(np.around(modelmean_eoi400[i], 3))
        e280 = np.str(np.around(modelmean_e280[i], 3))
        anom = np.str(np.around(modelmean_anomaly[i], 3))

        filetext.write((model + ',' + eoi400 + ',' + e280 +
                        ',' + anom + '\n'))
    filetext.write('MEAN,' +
                   np.str(np.around(mean_eoi400, 3)) + ',' +
                   np.str(np.around(mean_e280, 3)) +',' +
                   np.str(np.around(mean_anomaly, 3)) + '\n')

def write_nh_seascycle(filetext, seasanom):
    """
    write the NH seasonal anomaly from each model to a textfile
    """


    filetext.write("modelname, [jan feb mar apr may jun jul aug sep oct nov dec] \n")
    mean_anomaly = np.mean(seasanom, axis=0)

    for i, model in enumerate(MODELNAMES_FULL):
        anom = np.str(np.around(seasanom[i], 3))

        filetext.write((model + ',' + anom + '\n'))

    filetext.write('MEAN,' +  np.str(np.around(mean_anomaly, 3)) + '\n')

def write_landsea_means(filetext, landmean_eoi400_allmodels,
                        seamean_eoi400_allmodels, landmean_e280_allmodels,
                        seamean_e280_allmodels, region):
    """
    write the land sea averages to a text file
    """
    filetext.write("modelname" + region + "[mean_ocean_eoi400, meanocean_e280, meanocean_anom, " +
                   "mean_land_eoi400, mean_land_e280, mean_land_anom ] \n")

    eoi400_sea_mean = np.str(np.around(np.mean(seamean_eoi400_allmodels), 3))
    e280_sea_mean = np.str(np.around(np.mean(seamean_e280_allmodels), 3))
    eoi400_land_mean = np.str(np.around(np.mean(landmean_eoi400_allmodels), 3))
    e280_land_mean = np.str(np.around(np.mean(landmean_e280_allmodels), 3))
    sea_anom_mean = np.str(np.around(np.mean(seamean_eoi400_allmodels)
                                     - np.mean(seamean_e280_allmodels), 3))
    land_anom_mean = np.str(np.around(np.mean(landmean_eoi400_allmodels)
                                      - np.mean(landmean_e280_allmodels), 3))

    for i, model in enumerate(MODELNAMES_FULL):
        eoi400_sea = np.str(np.around(seamean_eoi400_allmodels[i], 3))
        e280_sea = np.str(np.around(seamean_e280_allmodels[i], 3))
        eoi400_land = np.str(np.around(landmean_eoi400_allmodels[i], 3))
        e280_land = np.str(np.around(landmean_e280_allmodels[i], 3))
        sea_anom = np.str(np.around((seamean_eoi400_allmodels[i]
                                     - seamean_e280_allmodels[i]), 3))
        land_anom = np.str(np.around((landmean_eoi400_allmodels[i]
                                      - landmean_e280_allmodels[i]), 3))


        filetext.write(model + ',' + eoi400_sea + ',' + e280_sea + ',' + sea_anom +
                       ',' + eoi400_land + ',' + e280_land + ',' + land_anom + '\n')

    filetext.write('MEAN,' + eoi400_sea_mean + ',' + e280_sea_mean + ',' + sea_anom_mean +
                   ',' + eoi400_land_mean + ',' + e280_land_mean + ',' + land_anom_mean + '\n')

def write_hemisphere_anom(filetext, NH_anomaly45,SH_anomaly45,
                          NH_anomaly60, SH_anomaly60):
    """
    writes the anomalies polewards of 45N and 45S and 60S and 60N
    """
    filetext.write("modelname, 45N-90N_anom, 45S-90S_anom, 60N-90N_anom, " +
                   "60S-90S_anom \n")
    for i, model in enumerate(MODELNAMES_FULL):

        filetext.write(model + ',' + np.str(np.around(NH_anomaly45[i], 3)) + ',' +
                       np.str(np.around(SH_anomaly45[i],3)) + ',' +
                       np.str(np.around(NH_anomaly60[i],3)) + ',' +
                       np.str(np.around(SH_anomaly60[i],3)) + '\n')
    filetext.write('MEAN,' + np.str(np.around(np.mean(NH_anomaly45), 3)) + ',' +
                   np.str(np.around(np.mean(SH_anomaly45),3)) + ',' +
                   np.str(np.around(np.mean(NH_anomaly60),3)) + ',' +
                   np.str(np.around(np.mean(SH_anomaly60),3)) +  '\n')

def main():

    # get data means and latitudinal anomalies
    modeldata = Getmodeldata(FIELDNAME, 'EOI400')
    (globalmean_eoi400_allmodels,
     latmean_eoi400_allmodels,
     landmean_eoi400_allmodels,
     seamean_eoi400_allmodels,
     landmean20_eoi400_allmodels,
     seamean20_eoi400_allmodels,
     NH45_eoi400_anomaly, SH45_eoi400_anomaly,
     NH60_eoi400_anomaly, SH60_eoi400_anomaly) = modeldata.extract_means()

    modeldata = Getmodeldata(FIELDNAME, 'E280')
    (globalmean_e280_allmodels,
     latmean_e280_allmodels,
     landmean_e280_allmodels,
     seamean_e280_allmodels,
     landmean20_e280_allmodels,
     seamean20_e280_allmodels,
     NH45_e280_anomaly, SH45_e280_anomaly,
     NH60_e280_anomaly, SH60_e280_anomaly) = modeldata.extract_means()

    # get seasonal cycle
    NH_seas_anom = get_nh_seascyc()


    # write to text

    filetext = open((FILESTART + '/means_for_' + FIELDNAME + '.txt'), "w+")
    write_global_means(filetext, globalmean_eoi400_allmodels,
                       globalmean_e280_allmodels)
    write_lat_means(filetext, latmean_eoi400_allmodels,
                    latmean_e280_allmodels)
    write_nh_seascycle(filetext, NH_seas_anom)
    write_landsea_means(filetext, landmean_eoi400_allmodels, # global land sea
                        seamean_eoi400_allmodels, landmean_e280_allmodels,
                        seamean_e280_allmodels, "global")
    write_landsea_means(filetext, landmean20_eoi400_allmodels, # global land sea
                        seamean20_eoi400_allmodels, landmean20_e280_allmodels,
                        seamean20_e280_allmodels, "20N-20S")
    write_hemisphere_anom(filetext, NH45_eoi400_anomaly - NH45_e280_anomaly,
                          SH45_eoi400_anomaly - SH45_e280_anomaly,
                          NH60_eoi400_anomaly - NH60_e280_anomaly,
                          SH60_eoi400_anomaly - SH60_e280_anomaly)

    filetext.close()


####################################
# definitions

LINUX_WIN = 'l'
#FIELDNAME = 'TotalPrecipitation'
FIELDNAME = 'NearSurfaceTemperature'
if LINUX_WIN == 'l':
    FILESTART = '/nfs/hera1/earjcti/PLIOMIP_old/'
else:
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\PLIOMIP1\\'

MODELNAMES = ['CCSM', 'COSMOS', 'GISS', 'HAD', 'IPSL',
              'MIROC', 'MRI', 'NOR']

CONVERT_MODELS = {'CCSM' : 'CCSM4',
                  'GISS' : 'GISS-E2-R',
                  'HAD'  : 'HadCM3',
                  'IPSL' : 'IPSLCM5A',
                  'MIROC': 'MIROC4m',
                  'MRI' : 'MRI2.3',
                  'NOR' : 'NORESM-L'}

MODELNAMES_FULL = []
for MOD in MODELNAMES:
    MODELNAMES_FULL.append(CONVERT_MODELS.get(MOD, MOD))

LATBANDS = [[-90., -60.], [-60., -30.], [-30., 0.], [0., 30],
            [30., 60.], [60., 90.]]


main()

