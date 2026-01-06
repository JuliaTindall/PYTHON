#!/usr/bin/env python3
"""
#NAME
#    ITCZ diagnostics
#PURPOSESS
#    This program will find the ITCZ in a given region based on the
#    Stanfield et al 2015 definition.  It will find the centerline width
#    and intensity and plot these by season for the Pliocene and the PI.
#
# search for 'main program' to find end of functions
# Julia August 2018
#
# Notes This program is like ITCZ_diagnostics_regional.  However it allows
# the ITCZ to be discontinuous.
"""

import os
os.environ['PROJ_LIB'] = 'C:/Users/julia/Miniconda3/envs/py3/Library/share'
#os.environ['PROJ_LIB'] = '/nfs/see-fs-02_users/earjcti/anaconda2/envs/py3/share/proj'
from mpl_toolkits.basemap import Basemap
import iris
import iris.plot as iplt
import iris.analysis.cartography
from iris.experimental.equalise_cubes import equalise_attributes
import numpy as np
#import matplotlib as mp
import matplotlib.pyplot as plt
import sys


def plotmean_newaxis(cube, modelno_):
    """
    add a new axis to cube for help with concatenation
    """
    tempcube = iris.util.new_axis(cube)
    tempcube.add_dim_coord(iris.coords.DimCoord(modelno_,
                                                standard_name='model_level_number', 
                                                long_name='model',
                                                var_name='model',
                                                units=None,
                                                bounds=None,
                                                coord_system=None,
                                                circular=False), 0)
    return tempcube

def recentre(cube, upper_bound_index, lower_bound_index,
             most_intense_precip_index):
    """
    recentre so that all the indexes go from -180 to 180
    """

    lons = cube.coord('longitude').points
    lons_recentre = np.where(lons < 180., lons, lons -360.)
    indexes = lons_recentre.argsort() # of the recentered array

    lats = cube.coord('latitude').points

    lat_upper = np.full((NMONTHS, len(lons)), np.nan)
    lat_lower = np.full((NMONTHS, len(lons)), np.nan)
    lat_max = np.full((NMONTHS, len(lons)), np.nan)

    for month in range(len(MONTHNAMES)):
        for i in range(len(lons)):
            if not upper_bound_index.mask[month, i]:
                lat_upper[month, i] = lats[upper_bound_index[month, i] - 1]
            if not lower_bound_index.mask[month, i]:
                lat_lower[month, i] = lats[lower_bound_index[month, i]]
            if not most_intense_precip_index.mask[month, i]:
                lat_max[month, i] = lats[most_intense_precip_index[month, i]]

        lat_upper[month, :] = lat_upper[month, indexes]
        lat_lower[month, :] = lat_lower[month, indexes]
        lat_max[month, :] = lat_max[month, indexes]

    return(lons_recentre[indexes], lat_upper, lat_lower, lat_max)


class GetITCZ:
    """
    the class object is the region name
    this class will deal with finding the itcz in this region (N/S bounds
    and maximum intensity)
    """
    def __init__(self, region):

        # l = land, b=both if not set then it is o=ocean
        landoceanind = {"IndianOceanLand" : "b",
                        "Globe" : "b",
                        "Africa" : "l",
                        "CentralAmerica" : "l",
                        "EastAtlantic" : "b",
                        "Indian" : "l", "EastAsia" : "l", "Indonesia" : "b",
                        "IndonesiaLand" : "l", "Australia" : "l",
                        "SouthAmerica" : "l", "WesternPacific" : "b",
                        "EasternPacific" : "b", "EastPacificExt" : "b"
                        }
        latmin = {"Globe" : -30., "Africa" : -30., "CentralAmerica" : -25.,
                  "EastAtlantic" : -25., "Indian" : 10.0,
                  "EastAsia" : 0.0, "Indonesia" : -12.0,
                  "IndonesiaLand" : -10.0, "Australia" : -30.0,
                  "SouthAmerica" : -20.0, "WesternPacific" : -30.0,
                  "EasternPacific" : -30.0, "EastPacificExt" : -10.0}

        latmax = {"IndianOcean" : 30, "IndianOceanLand" : 30,
                  "Globe" : 45.0, "Africa" : 35.0,
                  "CentralAmerica" : 25.0, "EastAtlantic" : 25.0,
                  "Indian" : 30.0, "EastAsia" : 50.0, "Indonesia" : 10.0,
                  "IndonesiaLand" : 0.0, "Australia" : -12.0,
                  "SouthAmerica" : 15.0, "WesternPacific" :30.0,
                  "EasternPacific" : 30.0, "EastPacificExt" : 30.0}

        lonmin = {"AtlanticOcean" : -60.0, "IndianOcean" : 30,
                  "IndianOceanLand" : 60.0, "Globe" : 0, "Africa" : -20.0,
                  "CentralAmerica" : -120.0, "EastAtlantic" : -120.0,
                  "Indian" : 70.0, "EastAsia" : 60.0, "Indonesia" : 110.0,
                  "IndonesiaLand" : 130.0, "Australia" : 110,
                  "SouthAmerica" : -80.0, "WesternPacific" : 150.0,
                  "EasternPacific" : -120.0, "EastPacificExt" : -120.0}

        lonmax = {"AtlanticOcean" : 20.0, "IndianOcean" : 120,
                  "IndianOceanLand" : 90.0, "Globe" : 360.0, "Africa" : 45.0,
                  "CentralAmerica" : -80.0, "EastAtlantic" : -80.0,
                  "Indian" : 110.0, "EastAsia" : 150.0, "Indonesia" : 150,
                  "IndonesiaLand" : 150.0, "Australia" : 160.0,
                  "SouthAmerica" : -30.0, "WesternPacific" : 200.0,
                  "EasternPacific" : -90.0, "EastPacificExt" : -60.0}

        reg_threshold = {"IndianOceanLand" : 4.0, "Globe" : 4.0,
                         "Africa" : 2.0, "CentralAmerica" : 0.0,
                         "EastAtlantic" : 4.0, "Indian" : 0.0,
                         "EastAsia" : 2.0, "Indonesia" : 0.0,
                         "IndonesiaLand" : 0.0, "Australia" : 1.0,
                         "SouthAmerica" : 4.0, "WesternPacific" : 4.0,
                         "EasternPacific" : 4.0, "EastPacificExt" : 4.0}

        self.regionname = region
        self.land_ocean_ind = landoceanind.get(region, 'o')
        self.latmin = latmin.get(region, -20.0)
        self.latmax = latmax.get(region, 20.0)
        self.latmax = self.latmax + 0.5 # correct grid
        self.latmin = self.latmin - 0.5
        self.lonmin = lonmin.get(region)
        self.lonmax = lonmax.get(region)
        self.threshold = THRESHOLD
        self.threshold[:] = reg_threshold.get(region, THRESHOLD)


    def get_region_AOI(self):
        """

        Returns
        -------
        The area of interest incase the calling program needs it

        """

        lons = np.arange(self.lonmin, self.lonmax, 1.0)
        lons = np.where(lons < 180., lons, lons-360.)
        return(np.min(lons), np.max(lons), self.latmin, self.latmax)


    def read_cube(self, period, model, modelno):
        """
        reads the cube and extracts the required region and equalises all the
        attributes
        """

        if period == 'E280':
            cubelsmo = iris.load_cube(E280LSM)
        if period == 'Eoi400':
            cubelsmo = iris.load_cube(EOI400LSM)

        cubegrid = iris.load_cube('one_lev_one_deg.nc')
        cube_lsm = cubelsmo.regrid(cubegrid, iris.analysis.Linear())

        # get precipitation data in mm/day

        filename = (FILESTART + 'regridded/' + model + '/' + period  +
                    '.TotalPrecipitation.mean_month.nc')
        cube_precip = iris.load_cube(filename)


        # mask out land or ocean as appropriate.
        if np.array_equal(cube_lsm.coord('longitude').points,
                          cube_precip.coord('longitude').points):
            if self.land_ocean_ind == 'l':
                cube_precip.data = cube_precip.data * cube_lsm.data
            if self.land_ocean_ind == 'o':
                cube_precip.data = cube_precip.data * np.abs(cube_lsm.data - 1.0)
        else:
            print('error lon/lat of land sea mask dont match')
            print('lsm long', cube_lsm.coord('longitude'))
            print('precip long', cube_precip.coord('longitude'))
            sys.exit()


        # decompose grid to that defined by the region

        sample_points = [('longitude',
                          np.arange(self.lonmin, self.lonmax + 1.0, 1.0)),
                         ('latitude',
                          np.arange(self.latmin, self.latmax + 1.0, 1.0))]

        cube_precip.coord('latitude').var_name = 'latitude'
        cube_precip.coord('longitude').var_name = 'longitude'
        cube_precip.coord('latitude').long_name = None
        cube_precip.coord('longitude').long_name = None

        cube_area = cube_precip.interpolate(sample_points, iris.analysis.Linear())
        cube_area.var_name = "precip"
        cube_area.standard_name = None
        cube_area.long_name = "precip"
          # remove auxillary coordinate year and time
        print(model, period, cube_area.aux_coords)
        for coord in cube_area.aux_coords:
            if coord.standard_name == None:
                name = coord.var_name
            else:
                name = coord.standard_name
            cube_area.remove_coord(name)
        # remove attributes from month
        cube_area.coord('month').attributes = None

        for coord in cube_area.coords():
            coord.points = coord.points.astype('float32')
        cube_area.data = cube_area.data.astype('float32')

        print(cube_area.coord('month').attributes)
        # add scalar coordinate for merge
        new_cube = plotmean_newaxis(cube_area, modelno)
        new_cube.cell_methods = None
        
        return new_cube


    def find_ITCZ(self, period, model, modelno):
        """
        #  to find the ITCZ we find precipitation in the given region.
        #  we then find the longest continuous stretch of precipitation above a certain
        #  threshold for each longitude.
        """

        print('in find_ITCZ')


        cube_AOI = self.read_cube(period, model, modelno)

        #


        # for each longitude find the longest continuous band of precipitation
        # above the threshold

        (upper_bound_index,
         lower_bound_index,
         most_intense_precip_index) = self.get_indexes(cube_AOI)



        # plot precipitation map for each month all on one page with bands shown
        # return some of the things we are plotting

        (lons_recentre,
         lat_upper,
         lat_lower,
         lat_max) = recentre(cube_AOI, upper_bound_index,
                             lower_bound_index,
                             most_intense_precip_index)



        return (cube_AOI, lons_recentre, lat_upper, lat_lower,
                lat_max)


    def get_indexes(self, cube):
        """
        the main part of the class which actually gets the indexes
        """

        def correct_undefined():
            """
            # where ITCZ is undefined set to NAN
            # if lower and upper are same set to undefined
            """
            nlats = len(lats)
            for j, upper_int in enumerate(upper_bound_index):
                for i, upper in enumerate(upper_int):
  #                  if upper == lower_bound_index[j,i] and upper == 0.0:
                    if upper == lower_bound_index[j, i]:
                        upper_bound_index[j, i] = np.ma.masked
                        lower_bound_index[j, i] = np.ma.masked
                        most_intense_precip_index[j, i] = np.ma.masked
                    if upper == 0 or lower_bound_index[j, i] == 0:
                        upper_bound_index[j, i] = np.ma.masked
                        lower_bound_index[j, i] = np.ma.masked
                        most_intense_precip_index[j, i] = np.ma.masked
                    if upper >= nlats-1 or lower_bound_index[j, i] >= nlats-1:
                        upper_bound_index[j, i] = np.ma.masked
                        lower_bound_index[j, i] = np.ma.masked
                        most_intense_precip_index[j, i] = np.ma.masked

            return upper_bound_index, lower_bound_index, most_intense_precip_index

        def correct_discontinuity(precip_AOI, bound_index, ns_ind):
            """
            # ns ind = 1.0 for upper bound, -1.0 for lower bound
            #  correct a very small (=< 4 gridboxes of longitude) discontinuity
            # a large discontinuity is > 6 degrees
            """

            for index, bound in np.ndenumerate(bound_index):
                mon = index[0]
                i = index[1]
                if bound_index.mask[mon, i]:
                    break

                bound_prev = bound_index[mon, i - 1]
                bound_diff = np.abs(bound - bound_prev)

                if bound_diff > 6 and bound_prev != 0:


                    # large difference will it return to previous position within
                    # n degrees (min diff = minumum difference over next 10 longitudes)
                    n = 11
                    min_diff = np.min(np.abs(
                            bound_index[mon, i + 1: i  + n]
                            - bound_prev))

                    if min_diff < bound_diff / 1.3:

                        # initially set to same as previous
                        bound_index[index] = bound_index[mon, i - 1]
                        # if new upper bound index has precipitation
                        # greater than threshold see if you can expand
                        # northwards if upper bound, southwards if lower bound
                        j1 = bound_index[index]
                        if precip_AOI[mon, bound_index[index], i] >= self.threshold[mon]:

                            j2 = len(AOI_lat)
                            if ns_ind < 0.0: j2 = 0
                            for j in range(j1, j2, ns_ind):
                                if precip_AOI[mon, j, i] >= self.threshold[mon]:
                                    bound_index[mon, i] = j # move northwards/southwards
                                else:
                                    break # no longer try and move northwards


                        # if new upper bound index has precipitation that is less
                        # then the threshold then you will have to move it
                        # southwards for the upper bound or northwards for the lower bound
                        else:
                            j2 = len(AOI_lat)
                            if ns_ind > 0: j2 = 0
                            for j in range(j1, j2, -1 * ns_ind):
                                if precip_AOI[mon, j, i] >= self.threshold[mon]:
                                    bound_index[mon, i] = j # move southwards/northwards
                                    break


            return bound_index



        def get_most_intense(precip_AOI):
            """
            gets the index of the most intense precipitation between the
            lower and upper bounds
            """
            most_intense_precip_index = np.ma.zeros(np.shape(lower_bound_index), dtype=int)
            for index, lower in np.ndenumerate(lower_bound_index):

                if (np.isfinite(lower) and np.isfinite(upper_bound_index[index])
                    and lower < upper_bound_index[index]):

                    max_val = 0.
                    for j in range(lower, upper_bound_index[index] + 1):
                        if precip_AOI[index[0], j, index[1]] > max_val:
                            max_val = precip_AOI[index[0], j, index[1]]
                            most_intense_precip_index[index] = j
                else:
                    most_intense_precip_index[index] = np.ma.masked


            return most_intense_precip_index



        print('in get_indexes')

        lons = cube.coord('longitude').points
        nlons = len(lons)
        lats = cube.coord('latitude').points
        nlats = len(lats)
        months = cube.coord('month').points
        nmonths = len(months)

        max_count_lons = np.zeros((nmonths, nlons), dtype=int)
        upper_bound_index = np.ma.zeros((nmonths, nlons), dtype=int)
        lower_bound_index = np.ma.zeros((nmonths, nlons), dtype=int)
        most_intense_precip_index = np.ma.zeros((nmonths, nlons), dtype=int)



        for mon in range(0, nmonths):
            smallcube = iris.util.squeeze(cube)
            AOI_precip = smallcube.data[mon, :, :]
            AOI_lat = cube.coord('latitude').points

            for i, lon in enumerate(lons):
                count_lons = 0
                jmax = len(lats) - 1
                jmin = 0
                # check it is not picking up mid latitude storm tracks in pacific
                # or atlantic by reducing range to -20 20 ignore for south america
                if ((lon > 180. or lon < 0) and (self.regionname != 'SouthAmerica')):
                    jmax = (np.abs(lats - 20.)).argmin()
                    jmin = (np.abs(lats + 20.)).argmin()
                # restrict to 20S in summer
                if 3 <= mon <= 8:
                    jmin = (np.abs(lats + 20.)).argmin()

                for j in range(jmax-1, jmin, -1):
                    # set up j-1 and j+1 for checking values either side
                    jmin1 = j - 1
                    jpl1 = j + 1
                    if jmin1 < 0:
                        jmin1 = 0
                    if jpl1 > nlats-1:
                        jpl1 = nlats-1
                    # check precipitation greater than threshold but include if there
                    # is only a single latitude blip

                    if (AOI_precip[j, i] >= self.threshold[mon] or
                            (AOI_precip[jmin1, i] > self.threshold[mon] and
                             AOI_precip[jpl1, i] > self.threshold[mon])):
                        count_lons = count_lons+1

                    # if the number of longitude is greater than the previous
                    # maximum then this is the new ITCZ.  However we are
                    # constraining the
                    # lower bound latitude to be less than 20N
                        if (count_lons > max_count_lons[mon, i]
                                and lats[j] < 20.):
                            max_count_lons[mon, i] = count_lons
                            lower_bound_index[mon, i] = j

                    else:
                        count_lons = 0

                #############################################
                # set up upper bound and get masimum intensity index
                upper_bound_index[mon, i] = lower_bound_index[mon, i] + max_count_lons[mon, i]


             ###########################################
            # CORRECTIONS TO THE ALGORITHM
            # make sure picks up upper bound in west pacific
            # don't use this as it is too complicated and doesn't make enough difference

            #(lower_bound_new,
            # upper_bound_new) = correct_west_pacific(lower_bound_index[mon, :],
            #                                         upper_bound_index[mon, :],
            #                                         mon)
            #lower_bound_index[mon, :] = lower_bound_new
            #upper_bound_index[mon, :] = upper_bound_new



        (upper_bound_index,
         lower_bound_index,
         most_intense_precip_index) = correct_undefined()

        upper_bound_index = correct_discontinuity(smallcube.data,
                                                  upper_bound_index, 1) # for upper bound

        lower_bound_index = correct_discontinuity(smallcube.data,
                                                  lower_bound_index, -1) # for upper bound

        most_intense_precip_index = get_most_intense(smallcube.data)


        #print('j1',upper_bound_index[1],lower_bound_index[1])
        return [upper_bound_index, lower_bound_index,
                most_intense_precip_index]


def plot_map(cubelist, lons, upper_lat, lower_lat,
             max_lat, period):
    """
    For each month plots the ICZ for each model.  Showing
    tropical precipitation and the upper and lower bounds

    """

    print('in plot map')
    fig = plt.figure(figsize=[12.0, 8.0])

    for month, monthname in enumerate(MONTHNAMES):
        fig = plt.figure(figsize=[12.0, 8.0])
        # plot all models for each month
        for i, modelname in enumerate(MODELNAMES):
            cube = cubelist[i]

            plt.subplot(4, 4, i + 1)
            cubeplot = iris.util.squeeze(cube)
            cubeplot = cubeplot[month, :, :]

            V = np.arange(0, 15, 1)
            cs = iplt.contourf(cubeplot, V, cmap='Blues', extend='max')
            plt.gca().coastlines()

            plt.plot(lons, upper_lat[i][month, :], color='red', linewidth=1)
            plt.plot(lons, lower_lat[i][month, :], color='red', linewidth=1)
            plt.plot(lons, max_lat[i][month, :], color='white', linewidth=1)
            plt.title(modelname)


        # sort out all subplots and print out to a file
        plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9, wspace=0.1, hspace=0.0)
        cb_ax = fig.add_axes([0.1, 0.1, 0.8, 0.05])
        cbarname = period + ' ' + monthname + ' precip (mm/day)'
        cbar = fig.colorbar(cs, cax=cb_ax, orientation="horizontal")
        cbar.set_label(cbarname, fontsize=15)

        fileoutstart = (FILESTART + 'ITCZ/' + REGION + '/' + period + '_'
                        + np.str(month) + monthname + '_precipitation.')

        if LINUX_WIN == 'w':
            plt.savefig(fileoutstart + 'png')
        else:
            plt.savefig(fileoutstart + 'pdf')

        plt.close()


    return

def plot_map_mmm(cubelist, lons, period):
    """
    For each month plots the multimodel mean ITCZ.  Showing
    tropical precipitation and the upper and lower bounds

    """

    print('in plot map')
    fig = plt.figure(figsize=[12.0, 8.0])

    # get average
    equalise_attributes(cubelist)
    cube = cubelist.concatenate_cube()

    cube_avg = cube.collapsed('model_level_number', iris.analysis.MEAN)

    # get indexfrom MMM
    Regioninfo = GetITCZ(REGION)

    (upper_bound_index,
     lower_bound_index,
     max_index) = Regioninfo.get_indexes(cube_avg)

    (lons_recentre, lat_upper,
     lat_lower, lat_max) = recentre(cube_avg, upper_bound_index,
                                    lower_bound_index, max_index)

    # plot
    fig = plt.figure(figsize=[12.0, 8.0])
    for month, monthname in enumerate(MONTHNAMES):

        plt.subplot(3, 4, month + 1)
        cubeplot = cube_avg[month, :, :]
        V = np.arange(0, 15, 1)
        cs = iplt.contourf(cubeplot, V, cmap='Blues', extend='max')
        plt.gca().coastlines()

        plt.plot(lons, lat_upper[month, :], color='red', linewidth=1)
        plt.plot(lons, lat_lower[month, :], color='red', linewidth=1)
        plt.plot(lons, lat_max[month, :], color='white', linewidth=1)
        plt.title(monthname)



    # sort out all subplots and print out to a file
    plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9, wspace=0.1, hspace=0.0)
    cb_ax = fig.add_axes([0.1, 0.1, 0.8, 0.05])
    cbarname = period + 'MMM precip (mm/day)'
    cbar = fig.colorbar(cs, cax=cb_ax, orientation="horizontal")
    cbar.set_label(cbarname, fontsize=15)

    fileoutstart = (FILESTART + 'ITCZ/' + REGION + '/MMM_' + period + '_'
                        + '_precipitation.')

    if LINUX_WIN == 'w':
        plt.savefig(fileoutstart + 'png')
    else:
        plt.savefig(fileoutstart + 'pdf')

    plt.close()


    # write means out so we can reaccess them
    write_info((FILESTART + 'ITCZ/' + REGION + '/data_MMM_'
                + period + '.txt'),
               lons, lat_upper, lat_lower, lat_max)

    fileout = (FILESTART + 'ITCZ/' + REGION + '_multimodelmean.nc')
    iris.save(cubeplot,fileout,netcdf_format='NETCDF3_CLASSIC',fill_value=1.0E20)



    return




def plot_lats(lons, lat_pi, lat_plio, title, lonmin, lonmax, latmin, latmax,
              filename):
    """
    will do a latitude vs longitude line on a map for both the pliocene and
    the preindustrial over 12 months
    """

    print(lons)
    print(lons[0])
    print(np.shape(lons))

    for monno, monname in enumerate(MONTHNAMES):
        fig = plt.figure(figsize=[12.0, 8.0])
        print(filename, monno)
        for i, model in enumerate(MODELNAMES):
            plt.subplot(4, 4, i+1)
            map = Basemap(llcrnrlon=lonmin, urcrnrlon=lonmax, llcrnrlat=latmin,
                          urcrnrlat=latmax, projection='cyl', resolution='c')
            map.drawcoastlines()
            plt.plot(lons[i], lat_pi[i][monno, :], label='pi')
            plt.plot(lons[i], lat_plio[i][monno, :], label='plio')
            plt.title(model)

        plt.legend()
        plt.subplots_adjust(left=0.1, bottom=0.2, right=0.9, top=0.9, wspace=0.1, hspace=0.0)
        fig.text(0.1, 0.9, title, fontsize=20)

        fileout = (FILESTART + '/ITCZ/' + REGION + '/' + filename +
                   np.str(monno) + monname)
        if LINUX_WIN == 'w':
            plt.savefig(fileout + '.png')
        else:
            plt.savefig(fileout + '.pdf')

        plt.close()

def write_info(filename, lons, upper, lower, maximum):
    """
    # write information to a file

    """
    f1 = open(filename, 'w+')
    f1.write('longitude, upper bound latitude, lower bound latitude, ' +
                 'max_precip latitude \n')

    for mon, month in enumerate(MONTHNAMES):
        f1.write('Data for ' + month +  '\n')
        for i, lon in enumerate(lons):
            f1.write(np.str(np.around(lon, 2)) + ',' +
                     np.str(np.around(upper[mon, i], 2)) + ',' +
                     np.str(np.around(lower[mon, i], 2)) + ',' +
                     np.str(np.around(maximum[mon, i], 2)) + '\n')
    f1.close()


def main():
    print('in main')


    cubelist_AOI_pi = iris.cube.CubeList([])
    lons_pi_allmods = []
    upper_bound_lat_pi_allmods = []
    lower_bound_lat_pi_allmods = []
    most_intense_lat_pi_allmods = []

    cubelist_AOI_plio = iris.cube.CubeList([])
    lons_plio_allmods = []
    upper_bound_lat_plio_allmods = []
    lower_bound_lat_plio_allmods = []
    most_intense_lat_plio_allmods = []

    Regioninfo = GetITCZ(REGION)
    lonmin, lonmax, latmin, latmax = Regioninfo.get_region_AOI()

    for modno, modname in enumerate(MODELNAMES):

        (cube_AOI_pi, lons_pi, upper_bound_lat_pi, lower_bound_lat_pi,
         most_intense_precip_lat_pi) = Regioninfo.find_ITCZ('E280',
                                                            modname, modno)

        cubelist_AOI_pi.append(cube_AOI_pi)
        lons_pi_allmods.append(lons_pi)
        upper_bound_lat_pi_allmods.append(upper_bound_lat_pi)
        lower_bound_lat_pi_allmods.append(lower_bound_lat_pi)
        most_intense_lat_pi_allmods.append(most_intense_precip_lat_pi)

        (cube_AOI_plio, lons_plio, upper_bound_lat_plio, lower_bound_lat_plio,
         most_intense_precip_lat_plio) = Regioninfo.find_ITCZ('Eoi400',
                                                              modname, modno)

        cubelist_AOI_plio.append(cube_AOI_plio)
        lons_plio_allmods.append(lons_pi)
        upper_bound_lat_plio_allmods.append(upper_bound_lat_plio)
        lower_bound_lat_plio_allmods.append(lower_bound_lat_plio)
        most_intense_lat_plio_allmods.append(most_intense_precip_lat_plio)


        write_info(FILESTART + 'ITCZ/' + REGION + '/data_' + modname + '_E280.txt',
                   lons_pi, upper_bound_lat_pi, lower_bound_lat_pi,
                   most_intense_precip_lat_pi)

        write_info(FILESTART + 'ITCZ/' + REGION + '/data_' + modname + '_Eoi400.txt',
                   lons_plio, upper_bound_lat_plio, lower_bound_lat_plio,
                   most_intense_precip_lat_plio)


    # plot a map with the itcz on for each month and the multimodel mean

    plot_map_mmm(cubelist_AOI_plio, lons_pi, 'Pliocene')

    plot_map_mmm(cubelist_AOI_pi, lons_pi, 'Preindustrial')

    plot_map(cubelist_AOI_pi, lons_pi, upper_bound_lat_pi_allmods,
             lower_bound_lat_pi_allmods,
             most_intense_lat_pi_allmods,
             'Preindustrial')

    plot_map(cubelist_AOI_plio, lons_plio, upper_bound_lat_plio_allmods,
             lower_bound_lat_plio_allmods,
             most_intense_lat_plio_allmods,
             'Pliocene')




    # plot the maximum intensity latitude for each longitude for pliocene
    # and preindustrial for each month

    plot_lats(lons_pi_allmods, most_intense_lat_pi_allmods,
              most_intense_lat_plio_allmods,
              'Latitude of maximum precipitation',
              lonmin, lonmax, latmin, latmax,
              'Lat_max_precip_')

    # plot upper and lower bounds
    plot_lats(lons_pi_allmods, upper_bound_lat_pi_allmods,
              upper_bound_lat_plio_allmods,
              'Northern boundary of ITCZ', lonmin, lonmax, latmin, latmax,
              'Northern_bound_ITCZ_')

    plot_lats(lons_pi_allmods, lower_bound_lat_pi_allmods,
              lower_bound_lat_plio_allmods,
              'Sourthern boundary of ITCZ', lonmin, lonmax, latmin, latmax,
              'Southern_bound_ITCZ_')











########     END OF MAIN ###########

MONTHNAMES = ['January', 'February', 'March', 'April', 'May', 'June',
              'July', 'August', 'September', 'October', 'November', 'December']
#MONTHNAMES = ['ja','fb','mr','ar','my','jn','jl']
NMONTHS = len(MONTHNAMES)
LINUX_WIN = 'w'
REGION = 'Globe'

THRESHOLD = np.zeros(12, dtype=float)
THRESHOLD[0:4] = 4.
THRESHOLD[4:12] = 6.

MODELNAMES = ['CESM2', 'IPSLCM6A', 'COSMOS',
              'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
              'MIROC4m', 'IPSLCM5A2', 'HadCM3',
              'GISS2.1G', 'CCSM4',
              'CCSM4-Utr', 'CCSM4-UoT',
              'NorESM-L', 'MRI2.3', 'NorESM1-F']
#MODELNAMES = ['HadCM3', 'CESM2']

if LINUX_WIN == 'w':
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
else:
    FILESTART = '/nfs/hera1/earjcti/'

EOI400LSM = FILESTART+'regridded/PlioMIP2_Boundary_conds/Plio_enh/Plio_enh/Plio_enh_LSM_v1.0.nc'
E280LSM = FILESTART+'regridded/PlioMIP2_Boundary_conds/Modern_std/Modern_std/Modern_std_LSM_v1.0.nc'




main()
