# -*- coding: utf-8 -*-
"""
Created on Sat Jun  6 08:43:12 2020

@author: julia
"""

#import sys
import iris
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np


def get_filenames(model):
    """
    gets the preindustrial and the pliocene filename for the given model
    """
    file_pi = (FILESTART + 'regridded/' + model +
               '/E280.TotalPrecipitation.mean_month.nc')
    file_plio = (FILESTART + 'regridded/' + model +
                 '/EOI400.TotalPrecipitation.mean_month.nc')

    return file_pi, file_plio

def get_zm_cube(file_pi, file_plio):
    """
    gets the zonal mean of the cube that is contained in the file
    """

    cube_pi = iris.load_cube(file_pi)
    cube_plio = iris.load_cube(file_plio)

    cube_zm_pi = cube_pi.collapsed(['longitude'], iris.analysis.MEAN)
    cube_zm_plio = cube_plio.collapsed(['longitude'], iris.analysis.MEAN)

    # remove auxillary coordinate year and time
    cube_pi.remove_coord('year')
    cube_plio.remove_coord('year')
    cube_pi.remove_coord('time')
    cube_plio.remove_coord('time')

    # remove attributes from month coordinate
    cube_pi.coord('month').attributes = None
    cube_plio.coord('month').attributes = None

    return cube_zm_pi, cube_zm_plio, cube_pi, cube_plio

def plot_cubes(cubelist_zm_pi, cubelist_zm_plio):
    """
    plots the preindustrial and pliocene zonal mean precipitation that
    is in the cubes
    """

    for mon, monthname in enumerate(MONTHNAMES):
        fig = plt.figure(figsize=[12.0, 8.0])
        for i, cube_pi in enumerate(cubelist_zm_pi):
            axis = plt.subplot(4, 4, i+1)
            cube_plio = cubelist_zm_plio[i]
            axis.plot(cube_pi[mon, :].data, cube_pi.coord('latitude').points,
                      color='blue')
            axis.plot(cube_plio[mon, :].data, cube_pi.coord('latitude').points,
                      color='orange')
            plt.text(5, 20, MODELNAMES[i], fontsize=10)
            plt.ylim(-30, 30)
            plt.xlim(0, 10)
            plt.axhline(y=0, linestyle='dashed', color='black')
            if i < 12: # only show xaxmis on bottom plot
                axis.axes.xaxis.set_visible(False)
            if i % 4 != 0: # only show yaxis on left plot
                axis.axes.yaxis.set_visible(False)
        fig.text(0, 0.95,
                 'Tropical precipitation (plio-orange, pi-blue) for:'
                 + monthname,
                 fontsize=25)
        fig.text(0.5, 0.05, 'mm/day', fontsize=20)
        fig.text(0.05, 0.5, 'latitude', fontsize=20, rotation=90)


        fileout = (FILESTART + 'ITCZ/Tropical_precip/absolute_values_'
                   +str(mon) + monthname + '.pdf')
        plt.savefig(fileout)
        fileout = (FILESTART + 'ITCZ/Tropical_precip/absolute_values_'
                   + str(mon) + monthname + '.png')
        plt.savefig(fileout)
        plt.close()


def plot_anomaly_cubes(cubelist_zm_pi, cubelist_zm_plio):
    """
    plots the preindustrial and pliocene zonal mean precipitation that
    is in the cubes
    """

    for mon, monthname in enumerate(MONTHNAMES):
        fig = plt.figure(figsize=[12.0, 8.0])
        for i, cube_pi in enumerate(cubelist_zm_pi):
            axis = plt.subplot(4, 4, i+1)
            cube_plio = cubelist_zm_plio[i]
            anomaly_data = cube_plio[mon, :].data - cube_pi[mon, :].data
            axis.plot(anomaly_data, cube_pi.coord('latitude').points)
            plt.text(0.5, 20, MODELNAMES[i], fontsize=10)
            plt.ylim(-30, 30)
            plt.xlim(-2, 2)
            plt.axhline(y=0, linestyle='dashed', color='black')
            plt.axvline(x=0, linestyle='dashed', color='black')
            if i < 12: # only show xaxmis on bottom plot
                axis.axes.xaxis.set_visible(False)
            if i % 4 != 0: # only show yaxis on left plot
                axis.axes.yaxis.set_visible(False)
        fig.text(0, 0.95, 'Tropical precipitation anomaly (plio - pi) for '  + monthname,
                 fontsize=25)
        fig.text(0.5, 0.05, 'mm/day', fontsize=20)
        fig.text(0.05, 0.5, 'latitude', fontsize=20, rotation=90)


        fileout = (FILESTART + 'ITCZ/Tropical_precip/anomalies_'
                   +str(mon) + monthname + '.pdf')
        plt.savefig(fileout)
        fileout = (FILESTART + 'ITCZ/Tropical_precip/anomalies_'
                   + str(mon) + monthname + '.png')
        plt.savefig(fileout)
        plt.close()


def plot_map_cubes(cubelist_pi, cubelist_plio):
    """
    plots the preindustrial and pliocene mean precipitation that
    is in the cubes for the globe
    """

    for mon, monthname in enumerate(MONTHNAMES):
        fig = plt.figure(figsize=[12.0, 6.0])
        for i, cube_pi in enumerate(cubelist_pi):
            plt.subplot(4, 4, i+1)
            cube_plio = cubelist_plio[i]
            anomaly_cube = cube_plio - cube_pi
            trop_anom_month = anomaly_cube.extract(iris.Constraint(month=mon+1,
                                                                   latitude = lambda cell: -30 < cell < 30))
            V = np.arange(-5, 6, 1)
            cs = iplt.contourf(trop_anom_month, levels=V, cmap='RdBu', extend='both')
            plt.gca().coastlines()
            plt.title(MODELNAMES[i], fontsize=10)

        # add colorbar and full plot title
        plt.subplots_adjust(left=0.0, bottom=0.2, right=0.9, top=0.9, wspace=0.1, hspace=0.0)
        cb_ax = fig.add_axes([0.10, 0.15, 0.7, 0.05])
        cbar = fig.colorbar(cs, cax=cb_ax, orientation='horizontal')
        cbar.set_label('mm/day', fontsize=15)
        fig.text(0, 0.95, 'Tropical precipitation anomaly (plio - pi) for '  + monthname,
                 fontsize=25)



        fileout = (FILESTART + 'ITCZ/Tropical_precip/global_anomalies_'
                   +str(mon) + monthname + '.pdf')
        plt.savefig(fileout)
        fileout = (FILESTART + 'ITCZ/Tropical_precip/global_anomalies_'
                   + str(mon) + monthname + '.png')
        plt.savefig(fileout)
        plt.close()




def main():
    """
    Will plot the tropical precipitation for the pliocene and the preindustrial
    from each of the models
    """
    cubelist_zm_e280 = iris.cube.CubeList([])
    cubelist_zm_eoi400 = iris.cube.CubeList([])
    cubelist_e280 = iris.cube.CubeList([])
    cubelist_eoi400 = iris.cube.CubeList([])
    for modname in MODELNAMES:
        filein_e280, filein_eoi400 = get_filenames(modname)
        (cube_zm_e280, cube_zm_eoi400,
         cube_e280, cube_eoi400) = get_zm_cube(filein_e280, filein_eoi400)
        cubelist_zm_e280.append(cube_zm_e280)
        cubelist_zm_eoi400.append(cube_zm_eoi400)
        cubelist_e280.append(cube_e280)
        cubelist_eoi400.append(cube_eoi400)

    #plot_cubes(cubelist_zm_e280, cubelist_zm_eoi400)
    #plot_anomaly_cubes(cubelist_zm_e280, cubelist_zm_eoi400)
    plot_map_cubes(cubelist_e280, cubelist_eoi400)




LINUX_WIN = 'w'
MODELNAMES = ['CESM2', 'IPSLCM6A', 'COSMOS',
              'EC-Earth3.3', 'CESM1.2', 'IPSLCM5A',
              'MIROC4m', 'IPSLCM5A2', 'HadCM3',
              'GISS2.1G', 'CCSM4',
              'CCSM4-Utr', 'CCSM4-UoT',
              'NorESM-L', 'MRI2.3', 'NorESM1-F'
             ]
#MODELNAMES = ['CESM2']
MONTHNAMES = ['January', 'February', 'March', 'April', 'May', 'June',
              'July', 'August', 'September', 'October', 'November', 'December']
if LINUX_WIN == 'w':
    FILESTART = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
else:
    FILESTART = '/nfs/hera1/earjcti/'

main()
