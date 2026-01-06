#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29.03.2022 by Julia

this is an experimental program to see what is in the MOSIS-CR***nc4 file
"""
import numpy as np
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import sys

def plot_MODIS_cloud_regiemes(MODIS_cube):
    """
    plots the MODIS cloud regiemes for reference
    """
    fig = plt.figure(figsize=[15,15])

    plotno=0
    for i in range(0,11):
        plotno=plotno+1
        CR_cube = MODIS_cube[:, :, i]
        CR_cube.coord('cloud top pressure dimension').rename('ctp')
        CR_cube.coord('cloud optical thickness dimension').rename('tau')
        fig.add_subplot(4,3,plotno)
        cs = plt.pcolormesh(CR_cube.data)
        cbar = plt.colorbar(cs)
        plt.yticks(ticks=np.arange(7),labels = CR_cube.coord('ctp').points)
        plt.xticks(ticks=np.arange(6),labels = CR_cube.coord('tau').points,rotation=45)
        plt.title('CR ' + str(i+1))
       
    plt.tight_layout()
    plt.savefig('MODIS_cloud_regiemes.eps')
   
def reformat_hist(orighist):
    """
    note the MODIS data has 6 categories of cloud optical thickness
    the cosp2 histograms have 7.  Combine the lower two catogeries so that
    they match.  
    """

    auxcoordpoints = np.copy(orighist.coord('atmosphere_optical_thickness_due_to_cloud').points)
    auxcoordpoints[0] = 0.635
    auxcoordpoints[1] = 0.635
  
    taucoord = iris.coords.AuxCoord(points = auxcoordpoints, long_name = 'tau')

    
    orighist.add_aux_coord(taucoord, data_dims = 1)
    
    newhist = orighist.aggregated_by(taucoord, iris.analysis.SUM)
    iris.util.promote_aux_coord_to_dim_coord(newhist,'tau')

    return newhist

def find_euclidian_distances(in_hist, MODIS_cube):

    """
    in_hist is shape 7,6
    MODIS_cube shape is 7, 6, 12
    finds 12 MODIS climate regiemes and finds the euclidian distance between
    the input histogram and each modis climate regieme
    """
   
    sum_squares_clear = np.sum(np.square(in_hist.data))
   
    sum_sq_CR = np.zeros(11)
    for i in range(0,11):
        CR_cube = MODIS_cube[:, :, i]
        anomaly = in_hist.data - CR_cube.data
        sum_sq_CR[i] = np.sum(np.square(anomaly))
       
    return sum_squares_clear, sum_sq_CR
 
def get_ind_cloud_regieme(cosp2hist, MODIS_CR_cube):
    """
    we have a histogram (hist) from cosp2
    we will calculate the euclidean distance from this histogram to 
    each of the histograms in the MODIS_CR.
    we will figure out which is the minimum euclidean distance and
    pass this back as the cloud regieme.  
    """
    cosp2R_hist = reformat_hist(cosp2hist) # reformat histogram for
                                           # comparing to MODIS

    # gets the euclidian distances between CR and each of the 12 MODIS
    # climate regiemes.  (Will return 12 numbers)
    (euclidian_distances_clear,
     euclidian_distance_CR)= find_euclidian_distances(cosp2R_hist, 
                                                      MODIS_CR_cube)

    #print('distances',euclidian_distances_clear)
    #for i, dist in enumerate(euclidian_distance_CR):
    #    print(i+1,dist)
    # add the cloud regieme with the lowest euclidian distance
    if euclidian_distances_clear == 0:
        CR_ass = 11 # clear sky
    else:
        CR_ass = np.argmin(euclidian_distance_CR)
   

    return CR_ass, cosp2R_hist


def get_all_cloud_regiemes(cosp2_test_cube, MODIS_CR_cube):
    """
    get the cloud regiemes.  Put into a cube showing all the cosp2 histograms in each cloud regieme
    """

    cloud_regieme_tally = np.zeros(12) # to store the number in each Cloud Regieme
    CR1_cubes = iris.cube.CubeList([])   
    CR2_cubes = iris.cube.CubeList([])
    CR3_cubes = iris.cube.CubeList([])   
    CR4_cubes = iris.cube.CubeList([])
    CR5_cubes = iris.cube.CubeList([])   
    CR6_cubes = iris.cube.CubeList([])
    CR7_cubes = iris.cube.CubeList([])   
    CR8_cubes = iris.cube.CubeList([])
    CR9_cubes = iris.cube.CubeList([])   
    CR10_cubes = iris.cube.CubeList([])
    CR11_cubes = iris.cube.CubeList([])   
    CR_clear_cubes = iris.cube.CubeList([])

    for histno in range(0,1236):
        hist = cosp2_test_cube[:, :, histno]
        cloud_regieme, newhist = get_ind_cloud_regieme(hist,MODIS_CR_cube)
  
        cloud_regieme_tally[cloud_regieme] = cloud_regieme_tally[cloud_regieme] + 1
    
        if cloud_regieme == 0:  # CR1
            CR1_cubes.append(newhist)
        elif cloud_regieme == 1:  # CR2
            CR2_cubes.append(newhist)
        elif cloud_regieme == 2:  # CR3
            CR3_cubes.append(newhist)
        elif cloud_regieme == 3:  # CR4
            CR4_cubes.append(newhist)
        elif cloud_regieme == 4:  # CR5
            CR5_cubes.append(newhist)
        elif cloud_regieme == 5:  # CR6
            CR6_cubes.append(newhist)
        elif cloud_regieme == 6:  # CR7
            CR7_cubes.append(newhist)
        elif cloud_regieme == 7:  # CR8
            CR8_cubes.append(newhist)
        elif cloud_regieme == 8:  # CR9
            CR9_cubes.append(newhist)
        elif cloud_regieme == 9:  # CR10
            CR10_cubes.append(newhist)
        elif cloud_regieme == 10:  # CR11
            CR11_cubes.append(newhist)
        elif cloud_regieme == 11:  # clear
            CR_clear_cubes.append(hist)
        else:
            print('missing cloud regieme',cloud_regieme)
            sys.exit(0)
   
    # save the cubelists for each cloud regiems
    if len(CR1_cubes) > 0:
        iris.save(CR1_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR1.nc')
    if len(CR2_cubes) > 0:
        iris.save(CR2_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR2.nc')
    if len(CR3_cubes) > 0:
        iris.save(CR3_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR3.nc')
    if len(CR4_cubes) > 0:
        iris.save(CR4_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR4.nc')
    if len(CR5_cubes) > 0:
        iris.save(CR5_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR5.nc')
    if len(CR6_cubes) > 0:
        iris.save(CR6_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR6.nc')
    if len(CR7_cubes) > 0:
        iris.save(CR7_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR7.nc')
    if len(CR8_cubes) > 0:
        iris.save(CR8_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR8.nc')
    if len(CR9_cubes) > 0:
        iris.save(CR9_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR9.nc')
    if len(CR10_cubes) > 0:
        iris.save(CR10_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR10.nc')
    if len(CR11_cubes) > 0:
        iris.save(CR11_cubes,MODIS_ISCCP_IND + '_map_cosp2test_CR11.nc')
  
  
    # check how many assigned to each cloud regieme
    for i,CR in enumerate(cloud_regieme_tally):
        print(str(i+1), CR)
   
def plot_avg_cloud_regieme():
    """
    plots the average of the histograms in each cloud regieme
    the histograms are from the cosp2 test data
    they are all in a file
    """
    
    fig = plt.figure(figsize=[15,15])
    plotno=0
 
    for CR in range(1,12):
        plotno=plotno+1
        
        try:
            cubes = iris.load(MODIS_ISCCP_IND + '_map_cosp2test_CR' + str(CR) + '.nc')
        except:
            continue

        # merge cubes for averaging
        cubes_2 = iris.cube.CubeList([])
        iris.util.equalise_attributes(cubes)
        for i,cube in enumerate(cubes):
            cube.attributes = None
            cube2 = iris.util.new_axis(cube,scalar_coord='loc')
            cube2.coord('loc').points = i
            cube2.coord('loc').var_name = None
            cube2.var_name = None
            cube2.cell_methods = None
            
            print(cube2.metadata)
            
            cubes_2.append(cube2)
        
        cube = cubes_2.concatenate_cube()
        mean_cube = cube.collapsed('loc',iris.analysis.MEAN)

        # plot it
        fig.add_subplot(4,3,plotno)
        cs = plt.pcolormesh(mean_cube.data)
        cbar = plt.colorbar(cs)
        plt.yticks(ticks=np.arange(7),labels = mean_cube.coord('air_pressure').points)
        plt.xticks(ticks=np.arange(6),labels = mean_cube.coord('tau').points,rotation=45)
        plt.title('CR :' + str(CR) + ' number:' + np.str(i+1))
       
    plt.tight_layout()
    plt.savefig(MODIS_ISCCP_IND + '_cosp2_avg_cloud_regiemes.eps')
        


########    START OF PROGRAM ########################
# modis cloud regieme has index ctp, cot, cr
MODIS_CR_cube = iris.load_cube('MODIS_cloud_regiemes.nc')
#plot_MODIS_cloud_regiemes(MODIS_CR_cube)

#isccp histogram has index ctp, cot, loc

MODIS_ISCCP_IND = 'Modis'  # options Modic, ISCCP

cubes = iris.load('/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/COLLABORATORS/TAMARA/CFMIP/COSPv2.0_juliatest/driver/data/outputs/UKMO/cosp2_output_um_julia_default.nc','cloud_area_fraction_in_atmosphere_layer')  
for cube in cubes:
    if MODIS_ISCCP_IND == 'Modis' and cube.var_name == 'clmodis':
        cosp2_test_cube = cube
    if MODIS_ISCCP_IND == 'ISCCP' and cube.var_name == 'clisccp':
        cosp2_test_cube = cube

####################################################################
# get the cloud regieme for each histogram and group the histograms
# this program will also write them out to a file

#get_all_cloud_regiemes(cosp2_test_cube, MODIS_CR_cube)


##################################################################
# plot the average histogram from the test data in each cloud regieme.
# it uses the files written out in the previous steop

plot_avg_cloud_regieme()
