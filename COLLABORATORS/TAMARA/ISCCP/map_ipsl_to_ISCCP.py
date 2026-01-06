#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29.03.2022 by Julia

this will attempt to map the IPSL data onto the ISCCP weather states
"""
import numpy as np
import iris
import iris.plot as iplt
from iris.cube import CubeList as CubeList
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy
import cartopy.crs as ccrs
import sys

def plot_ISCCP_cloud_regiemes(ISCCP_cube):
    """
    plots the ISCCP cloud regiemes for reference
    """
    fig = plt.figure(figsize=[12,12])

    print(ISCCP_cube)
    plotno=0
    for i in range(0,10):
        plotno=plotno+1
        CR_cube = ISCCP_cube[:, :, i]
        CR_cube.coord('cloud top pressure dimension').rename('ctp')
        CR_cube.coord('cloud optical thickness dimension').rename('tau')
        fig.add_subplot(4,3,plotno)
        cs = plt.pcolormesh(CR_cube.data)
        cbar = plt.colorbar(cs)
        plt.yticks(ticks=np.arange(7),labels = CR_cube.coord('ctp').points)
        plt.xticks(ticks=np.arange(6),labels = CR_cube.coord('tau').points,rotation=45)
        plt.title('WS ' + str(i+1))
       
    plt.tight_layout()
    plt.savefig('ISCCP_weather_states.png')
   
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

def find_euclidian_distances(in_hist, ISCCP_cube):

    """
    in_hist is shape 7,6
    ISCCP_cube shape is 7, 6, 10
    finds 10 ISCCP climate regiemes and finds the euclidian distance between
    the input histogram and each modis climate regieme
    """
   
    sum_squares_clear = np.sum(np.square(in_hist.data))
   
    sum_sq_CR = np.zeros(10)
    for i in range(0,10):
        CR_cube = ISCCP_cube[:, :, i]
        anomaly = in_hist.data - CR_cube.data
        sum_sq_CR[i] = np.sum(np.square(anomaly))
       
    return sum_squares_clear, sum_sq_CR
 
def get_ind_cloud_regieme(hist, ISCCP_CR_cube):
    """
    we have a histogram (hist) from ipsl
    we will calculate the euclidean distance from this histogram to 
    each of the histograms in the ISCCP_CR.
    we will figure out which is the minimum euclidean distance and
    pass this back as the cloud regieme.  
    """
    R_hist = reformat_hist(hist) # reformat histogram for
                                 # comparing to ISCCP
   
    # gets the euclidian distances between CR and each of the 11 ISCCP
    # climate regiemes.  (Will return 11 numbers)
    (euclidian_distances_clear,
     euclidian_distance_CR)= find_euclidian_distances(R_hist, ISCCP_CR_cube)


    
    #cs = plt.pcolormesh(R_hist.data)
    #cbar = plt.colorbar(cs)
    #plt.yticks(ticks=np.arange(7),labels = R_hist.coord('pressure2').points)
    #plt.xticks(ticks=np.arange(6),labels = R_hist.coord('tau').points,rotation=45)

    #print('distances',euclidian_distances_clear)
    #for i, dist in enumerate(euclidian_distance_CR):
    #    print(i+1,dist)
    # add the cloud regieme with the lowest euclidian distance
    if euclidian_distances_clear == 0:
        CR_ass = 10 # clear sky
    else:
        CR_ass = np.argmin(euclidian_distance_CR)

    print('assigned cloud regieme is',CR_ass+1)

    #plt.show()
    return CR_ass, R_hist

def get_mean_histogram(histogram_cubelist, CR_no, dummy_cube):
    """
    we are passed a list of all the histograms that are in this cubelist
    we want to find the mean and pass it back
    """
    # first check that the histogram is not empty
    if len(histogram_cubelist) <= 0:
        mean_cube = dummy_cube.copy()
    elif len(histogram_cubelist) == 1:
        mean_cube = histogram_cubelist[0]
    else:
        cubes_2 = CubeList([])
        for i,cube in enumerate(histogram_cubelist):
            cube.coord('time').points = i
            cubes_2.append(cube)
        iris.util.equalise_attributes(cubes_2)
        iris.util.unify_time_units(cubes_2)
        cube = cubes_2.merge_cube()
        print(cube)
        mean_cube = cube.collapsed('time',iris.analysis.MEAN)

    print(CR_no)
    mean_cube.long_name = 'mean histogram for CR ' + str(CR_no)

    return mean_cube
        
def get_all_cloud_regiemes(ipsl_cube, ISCCP_CR_cube):
    """
    get the cloud regiemes.  Put into a cube showing all the ipsl histograms in each cloud regieme
    """


    ipsl_cube.coord('atmosphere_optical_thickness_due_to_cloud').bounds = None
    ipsl_cube.coord('pressure2').points = ipsl_cube.coord('pressure2').points / 100.

    iris.util.promote_aux_coord_to_dim_coord(ipsl_cube,'atmosphere_optical_thickness_due_to_cloud')
    lats = ipsl_cube.coord('latitude').points
    lons = ipsl_cube.coord('longitude').points
    for time in range(ndays):
        cloud_regieme_tally = np.zeros(12) # to store the number in each CR
        CR1_cubes = CubeList([])   
        CR2_cubes = CubeList([])
        CR3_cubes = CubeList([])   
        CR4_cubes = CubeList([])
        CR5_cubes = CubeList([])   
        CR6_cubes = CubeList([])
        CR7_cubes = CubeList([])   
        CR8_cubes = CubeList([])
        CR9_cubes = CubeList([])   
        CR10_cubes = CubeList([])
        CR11_cubes = CubeList([])   
        CR_clear_cubes = CubeList([])

        CR_map_array = np.zeros((len(lats),len(lons)))
        for lat in range(0,len(lats)):
      #  for lat in range(120,121):
            print(lats[lat])
            for lon in range(0,len(lons)):
                hist = ipsl_cube[time, :, :, lat,lon]
                hist.attributes = None
                hist.cell_methods = None
                hist.remove_coord('latitude')
                hist.remove_coord('longitude')

                cloud_amount = np.sum(hist.data)
                if np.isfinite(cloud_amount):
                   cloud_regieme, newhist = get_ind_cloud_regieme(hist,ISCCP_CR_cube)
                   cloud_regieme_tally[cloud_regieme] = cloud_regieme_tally[cloud_regieme] + 1
                   CR_map_array[lat,lon] = cloud_regieme + 1
                else:
                   cloud_regieme=-99999       
                                
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
                elif cloud_regieme == 10:  # clear
                    CR_clear_cubes.append(hist)

        dummy_hist = newhist.copy(data = newhist.data * 0.0)
        mean_CR1 = get_mean_histogram(CR1_cubes,1, dummy_hist)
        mean_CR2 = get_mean_histogram(CR2_cubes,2, dummy_hist)
        mean_CR3 = get_mean_histogram(CR3_cubes,3, dummy_hist)
        mean_CR4 = get_mean_histogram(CR4_cubes,4, dummy_hist)
        mean_CR5 = get_mean_histogram(CR5_cubes,5, dummy_hist)
        mean_CR6 = get_mean_histogram(CR6_cubes,6, dummy_hist)
        mean_CR7 = get_mean_histogram(CR7_cubes,7, dummy_hist)
        mean_CR8 = get_mean_histogram(CR8_cubes,8, dummy_hist)
        mean_CR9 = get_mean_histogram(CR9_cubes,9, dummy_hist)
        mean_CR10 = get_mean_histogram(CR10_cubes,10, dummy_hist)
        mean_CR11 = get_mean_histogram(CR11_cubes,11, dummy_hist)
    
        CR_map_array_cube =  ipsl_cube[0, 0, 0, :,:].copy(data=CR_map_array)

        cubelist = CubeList([mean_CR1,mean_CR2, mean_CR3, mean_CR4,
                             mean_CR5,mean_CR6, mean_CR7, mean_CR8,
                             mean_CR9,mean_CR10, mean_CR11,
                             CR_map_array_cube])

        filename = ('/home/users/jctindall/cloud_regiemes/ISCCP/'+ preind_plio_ind  + '_' + np.str(daystart + time) + '.nc')


        iris.save(cubelist,filename, netcdf_format = "NETCDF3_CLASSIC")

                                        
    # check how many assigned to each cloud regieme
    for i,CR in enumerate(cloud_regieme_tally):
        print(str(i+1), CR, CR *100. / np.sum(cloud_regieme_tally))


def get_CR_details(CRday_cubes):
    """
    This subroutine has details about cloud regiemes for each day.
    CRdaycubes contains
    a) a map showing which cloud regieme is dominant for each gridbox
    b) a picture of what each cloud regieme looks like for this day 
        (averaged over the globe)
    Firstly find how many gridboxes are for each cloud regieme
    Next extract each cloud regieme and put in a single cube
    """

    for cube in CRday_cubes:
        print(cube.long_name)
        if cube.long_name == 'mean histogram for CR 1':
            CR1_cube = cube
        if cube.long_name == 'mean histogram for CR 2':
            CR2_cube = cube
        if cube.long_name == 'mean histogram for CR 3':
            CR3_cube = cube
        if cube.long_name == 'mean histogram for CR 4':
            CR4_cube = cube
        if cube.long_name == 'mean histogram for CR 5':
            CR5_cube = cube
        if cube.long_name == 'mean histogram for CR 6':
            CR6_cube = cube
        if cube.long_name == 'mean histogram for CR 7':
            CR7_cube = cube
        if cube.long_name == 'mean histogram for CR 8':
            CR8_cube = cube
        if cube.long_name == 'mean histogram for CR 9':
            CR9_cube = cube
        if cube.long_name == 'mean histogram for CR 10':
            CR10_cube = cube
        if cube.long_name == 'mean histogram for CR 11':
            CR11_cube = cube
        if cube.long_name == 'ISCCP Cloud Area Fraction':
            CR_map_cube = cube

   # qplt.contourf(CR_map_cube)
    CR_map_data = CR_map_cube.data
    cr_no = np.zeros(13)
    for CR in range(0,12):
        cr_no[CR] = np.sum(np.where(CR_map_data == CR+1, 1.0, 0.0))

    CR_cubelist = CubeList([CR1_cube, CR2_cube, CR3_cube, CR4_cube,
                                              CR5_cube, CR6_cube, CR7_cube, CR8_cube,
                                              CR9_cube, CR10_cube, CR11_cube])
 
    return (cr_no, CR_cubelist)

def plot_avg_cloud_regieme():
    """
    plots the average of the histograms in each cloud regieme
    the histograms are from the cosp2 test data
    they are all in a file
    """

    filestart = '/home/users/jctindall/cloud_regiemes/ISCCP/' + preind_plio_ind
    fig = plt.figure(figsize=[12,12])
    plotno=0

    total_in_CR = np.zeros(13)
    firsttimein = 'y'
    sum_of_CR_cubelist = CubeList([])

    for day in range (daystart,daystart + ndays):
        filename = filestart + '_' + str(day) + '.nc'
        try:
            CRday_cubes = iris.load(filename)
            print('found ' + filename)
            
        except:
            print('not found ' + filename)
            continue

        number_in_CR, day_CR_cubelist = get_CR_details(CRday_cubes)
        print(day_CR_cubelist)

        # set up for averaging
        total_in_CR = total_in_CR + number_in_CR

        for CR in range(0,11):
            cube = day_CR_cubelist[CR]
            cubetot = cube.copy(data=cube.data * number_in_CR[CR])
            if firsttimein == 'y':
                sum_of_CR_cubelist.append(cubetot)
                print('appending CR', len(sum_of_CR_cubelist))
            else:
                if number_in_CR[CR] > 0:
                    sum_of_CR_cubelist[CR].data = sum_of_CR_cubelist[CR].data + cubetot.data
        firsttimein = 'n'
    
    # find average of each cloud regieme
    CR_cubelist = CubeList([])
    print(len(sum_of_CR_cubelist))
    for CR, cube in enumerate(sum_of_CR_cubelist):
        newcube = cube.copy(data = cube.data / total_in_CR[CR])
        CR_cubelist.append(newcube)

    print(len(CR_cubelist),len(sum_of_CR_cubelist))
       
    for CR,CR_cube in enumerate(CR_cubelist):
        percentage = np.round(total_in_CR[CR] * 100.  / np.sum(total_in_CR),2)
        print(CR+1, total_in_CR[CR], percentage)
        # plot it
        fig.add_subplot(4,3,CR+1)
        cs = plt.pcolormesh(CR_cube.data)
        cbar = plt.colorbar(cs)
        plt.yticks(ticks=np.arange(7),labels =CR_cube.coord('pressure2').points)
        plt.xticks(ticks=np.arange(6),labels = CR_cube.coord('tau').points,rotation=45)
        plt.title('WS :' + str(CR+1) + ' number:' + str(np.round(total_in_CR[CR])) + ' (' + str(percentage) + '%)')

    plt.tight_layout()
    plt.savefig('figures/isccp_WS_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.eps')
    plt.savefig('figures/isccp_WS_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.pdf')

def plot_avg_cloud_regieme_Tselioudis():
    """
    plots the average of the histograms in each cloud regieme
    for comparison with tseulioudis we group (3 and 6) and (9 and 10)
    """
    filestart = '/home/users/jctindall/cloud_regiemes/ISCCP/' + preind_plio_ind
    fig = plt.figure(figsize=[12,12])
    plotno=0

    total_in_CR = np.zeros(9)
    sum_of_CR_cubelist = CubeList([])
        
    for day in range (daystart,daystart + ndays):
        filename = filestart + '_' + str(day) + '.nc'
        try:
            CRday_cubes = iris.load(filename)
            print('found ' + filename)
            
        except:
            print('not found ' + filename)
            continue

      
        # this gets the number in each cloud regieme and the
        # average histogram
        number_in_CR, day_CR_cubelist = get_CR_details(CRday_cubes)
        if day == daystart:
            # set up empty cubelist for storing 9 new weather states
            zerocube = day_CR_cubelist[0] *0.0
            sum_of_CR_cubelist = CubeList([zerocube.copy(),zerocube.copy()
                        ,zerocube.copy(),zerocube.copy(),zerocube.copy()
                        ,zerocube.copy(),zerocube.copy(),zerocube.copy()
                                           ,zerocube.copy()])
    
        
        # convert from original weather states to tesioudis et al weather
        # states
        total_in_CR[0] = total_in_CR[0] + number_in_CR[0] # new WS1
        total_in_CR[1] = total_in_CR[1] + number_in_CR[1] # new WS2
        total_in_CR[2] = (total_in_CR[2] + number_in_CR[2] +
                          number_in_CR[5]) # new WS3
        total_in_CR[3] = total_in_CR[3] + number_in_CR[3] # new WS4
        total_in_CR[4] = total_in_CR[4] + number_in_CR[4] # new WS5
        total_in_CR[5] = total_in_CR[5] + number_in_CR[6] # new WS6
        total_in_CR[6] = total_in_CR[6] + number_in_CR[7] # new WS7
        total_in_CR[7] = (total_in_CR[7] + number_in_CR[8]
                          + number_in_CR[9]) # new WS8
        total_in_CR[8] = total_in_CR[8] + number_in_CR[10] # new WS9 - clear

        # get the cloud regiemes in the new weather states
        # new WS1 and WS2 and ws4 and ws 5
        for i in (0,1,3,4):
            cube = day_CR_cubelist[i]
            cubetot = cube.copy(data=cube.data * number_in_CR[i])  
            sum_of_CR_cubelist[i].data = sum_of_CR_cubelist[i].data + cubetot.data
        # new WS3 (old ws3+ws6)
        cube2 = day_CR_cubelist[2]
        cubetot2 = cube2.copy(data=cube2.data * number_in_CR[2])  
        cube5 = day_CR_cubelist[5]
        cubetot5 = cube5.copy(data=cube5.data * number_in_CR[5])  
        sum_of_CR_cubelist[2].data = (sum_of_CR_cubelist[2].data +
                                          cubetot2.data + cubetot5.data)

        # new WS6 and new ws7
        for i in (6,7):
            cube = day_CR_cubelist[i]
            cubetot = cube.copy(data=cube.data * number_in_CR[i])  
            sum_of_CR_cubelist[i-1].data = sum_of_CR_cubelist[i-1].data + cubetot.data
      
        # new WS8 (old ws9+ws10)
        cube8 = day_CR_cubelist[8]
        cubetot8 = cube8.copy(data=cube8.data * number_in_CR[8])  
        cube9 = day_CR_cubelist[9]
        cubetot9 = cube9.copy(data=cube9.data * number_in_CR[9])  
        sum_of_CR_cubelist[7].data = (sum_of_CR_cubelist[7].data +
                                          cubetot8.data + cubetot9.data)

        # new WS9 clear sky (old WS11)
        cube10 = day_CR_cubelist[10]
        cubetot10 = cube10.copy(data=cube10.data * number_in_CR[10])  
        sum_of_CR_cubelist[8].data = (sum_of_CR_cubelist[8].data +
                                          cubetot10.data)

        
    # find average of each cloud regieme
    CR_cubelist = CubeList([])
    print(len(sum_of_CR_cubelist))
    for CR, cube in enumerate(sum_of_CR_cubelist):
        newcube = cube.copy(data = cube.data / total_in_CR[CR])
        CR_cubelist.append(newcube)

    print(len(CR_cubelist),len(sum_of_CR_cubelist))
       
    for CR,CR_cube in enumerate(CR_cubelist):
        percentage = np.round(total_in_CR[CR] * 100.  / np.sum(total_in_CR),2)
        print(CR+1, total_in_CR[CR], percentage)
        # plot it
        fig.add_subplot(3,3,CR+1)
        cs = plt.pcolormesh(CR_cube.data)
        cbar = plt.colorbar(cs)
        plt.yticks(ticks=np.arange(7),labels =CR_cube.coord('pressure2').points)
        plt.xticks(ticks=np.arange(6),labels = CR_cube.coord('tau').points,rotation=45)
        plt.title('WS :' + str(CR+1) + ' number:' + str(np.round(total_in_CR[CR])) + ' (' + str(percentage) + '%)')

    plt.tight_layout()
    plt.savefig('figures/isccp_Tse_WS_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.eps')
    plt.savefig('figures/isccp_Tse_WS_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.png')
    plt.savefig('figures/isccp_Tse_WS_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.pdf')


def plot_map_RFO_cloud_regieme():
    """ 
    plots a map for each cloud regieme showing where it occurs
    """

    filestart = '/home/users/jctindall/cloud_regiemes/ISCCP/' + preind_plio_ind
    fig = plt.figure(figsize=[12,12])
    plotno=0

    firsttimein = 'y'
    CR_cubelist = CubeList([])

    daysused=0.
    for day in range (daystart,daystart + ndays):
        filename = filestart + '_' + str(day) + '.nc'
        try:
            CRday_cubes = iris.load_cube(filename,'ISCCP Cloud Area Fraction')
            daysused=daysused+1.
            
        except:
            print('not found ' + filename)
            continue

        # we have the data.  Now create a cube for each cloud regieme
        CRday_data = CRday_cubes.data
        for CR in range(1,13):
            crdata = np.where(CRday_data==CR, 1.0, 0.0)
            crcube = CRday_cubes.copy(data=crdata)
            crcube.long_name = 'cloud_regieme_' + str(CR) + '_cloud_area_percentage'

            if firsttimein == 'y':
                CR_cubelist.append(crcube)
            else:
                CR_cubelist[CR-1].data = CR_cubelist[CR-1].data + crcube.data
        firsttimein = 'n'

    # CR_cubelist contains 12 cubes each showing a map of where the cloud
    # regiemes are

    print('daysused',daysused)
   # jucmap= ListedColormap(["white", "mistyrose", "plum", "lightpink", "plum", "mediumpurple","cornflowerblue","skyblue","aqua","cyan","turquoise","springgreen","greenyellow","yellow","orange","red"])
    jucmap= ListedColormap(["white","lavenderblush", "pink", "lightpink", "plum","mediumslateblue","turquoise","lightskyblue","aquamarine","springgreen","darkgreen","green", "yellow","gold","orange","red"])
    V = [0,1,2,3,4,6,8,10,13,16,20,25,30,40,50,70,90]
    V = np.arange(0,93.5,5.5)
    fig, axs = plt.subplots(nrows=4,ncols=3,
                           #subplot_kw={'projection': ccrs.Mollweide(central_longitude=180.)},
                           subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180.)},
                           figsize=(11.0,11.0))
    axs=axs.flatten() # axs is 2d array flatten to 1d array.

    for CR,cube in enumerate(CR_cubelist):
        cube.data = cube.data * 100 / daysused

        # plot a map of each cloud regieme

        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        print('plotting CR',CR)
        cs = axs[CR].contourf(lons,lats,cube.data,cmap=jucmap,levels=V,
                              extend='max',transform=ccrs.PlateCarree(),
                             norm=colors.BoundaryNorm(boundaries=V,ncolors=17))
        #cbar = axs[CR].colorbar(cs,orientation='horizontal',ticks=V)
        axs[CR].coastlines()
        axs[CR].set_title('WS :' + str(CR+1))

    fig.subplots_adjust(bottom=0.2)
    cbar_ax = fig.add_axes([0.10, 0.15, 0.80, 0.03])
    fig.colorbar(cs, cax=cbar_ax,orientation='horizontal',ticks=V)
        
   # plt.tight_layout()
    plt.savefig('figures/isccp_WS_RFO_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.eps')
    plt.savefig('figures/isccp_WS_RFO_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.pdf')
    iris.save(CR_cubelist,'figures/isccp_WS_RFO_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.nc',netcdf_format="NETCDF3_CLASSIC")


def plot_map_RFO_cloud_regieme_Tselioudis():
    """ 
    plots a map for each cloud regieme showing where it occurs
    for comparison with Tselioudis we merge CR 3-6 and CR9-10
    """

    filestart = '/home/users/jctindall/cloud_regiemes/ISCCP/' + preind_plio_ind
    fig = plt.figure(figsize=[12,10])
    plotno=0

    CR_cubelist = CubeList([])

    for day in range (daystart,daystart + ndays):
   # for day in range(daystart, daystart+3):
        print('day',day)
        filename = filestart + '_' + str(day) + '.nc'
        try:
            CRday_cubes = iris.load_cube(filename,'ISCCP Cloud Area Fraction')
            print('found ' + filename)
            
        except:
            print('not found ' + filename)
            continue

        if day == daystart:
            # set up empty cubelist for storing 9 new weather states
            zerocube = CRday_cubes *0.0
            CR_cubelist = CubeList([zerocube.copy(),zerocube.copy()
                        ,zerocube.copy(),zerocube.copy(),zerocube.copy()
                        ,zerocube.copy(),zerocube.copy(),zerocube.copy()
                                           ,zerocube.copy()])
    
       
        # we have the data.  Now create a cube for each cloud regieme
        CRday_data = CRday_cubes.data
        for CR in range(1,13):
            crdata = np.where(CRday_data==CR, 1.0, 0.0)
            crcube = CRday_cubes.copy(data=crdata)
            crcube.long_name = 'cloud_regieme_' + str(CR) + '_cloud_area_percentage'
            if CR in (1,2,4,5):
                CR_cubelist[CR-1].data = CR_cubelist[CR-1].data + crdata
            if CR == 3 or CR==6:  # merge 3 and 6
                CR_cubelist[2].data = CR_cubelist[2].data + crdata
            if CR in (7,8): # shift one weather state down
                CR_cubelist[CR-2].data = CR_cubelist[CR-2].data + crdata
            if CR == 9 or CR==10:  # merge 9 and 10
                CR_cubelist[7].data = CR_cubelist[7].data + crdata
            if CR == 11: # clear sky now CR9
                CR_cubelist[8].data = CR_cubelist[8].data + crdata
         
   # jucmap= ListedColormap(["white", "mistyrose", "plum", "lightpink", "plum", "mediumpurple","cornflowerblue","skyblue","aqua","cyan","turquoise","springgreen","greenyellow","yellow","orange","red"])
    jucmap= ListedColormap(["white", "lavenderblush", "pink", "lightpink", "plum","mediumslateblue","turquoise","lightskyblue","aquamarine","springgreen","darkgreen","green", "yellow","gold","orange","red"])
    V = [0,1,2,3,4,6,8,10,13,16,20,25,30,40,50,70,90]
    V = np.arange(0,93.5,5.5)
    fig, axs = plt.subplots(nrows=3,ncols=3,
                           #subplot_kw={'projection': ccrs.Mollweide(central_longitude=180.)},
                           subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180.)},
                           figsize=(11.0,11.0))
    axs=axs.flatten() # axs is 2d array flatten to 1d array.

    for CR,cube in enumerate(CR_cubelist):
        cube.data = cube.data * 100 / ndays
        cube.long_name = 'cloud_regieme_' + str(CR) + '_cloud_area_percentage'
        
        # plot a map of each cloud regieme

        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        print('plotting CR',CR)
        cs = axs[CR].contourf(lons,lats,cube.data,cmap=jucmap,levels=V,
                              extend='max',transform=ccrs.PlateCarree(),
                             norm=colors.BoundaryNorm(boundaries=V,ncolors=17))
        #cbar = axs[CR].colorbar(cs,orientation='horizontal',ticks=V)
        axs[CR].coastlines()
        axs[CR].set_title('WS :' + str(CR+1))

    fig.subplots_adjust(bottom=0.2)
    cbar_ax = fig.add_axes([0.10, 0.15, 0.80, 0.03])
    fig.colorbar(cs, cax=cbar_ax,orientation='horizontal',ticks=V)
        
   # plt.tight_layout()
    plt.savefig('figures/isccp_Tse_WS_RFO_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.eps')
    plt.savefig('figures/isccp_Tse_WS_RFO_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.png')
    plt.savefig('figures/isccp_Tse_WS_RFO_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.pdf')
    iris.save(CR_cubelist,'figures/isccp_Tse_WS_RFO_' + preind_plio_ind +'_' + str(daystart) + '_' + str(daystart + ndays) + '.nc',netcdf_format="NETCDF3_CLASSIC")



def plot_map_RFO_anomaly():
    """ 
    plots a map for each cloud regieme anomaly.  (To run this section
    you will previously need to have run the plot_map_RFO_cloud_regieme
    section for both the Pliocene and the preindustrial.  This will produce
    the required files
    """

    pifile = 'figures/isccp_Tse_WS_RFO_PREIND_' + str(daystart) + '_' + str(daystart + ndays) + '.nc'
    pliofile = 'figures/isccp_Tse_WS_RFO_PLIO_' + str(daystart) + '_' + str(daystart + ndays) + '.nc'

    picubelist = iris.load(pifile)
    pliocubelist = iris.load(pliofile)
    print(pliocubelist)
    # check all percentages add to 100
    for i,cube in enumerate(picubelist):
        if i==0:
            pi_data = picubelist[i].data
        else:
            pi_data = pi_data + picubelist[i].data
    pi_sumcube = picubelist[0].copy(data=pi_data)
    for i,cube in enumerate(pliocubelist):
        if i==0:
            plio_data = pliocubelist[i].data
        else:
            plio_data = plio_data + pliocubelist[i].data
    plio_sumcube = pliocubelist[0].copy(data=plio_data)

    #plt.subplot(2,1,1)
    #qplt.contourf(pi_sumcube)
    #plt.subplot(2,1,2)
    #qplt.contourf(plio_sumcube)
    #plt.show()
    
    anomcubelist = CubeList([])

    for CR in range(1,10):

        namereq = 'cloud_regieme_' + str(CR-1) + '_cloud_area_percentage'
        
        for cube in picubelist:
            if cube.long_name ==namereq:
                picube = cube
        for cube in pliocubelist:
            print(cube.var_name,namereq)
            if cube.long_name ==namereq:
                pliocube = cube
        anomcube = pliocube - picube

        print(CR, anomcube.data[40,40],pliocube.data[40,40],picube.data[40,40],
              picube.coord('latitude').points[40],
              picube.coord('longitude').points[40])
        anomcubelist.append(anomcube)
    fig = plt.figure(figsize=[12,10])
    plotno=0


    # CR_cubelist contains 12 cubes each showing a map of where the cloud
    # regiemes are

   
    jucmap= ListedColormap(["white", "lavenderblush", "pink", "lightpink", "plum", "mediumpurple","mediumslateblue","turquoise","lightskyblue","aquamarine","springgreen","lime","palegreen", "greenyellow","gold","orange","red"])
    #V = [0,1,2,3,4,6,8,10,13,16,20,25,30,40,50,70,90]
    V = np.arange(-10,11,1)

    fig, axs = plt.subplots(nrows=3,ncols=3,
                           #subplot_kw={'projection': ccrs.Mollweide(central_longitude=180.)},
                           subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180.)},
                           figsize=(11.0,11.0))
    axs=axs.flatten() # axs is 2d array flatten to 1d array.

    for CR,cube in enumerate(anomcubelist):
        # plot a map of each cloud regieme

        lats = cube.coord('latitude').points
        lons = cube.coord('longitude').points
        print('plotting CR',CR)
        cs = axs[CR].contourf(lons,lats,cube.data,cmap='RdBu_r',levels=V,
                              extend='both',transform=ccrs.PlateCarree())
#                             norm=colors.BoundaryNorm(boundaries=V,ncolors=16))
        #cbar = axs[CR].colorbar(cs,orientation='horizontal',ticks=V)
        axs[CR].coastlines()
        axs[CR].set_title('WS :' + str(CR+1))

    fig.subplots_adjust(bottom=0.2)
    cbar_ax = fig.add_axes([0.10, 0.15, 0.80, 0.03])
    fig.colorbar(cs, cax=cbar_ax,orientation='horizontal',ticks=V)
        
   # plt.tight_layout()
    plt.savefig('figures/isccp_WS_RFO_anom_' + str(daystart) + '_' + str(daystart + ndays) + '.eps')
    plt.savefig('figures/isccp_WS_RFO_anom_' + str(daystart) + '_' + str(daystart + ndays) + '.pdf')


###################################################
def create_netcdf_1year():
    """
    find the time index where the year starts and ends
    create the file for all the weather states
    """

    # find all the dates and the indexes we need for this year
    cube=iris.load_cube(filename)
    time = cube.coord('time')
    print(time)
    dates = time.units.num2date(time.points)
    dayreq=[]
    datesreq=[]
    for i,date in enumerate(dates):
        year = date.year
        month= date.month
        day = date.day
   
        if year == YEARREQ:
           dayreq.append(i)
           datesreq.append(date)
    # read in the weather state for those dates
    oneyearcubelist = CubeList([])
    filestart = '/home/users/jctindall/cloud_regiemes/ISCCP/' + preind_plio_ind
    for fileno in dayreq:
        print(fileno)
        file = filestart + '_' + str(fileno) + '.nc'
        cube = iris.load_cube(file,'ISCCP Cloud Area Fraction')
        ncube=iris.util.new_axis(cube, scalar_coord='time')
        ncube.coord('time').points=time.points[fileno]
        ncube.coord('time').bounds=None
        ncube.coord('time').var_name = 'time'
        oneyearcubelist.append(ncube)
      
   
    #catcube=oneyearcubelist.concatenate()
    #print(catcube)
    oneyearcube = oneyearcubelist.concatenate_cube()
    oneyearcube.long_name='ISCCP weather state'
    print(oneyearcube)

    fileout = ('/home/users/jctindall/cloud_regiemes/ISCCP/'+ preind_plio_ind  + '_year_' + str(YEARREQ) + '.nc')


    #iris.save(oneyearcube,fileout, netcdf_format = "NETCDF3_CLASSIC")
    iris.save(oneyearcube,fileout)
    print('saved cube')
       
        
########    START OF PROGRAM ########################
# modis cloud regieme has index ctp, cot, cr
ISCCP_CR_cube = iris.load_cube('ISCCP_weather_states.nc')
#plot_ISCCP_cloud_regiemes(ISCCP_CR_cube)

#isccp histogram has index ctp, cot, loc

preind_plio_ind = 'PLIO' #  optionS PREIND PLIO
ndays=365
daystart = 0
YEARREQ=2075  # plio first year - 2025, preind first year = 2075

filestart = '/home/users/jctindall/cloud_regiemes/clisccp_CFday_IPSL-CM6A-LR_'
if preind_plio_ind == 'PREIND':
    filename = filestart + 'piControl_r1i2p1f1_gr_20750101-20981231.nc'
if preind_plio_ind == 'PLIO':
    filename = filestart + 'midPliocene-eoi400_r1i1p1f1_gr_20250101-20491231.nc'

cubes = iris.load(filename)
cube = cubes[0]
ipsl_cube = cube[daystart:daystart + ndays,:,:,:,:]
cube = 0


####################################################################
# get the cloud regieme for each histogram and group the histograms
# this program will also write them out to a file

#get_all_cloud_regiemes(ipsl_cube, ISCCP_CR_cube)
#sys.exit(0)

##################################################################
# plot the average histogram from the test data in each cloud regieme.
# it uses the files written out in the previous steop

#plot_avg_cloud_regieme()
#plot_avg_cloud_regieme_Tselioudis()
#plot_map_RFO_cloud_regieme()
plot_map_RFO_cloud_regieme_Tselioudis()
plot_map_RFO_anomaly()
#sys.exit(0)

###################################################################
# creates a netcdf file containing all of the weather states for 1 year

#create_netcdf_1year()
