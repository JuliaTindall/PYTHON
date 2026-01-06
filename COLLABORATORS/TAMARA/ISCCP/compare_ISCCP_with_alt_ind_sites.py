#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 21.12.2022 by Julia

we have two ways of mapping the ipsl clouds onto ISCCP.  These are
1. map_ipsl_to_ISCCP.py  which minimises the euclidean distance between
   the histograms
2. map_ipsl_to_ISCCP_alternative_williams_webb.py  which finds the centroids
   of tcc, ctp, albedo and minimises the euclidean distance between the
   centroids

This program will read the westher states from both method at a site.
It will show a graph of the weather state by month
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

def get_weather_states():
    """
    reads in the weather states from the files
    """
    filestart = '/home/users/jctindall/cloud_regiemes/'
    fileend = preind_plio_ind + '_year_' + str(YEAR) + '.nc'
    orig_ws_file = filestart + '/ISCCP/' + fileend
   # orig_ws_file = filestart + '/ISCCP/PREIND_year_2075.nc'
    alt_ws_file = filestart + '/ISCCP_alternative/' + fileend


    orig_cube = iris.load_cube(orig_ws_file,'cloud_area_fraction_in_atmosphere_layer')
    alt_cube = iris.load_cube(alt_ws_file,'cloud_area_fraction_in_atmosphere_layer')

    lats = orig_cube.coord('latitude').points
    latix = np.abs(lats - SITELAT).argmin()
    latuse = lats[latix]
    
    lons = orig_cube.coord('longitude').points
    lonix=np.abs(lons - SITELON).argmin()
    lonuse = lons[lonix]


    # extract the data
    lonconstraint = iris.Constraint(longitude=lonuse)
    latconstraint = iris.Constraint(latitude=latuse)
    orig_at_site_cube = orig_cube.extract(lonconstraint & latconstraint)
    alt_at_site_cube = alt_cube.extract(lonconstraint & latconstraint)
    print(orig_at_site_cube.data)
    
    return orig_at_site_cube, alt_at_site_cube


def plot_weather_states(orig_ws_cube, alt_ws_cube):
    """
    plots the weather states so we can compare
    """
    plt.subplot(2,1,1)
    times = orig_ws_cube.coord('time').points
    times = times - times[0]
    plt.bar(times, orig_ws_cube.data, align='center')
    plt.subplot(2,1,2)
    plt.bar(times, alt_ws_cube.data, align='center')
 
    return times

def weather_states_statistics(times,orig_method, alt_method):
    """
    calculates some statistics for the weather states
    1. How many are 0 # weather states not found
    2. How many are 7 in the alt method - these will remain as 7
    3. How many are left after 1 and 2 above
    4. See how many of the remaining ones agree / disagree
    5. Count how many of each weather state are in the agree/disagree band
    6. Plot when we are getting our 7's, our agreements and our disagreements
    """

    # 1.2.3.
    zeromask = np.where(orig_method == 0, 1.0, 0.0)
    sevenmask = np.where(alt_method == 7, 1.0, 0.0)
    remainingmask = np.where(zeromask+sevenmask == 0.0, 1.0, 0.0)
    print(zeromask + sevenmask, np.sum(zeromask), np.sum(sevenmask), np.sum(remainingmask))

    #4  see how many of the remaining ones agree
    agreemask = np.where((remainingmask == 1.0) & (orig_method == alt_method), 1.0, 0.0)
    disagreemask = np.where((remainingmask == 1.0) & (orig_method != alt_method), 1.0, 0.0)
    print('agree ',np.sum(agreemask), 'disagree',np.sum(disagreemask))

    # 5 count how many of each weather state are in the agree/disagree band
    WS_agree = np.zeros(12)
    WS_disagree_origmethod = np.zeros(12)
    WS_disagree_altmethod = np.zeros(12)
    for ws in range(1,12):
        WS_agree[ws-1] =  np.sum(np.where((orig_method * agreemask) == ws, 1.0, 0.0))
        WS_disagree_origmethod[ws-1] =  np.sum(np.where((orig_method * disagreemask) == ws, 1.0, 0.0))
        WS_disagree_altmethod[ws-1] =  np.sum(np.where((alt_method * disagreemask) == ws, 1.0, 0.0))
        print('WS',ws,':agree',WS_agree[ws-1],'disagree: orig=',WS_disagree_origmethod[ws-1],'disagree: alt=',WS_disagree_altmethod[ws-1])
    
    # 6 plot when we are getting stuff
    times = orig_ws_cube.coord('time').points
    times = times - times[0]
    plt.subplot(3,1,1)
    plt.bar(times, alt_method * sevenmask, align='center')
    plt.title('sevens')
    plt.subplot(3,1,2)
    plt.bar(times, alt_method * agreemask, align='center')
    plt.title('agree')
    plt.subplot(3,1,3)
    plt.bar(times, alt_method * disagreemask, align='center')
    plt.title('disagree - alternative regieme')
    plt.show()
    sys.exit(0)



    

########    START OF PROGRAM ########################

preind_plio_ind = 'PLIO' #  optionS PREIND PLIO
sitename = 'BeaverPond' # 79N, 82W
SITELAT = 79.0
SITELON = 278.0
YEAR = 2025

orig_ws_cube, alt_ws_cube = get_weather_states() # get w states at location
times = orig_ws_cube.coord('time').points
times = times - times[0]
#times = plot_weather_states(orig_ws_cube, alt_ws_cube)
weather_states_statistics(times,orig_ws_cube.data, alt_ws_cube.data)

#plt.show()
sys.exit(0)
# find where orig_ws and alt_ws agree
orig_ws_data = orig_ws.data
alt_ws_data = alt_ws.data
agree_data = np.ma.where(orig_ws_data == alt_ws_data, orig_ws_data,-99999.)
agree_cube = orig_ws.copy(data=agree_data)
agree_cube.long_name = 'both_methods_agree'
number_agree = np.sum(np.where(agree_data > 0, 1.0, 0.0))
number_zero = np.sum(np.where(agree_data == 0, 1.0, 0.0))
print('number in agreement',number_agree,'total=',144*143,'zero=',number_zero)
print('percentage in disagreement',((144.*143.) - np.float64(number_agree) - np.float64(number_zero)) * 100. / (144.*143.))

# find where they agree or where alt_ws is 7
orig_7_data = np.ma.where(alt_ws_data==7, alt_ws_data,orig_ws_data)
agree_or_7_data = np.ma.where(orig_7_data == alt_ws_data, alt_ws_data,-99999.)
agree_or_7_cube = orig_ws.copy(data=agree_or_7_data)
agree_or_7_cube.long_name = 'both_methods_agree or low_cloud_fraction'
number_agree_7 = np.sum(np.where(agree_or_7_data > 0, 1.0, 0.0))
print('number in agreementor 7',number_agree_7)
print('percentage in disagreement excl 7',((144.*143.) - np.float64(number_agree_7) - np.float64(number_zero)) * 100. / (144.*143.))


cubelist = [agree_cube, agree_or_7_cube]
iris.save(cubelist,'temporary.nc',fill_value=-99999.,netcdf_format="NETCDF3_CLASSIC")
