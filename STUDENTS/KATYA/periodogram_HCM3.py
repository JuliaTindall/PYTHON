#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 08/09/2020

@author: Katya

This is katyas periodogram code that she wants help with

"""
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from matplotlib.dates import YearLocator, DateFormatter
import matplotlib.dates as mdates
import matplotlib.patches as mpatches
#import seaborn as sns

#import mtspec

from scipy.signal import butter,filtfilt
import statistics
from scipy import signal
import iris
import iris.quickplot as qplt 
import iris.plot as iplt
#iris.FUTURE.netcdf_promote=True
import iris.coord_categorisation


# Upload Eoi400
sst_constraint = iris.Constraint(cube_func=lambda cube: cube.var_name == 'temp')
sst_constraint1 = iris.Constraint(cube_func=lambda cube: cube.var_name == 'tos')
sst_constraint2 = iris.Constraint(cube_func=lambda cube: cube.var_name == 'THO')
sst_constraint3 = iris.Constraint(cube_func=lambda cube: cube.var_name == 'TS')
sst_constraint4 = iris.Constraint(cube_func=lambda cube: cube.var_name == 'sst')
temporary = iris.load_cube('/nfs/hera1/ee22kvcs/EOI400.SST.2001-2501_timeseries_HadCM3_no_ann_cycle.nc', sst_constraint)
HadCM3_Eoi400_SST = temporary[0:6000,:,:]

# Upload E280
temporary = iris.load_cube('/nfs/hera1/ee22kvcs/E280.SST.2411-2911_timeseries_HadCM3_no_ann_cycle.nc', sst_constraint)
HadCM3_E280_SST = temporary[0:6000,:,:]

# Delimit the region for the cube
min_lat = -5.0
max_lat = 5.0
min_lat12 = -10.0
max_lat12 = 0.0
min_lon34 = 190.0
max_lon34 = 240.0
min_lon12 = 270.0
max_lon12 = 280.0

# Functions to identify latitudes and longitudes we want to subset to
def nino_lat(input):
    return min_lat  <= input <= max_lat 
def nino_lat12(input):
    return min_lat12  <= input <= max_lat12 
def nino_lon34(input):
    return min_lon34  <= input <= max_lon34
def nino_lon12(input):
    return min_lon12  <= input <= max_lon12 

nino_34con = iris.Constraint(latitude = nino_lat, longitude = nino_lon34)
nino_12con = iris.Constraint(latitude = nino_lat12, longitude = nino_lon12)

# Subset the cube to the location of interest
MODELNAMES_MP = [HadCM3_Eoi400_SST]
MODELNAMES_PI = [HadCM3_E280_SST]

nino34_mp_subsets = []
nino12_mp_subsets = []
for i in MODELNAMES_MP:
    nino34_mp_sub = nino_34con.extract(i)
    nino12_mp_sub = nino_12con.extract(i)
    nino34_mp_subsets.append(nino34_mp_sub)
    nino12_mp_subsets.append(nino12_mp_sub)

nino34_pi_subsets = []
nino12_pi_subsets = []
for j in MODELNAMES_PI:
    nino34_pi_sub = nino_34con.extract(j)
    nino12_pi_sub = nino_12con.extract(j)
    nino34_pi_subsets.append(nino34_pi_sub)
    nino12_pi_subsets.append(nino12_pi_sub)
    
## 3-month running mean
window_size = 3
nino34_mp_running_mean = []
nino34_pi_running_mean = []
nino12_mp_running_mean = []
nino12_pi_running_mean = []
# Calculate the running mean
for e in range(len(MODELNAMES_MP)):
    running_mean_nino34_mp = nino34_mp_subsets[e].rolling_window('time', iris.analysis.MEAN, window=window_size)
    running_mean_nino34_pi = nino34_pi_subsets[e].rolling_window('time', iris.analysis.MEAN, window=window_size)
    running_mean_nino12_mp = nino12_mp_subsets[e].rolling_window('time', iris.analysis.MEAN, window=window_size)
    running_mean_nino12_pi = nino12_pi_subsets[e].rolling_window('time', iris.analysis.MEAN, window=window_size)
    nino34_mp_running_mean.append(running_mean_nino34_mp)
    nino34_pi_running_mean.append(running_mean_nino34_pi)
    nino12_mp_running_mean.append(running_mean_nino12_mp)
    nino12_pi_running_mean.append(running_mean_nino12_pi)
    
# Extract the running mean values
print(nino34_mp_running_mean[0].shape)
print(nino34_pi_running_mean[0].shape)

# Spatial average
mean_mp34 = []
mean_mp12 = []
mean_pi34 = []
mean_pi12 = []
for d in range(len(MODELNAMES_MP)):
    equator_mean_mp34 = nino34_mp_running_mean[d].collapsed(['latitude', 'longitude'],iris.analysis.MEAN)
    equator_mean_mp12 = nino12_mp_running_mean[d].collapsed(['latitude', 'longitude'],iris.analysis.MEAN)
    equator_mean_pi34 = nino34_pi_running_mean[d].collapsed(['latitude', 'longitude'],iris.analysis.MEAN)
    equator_mean_pi12 = nino12_pi_running_mean[d].collapsed(['latitude', 'longitude'],iris.analysis.MEAN)
    # Store
    mean_mp34.append(equator_mean_mp34)
    mean_mp12.append(equator_mean_mp12)
    mean_pi34.append(equator_mean_pi34)
    mean_pi12.append(equator_mean_pi12)

   
## Normalise and filter the data
filtered_norm_mp34 = []
filtered_norm_pi34 = []
filtered_norm_mp12 = []
filtered_norm_pi12 = []

for mod in range(1):
    mp34_stdev = mean_mp34[mod].collapsed('time', iris.analysis.STD_DEV)
    pi34_stdev = mean_pi34[mod].collapsed('time', iris.analysis.STD_DEV)
    mp12_stdev = mean_mp12[mod].collapsed('time', iris.analysis.STD_DEV)
    pi12_stdev = mean_pi12[mod].collapsed('time', iris.analysis.STD_DEV)
    mp34_mean = mean_mp34[mod].collapsed('time', iris.analysis.MEAN)
    pi34_mean = mean_pi34[mod].collapsed('time', iris.analysis.MEAN)
    mp12_mean = mean_mp12[mod].collapsed('time', iris.analysis.MEAN)
    pi12_mean = mean_pi12[mod].collapsed('time', iris.analysis.MEAN)
    
    norm_mp34 = (mean_mp34[mod].data-mp34_mean.data)/mp34_stdev.data
    norm_pi34 = (mean_pi34[mod].data-pi34_mean.data)/pi34_stdev.data
    norm_mp12 = (mean_mp12[mod].data-mp12_mean.data)/mp12_stdev.data
    norm_pi12 = (mean_pi12[mod].data-pi12_mean.data)/pi12_stdev.data
    ## Normalisation
    # Define the sampling rate
    sampling_rate = 1
    cutoff_freq = 1/18
    # Convert cutoff frequency to a normalized frequency
    nyquist_freq = 0.5 * sampling_rate
    normalized_cutoff_freq = cutoff_freq / nyquist_freq
    # Butterworth low-pass filter, order of 4
    order = 4  # changed from 4
    b_butter, a_butter = butter(order, normalized_cutoff_freq, btype='low', analog=False)
    # Apply the filter
    filtered_mp34 = filtfilt(b_butter, a_butter, norm_mp34)
    filtered_pi34 = filtfilt(b_butter, a_butter, norm_pi34)
    filtered_mp12 = filtfilt(b_butter, a_butter, norm_mp12)
    filtered_pi12 = filtfilt(b_butter, a_butter, norm_pi12)
    filtered_norm_mp34.append(filtered_mp34)
    filtered_norm_pi34.append(filtered_pi34)
    filtered_norm_mp12.append(filtered_mp12)
    filtered_norm_pi12.append(filtered_pi12)

    print('julia1')
    
def periodogram(model, model_name):
    print('julia2',model,model_name)
    # a and c are frequency
    # d and d are power
    a, b = signal.periodogram(filtered_norm_mp34[0].data)
    c, d = signal.periodogram(filtered_norm_pi34[0].data)
    # Transform the frequency to period
    a_period = 1 / (a * 12)
    c_period = 1 / (c * 12)
    # Try to normalise the 
    norm_b = b/np.std(filtered_norm_mp34[model].data)
    norm_d = d/np.std(filtered_norm_pi34[model].data)
    
    fig = plt.figure(figsize=(12,12))
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.plot(a_period, norm_b, linewidth=2, label="Eoi$^{400}$")
    ax1.plot(c_period, norm_d, linewidth=2, alpha=0.7, label="E$^{280}$")
    
    plt.xlim(0, 500)
    ax1.set_xticks(np.arange(0, 15 + 1, 1))

    plt.ylabel('Power', fontsize=19)
    plt.xlabel('Period (years)', fontsize=18)
    ax1.xaxis.set_tick_params(labelsize=15)
    ax1.yaxis.set_tick_params(labelsize=15)

    plt.title('Periodogram of El Niño 3.4 index - '+(model_name), fontsize=20)
    #plt.savefig('/nfs/see-fs-01_teaching/ee22kvcs/Task2/Periodogram_Nino34_' + model_name+ '.png',  bbox_inches = 'tight', dpi=300, format='png')
    ax1.legend(fontsize=16)
  

    ax2 = fig.add_subplot(2, 2, 2)
    ax2.plot(a, b, linewidth=2, label="Eoi$^{400}$")
    ax2.plot(c, d, linewidth=2, alpha=0.7, label="E$^{280}$")
    ax2.plot([0., 0.5],[1.0, 1.0],color='red')

    # calculate total power
    freqdiff = a[1]-a[0]
    sumpower = np.sum(b *freqdiff)
    sumpower2 = np.sum(d *freqdiff)
    print('sumpower',sumpower,sumpower2)
    #sys.exit(0)

    plt.xlim(0, 0.1)
    #ax2.set_xticks(np.arange(0,  + 1, 1))

    plt.ylabel('Power', fontsize=19)
    plt.xlabel('frequency', fontsize=18)
    ax2.xaxis.set_tick_params(labelsize=15)
    ax2.yaxis.set_tick_params(labelsize=15)

    plt.title('nonnormalised', fontsize=20)
    #plt.savefig('/nfs/see-fs-01_teaching/ee22kvcs/Task2/Periodogram_Nino34_' + model_name+ '.png',  bbox_inches = 'tight', dpi=300, format='png')
    ax1.legend(fontsize=16)
    plt.show()
    return(c_period)
    
## Run the function
MODELNAMES = ['CCSM4-Utr','CESM1.2','COSMOS','HadCM3','IPSLCM6A','MIROC4m','NorESM-L']

PLOTS = periodogram(0, 'HadCM3')
