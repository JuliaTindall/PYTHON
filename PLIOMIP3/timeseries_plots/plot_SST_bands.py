#python program
#
# this program will plot the salinity SH from files created from SH_salinity.py
import matplotlib.pyplot as plt
import numpy as np
import sys

exptnames = ['xpsic','xpsid','xpsij','xpsig','xpsie']
#exptnames = ['xpsig']

# Function to read data from a file, skipping the first line
def read_data(filename,band_required):

    band_index={'30S-30N':1,'30N-60N':2,'60N-90N':3,'30N-90N':4,
                '60S-30S':5,'90S-60S':6,'90S-30S':7}
    ix=band_index.get(band_required)
  
    x_vals = []
    y_vals = []
    with open(filename, 'r') as file:
        next(file)  # Skip the title line
        for line in file:
            try:
                values = line.strip().split(',')
                x_vals.append(float(values[0]))
                y_vals.append(float(values[ix]))
   
            except ValueError:
                print(f"Skipping invalid line in {filename}: {line.strip()}")
    return (x_vals, y_vals)
   

# Function to compute running mean**
def running_mean(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


def process_band(band,exptnames):
    """
    this will plot everything for one of the latitude bands
    """
    plt.figure(figsize=(12, 9))

    for expt in exptnames:
        filestart = ('/home/earjcti/um/' + expt +
                     '/timeseries/SST_bands_' + expt)
   
        if expt == 'xpsij':
            filename = filestart + '_1991_2999.tex'
        else:
            filename = filestart + '_12_2999.tex'
  
        (x, y) = read_data(filename,band)
 
        # Plot the data
        plt.subplot(211)
        plt.plot(x, y, label=expt + ': ' + period.get(expt) )
        plt.title('SST ' + band)

    plt.xlabel('year')
    plt.ylabel('degC')
    plt.xlim(0,1000)
    plt.legend()
    plt.grid(True)
  
    #do a 30 year running mean
    for expt in exptnames:
        filestart = ('/home/earjcti/um/' + expt +
                     '/timeseries/SST_bands_' + expt)
   
        if expt == 'xpsij':
            filename = filestart + '_1991_2999.tex'
        else:
            filename = filestart + '_12_2999.tex'
  
        x, y = read_data(filename,band)
        # calculate running mean
        window=30
        y_mean=running_mean(y,window)
        x_mean=x[window-1:]

  
        # Plot the data
        plt.subplot(212)
        plt.plot(x_mean, y_mean, label=expt + ': ' + period.get(expt) )
        plt.title('SST 60S-90S (30 year running mean)')
        plt.title('SST ' + band + '(30 year running mean)')

    plt.xlabel('year')
    plt.ylabel('degC')
    plt.legend()
    plt.xlim(1000,2000)
    plt.grid(True)

 
    # Save the plot to a file
    plt.savefig('SST_'+band+'.png')
    print("Plot saved to 'combined_plot.png'.")



##################################################################
period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400',
          'xpsig':'EP','xpsic':'PI'}

bands_required = ['30S-30N','30N-60N','60N-90N','30N-90N',
                  '60S-30S','90S-60S','90S-30S']

for band in bands_required:
    process_band(band,exptnames)

