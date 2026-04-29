#python program
#
# this program will plot the sea ice area from files created from seaicearea.py
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import sys

exptnames = ['xpsid','xpsig','xpsie']

# Function to read data from a file, skipping the first line
def read_data(filename):
    x_vals = []
    y_vals = []
    print(filename)
    with open(filename, 'r') as file:
        next(file)  # Skip the title line
        for line in file:
            try:
                x, y = map(float, line.strip().split(','))
                x_vals.append(x)
                y_vals.append(y)
            except ValueError:
                print(f"Skipping invalid line in {filename}: {line.strip()}")
    return x_vals, y_vals

# Function to compute running mean**
def running_mean(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


# Function to read density from a file, skipping the first line
def read_dens(filename):
    print(filename)
    x_vals = []
    y5_vals = []
    y25_vals = []
    y47_vals = []
    y203_vals = []
    y447_vals = []
    y666_vals = []
    y995_vals = []
    y1500_vals = []
   
    with open(filename, 'r') as file:
        next(file)  # Skip the title line
        for line in file:
            try:
                (x, y5, y25,y47,y203,
                 y447,y666,y995,y1500) = map(float, line.strip().split(','))
                x_vals.append(x)
                y5_vals.append(y5)
                y25_vals.append(y25)
                y47_vals.append(y47)
                y203_vals.append(y203)
                y447_vals.append(y447)
                y666_vals.append(y666)
                y995_vals.append(y995)
                y1500_vals.append(y1500)
     
            except ValueError:
                print(f"Skipping invalid line in {filename}: {line.strip()}")


    ydiff=np.array(y666_vals) - np.array(y203_vals)
    
    
                
    return (x_vals, ydiff)


##################################################################
# set up figure and subplot
fig = plt.figure(figsize=[10.5,8])

spec = gridspec.GridSpec(ncols=2, nrows=3,
                         width_ratios=[2, 1], wspace=0.5,
                         hspace=0.5, height_ratios=[1, 1, 1])


period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400',
          'xpsig':'EP','xpsic':'PI'}

#plot data from all the seaice
ax0=fig.add_subplot(spec[0])
for expt in exptnames:
    if expt == 'xpsij':
        filename = 'seaice/seaice_area_' + expt + '_SH_1991_2999.tex'
    else:
        filename = 'seaice/seaice_area_' + expt + '_SH_2_2999.tex'
  
    x, y = read_data(filename)
    window=30
    y_mean=running_mean(y,window)
    x_mean=x[window-1:]
   
    # Plot the data
    if expt == 'xpsie':
        ax0.plot(x_mean, y_mean /1.0E6, label=period.get(expt),linestyle='--',
                 linewidth=3)
    else:
        ax0.plot(x_mean, y_mean /1.0E6, label=period.get(expt) )
    #plt.title('SH Sea ice area (30 year running mean)')
ax0.set_title('a) Southern Hemisphere Sea Ice',fontsize=15)
plt.ylabel('million km$^2$',fontsize=15)
plt.yticks(fontsize=15)
plt.legend(fontsize=20)
plt.xlim(0,3000)
plt.grid(True)
plt.axvline(x=1803,color='red',linestyle='--')

# Hide bottom spine and all x ticks/labels
#ax.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax0.set_xlabel(None)


#plot data from all the salinity
ax1=fig.add_subplot(spec[2])
for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/salinity_' + expt
    if expt == 'xpsij':
        filename = filestart + '_60.0S-90S_1991_2999.tex'
    else:
        filename = filestart + '_60.0S-90S_2_2999.tex'
  
  
    x, y = read_data(filename)
    window=30
    y_mean=running_mean(y,window)
    x_mean=x[window-1:]
   
    # Plot the data
    if expt == 'xpsie':
        ax1.plot(x_mean, y_mean, label=period.get(expt),linestyle='--',
                 linewidth=3)
    else:
        ax1.plot(x_mean, y_mean, label=period.get(expt) )
ax1.set_title('b) Sea Surface Salinity (60\N{degree sign}S - 90\N{degree sign}S)',fontsize=15)
plt.ylabel('psu',fontsize=18)
plt.yticks(fontsize=15)
plt.xlim(0,3000)
plt.grid(True)
plt.axvline(x=1803,color='red',linestyle='--')

# Hide bottom spine and all x ticks/labels
#ax.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax1.set_xlabel(None)



# do density
ax2=fig.add_subplot(spec[4])

for expt in exptnames:
    filestart = ('/home/earjcti/um/' + expt + '/timeseries/' + expt +
                 '_Pacific_density_levs')
    if expt == 'xpsij':
        filename = filestart + '1991_2999_-65.txt'
    else:
        filename = filestart + '12_2999_-65.txt'
  
    x, y = read_dens(filename)
    window=30
    y_mean=running_mean(y,window)
    x_mean=x[window-1:]

    # Plot the data
    if expt == 'xpsie':
        ax2.plot(x_mean, y_mean, label=period.get(expt),
                 linewidth=3,linestyle='--')
    else:
        ax2.plot(x_mean, y_mean, label=period.get(expt) )
    ax2.set_title('c) density anomaly (intermediate - near surface)',fontsize=15)

plt.xlabel('year',fontsize=15)
plt.ylabel('kg/m$^3$',fontsize=15)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.xlim(0,3000)
#plt.legend(fontsize=14)
plt.axvline(x=1803,color='red',linestyle='--')
plt.grid(True)


plt.show()
sys.exit(0)



#MARGINS = dict(left=0.14, right=0.98, bottom=0.16, top=0.96)  # tweak to suit labels
#plt.subplots_adjust(**MARGINS)


# Save the plot to a file
#plt.savefig('seaice/sea-ice_30yr_2_2999.png')
