#
#python program
#
# this program will plot the spinup from the experiments
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import iris
#import iris.quickplot as qplt
import iris.plot as iplt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys

exptnames = ['xpsid','xpsig','xpsie','xpsij']
exptnames = ['d','g','e','j']

# Function to read data from a file, skipping the first line
def read_data(filename,deliminator):
    x_vals = []
    y_vals = []
    print(filename)
    with open(filename, 'r') as file:
        next(file)  # Skip the title line
        for line in file:
            try:
                x, y = map(float, line.strip().split(deliminator))
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
fig = plt.figure(figsize=[10.5,12])
spec = fig.add_gridspec(ncols=1, nrows=5, height_ratios=[1,1,1,1,0.4])

period = {'d':'LP','j':'LP$_{490}$','e':'EP$_{400}$',
          'g':'EP','c':'PI'}

#plot TOA radiative inbalance
ax0=fig.add_subplot(spec[0])
for expt in exptnames:
    filename = ('/home/earjcti/um/xpsi' + expt + '/timeseries/xpsi'
                + expt + '_TOAinbal.tex')
  
    x, y = read_data(filename,' ')
    window=100
    y_mean=running_mean(y,window)
    x_mean=x[window//2 : -(window//2)+1]
   
    # Plot the data
    #if expt == 'xpsie':
    #    ax0.plot(x_mean, y_mean /1.0E6, label=period.get(expt),linestyle='--',
    #             linewidth=3)
    #else:
    ax0.plot(x_mean, y_mean, label=period.get(expt),linestyle='--',linewidth=3 )
    #plt.title('SH Sea ice area (30 year running mean)')
ax0.set_title('a) TOA radiative inbalance (100 year running mean)',fontsize=15)
plt.ylabel('W / m$^2$',fontsize=15)
plt.yticks(fontsize=15)
#plt.legend(fontsize=20)
handles,labels = ax0.get_legend_handles_labels()
plt.xlim(0,4000)
plt.grid(True)


# Hide bottom spine and all x ticks/labels
#ax0.spines['bottom'].set_visible(False)
ax0.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax0.set_xlabel(None)

#legend
axleg=fig.add_subplot(spec[4])
plt.legend(handles,labels,ncol=4,loc='lower center',fontsize=15)
axleg.spines[:].set_visible(False)
axleg.tick_params(axis='both', which='both', bottom=False, top=False,
                  right=False, left=False,labelbottom=False,
                  labeltop=False,labelleft=False,labelright=False)

##############################
#plot data from temperature top levels
ax1=fig.add_subplot(spec[1])
for expt in exptnames:
    filename = ('/home/earjcti/um/xpsi' + expt + '/timeseries/tdrift_xpsi'
                + expt + '_lev1_10.txt')
    
    x, y = read_data(filename,',')
    y_mean=running_mean(y,window)
    x_mean=x[window//2 : -(window//2)+1]

   
    # Plot the data
    ax1.plot(x_mean, y_mean, label=period.get(expt),linestyle='--',
                 linewidth=3)

    ax1.set_title('b) mean ocean temperature (0 - 360 m; 100 year running mean)',fontsize=15)
plt.ylabel('$^\circ$C',fontsize=15)
ax1.set_yticks(np.arange(15.4,17.1,0.3))
plt.yticks(fontsize=15)
plt.xlim(0,4000)
plt.grid(True)
#plt.axvline(x=1803,color='red',linestyle='--')

# Hide bottom spine and all x ticks/labels
#ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax1.set_xlabel(None)

###############################
# plot data from temperature all levels
ax1=fig.add_subplot(spec[2])
for expt in exptnames:
    filename = ('/home/earjcti/um/xpsi' + expt + '/timeseries/tdrift_xpsi'
                + expt + '_alllevs.txt')
    
    x, y = read_data(filename,',')
    y_mean=running_mean(y,window)
    x_mean=x[window//2 : -(window//2)+1]

   
    # Plot the data
    ax1.plot(x_mean, y_mean, label=period.get(expt),linestyle='--',
                 linewidth=3)

    ax1.set_title('c) global mean ocean temperature (100 year running mean)',fontsize=15)
plt.ylabel('$^\circ$C',fontsize=15)
#ax1.set_yticks([33.2,33.6,34.0])
plt.yticks(fontsize=15)
plt.xlim(0,4000)
plt.grid(True)
#plt.axvline(x=1803,color='red',linestyle='--')

# Hide bottom spine and all x ticks/labels
#ax1.spines['bottom'].set_visible(False)
ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax1.set_xlabel(None)

###########################
#plot data from salinity
# note all levels is all levels. 1-10 is levels 1-10 which is 0-360 m;  1-5 is 0-55m

ax1=fig.add_subplot(spec[3])
for expt in exptnames:
    filename = ('/home/earjcti/um/xpsi' + expt + '/timeseries/drift_sal_xpsi'
                + expt + '_lev1_10.txt')
    
    x, y = read_data(filename,',')
    y_mean=running_mean(y,window)
    x_mean=x[window//2 : -(window//2)+1]

   
    # Plot the data
    ax1.plot(x_mean, y_mean,label=period.get(expt),linestyle='--',
                 linewidth=3)

    ax1.set_title('d) mean ocean salinity (0 - 360 m; 100 year running mean)',fontsize=15)
plt.ylabel('psu',fontsize=15)
plt.xlabel('model year',fontsize=15)
#ax1.set_yticks([33.2,33.6,34.0])
plt.yticks(fontsize=15)
plt.xticks(fontsize=15)
plt.xlim(0,4000)
plt.grid(True)
#plt.axvline(x=1803,color='red',linestyle='--')

# Hide bottom spine and all x ticks/labels
#ax.spines['bottom'].set_visible(False)
#ax1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
#ax1.set_xlabel(None)

plt.tight_layout()
plt.savefig('spinup.png')


# Save the plot to a file
#plt.savefig('seaice/sea-ice_30yr_2_2999.png')
