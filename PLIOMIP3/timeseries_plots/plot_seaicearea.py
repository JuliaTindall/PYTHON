#python program
#
# this program will plot the sea ice area from files created from seaicearea.py
import matplotlib.pyplot as plt
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


##################################################################
plt.figure(figsize=(12, 9))
period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400',
          'xpsig':'EP','xpsic':'PI'}

#plot data from all the files
for expt in exptnames:
    if expt == 'xpsij':
        filename = 'seaice_area_' + expt + '_SH_1991_2999.tex'
    else:
        filename = 'seaice_area_' + expt + '_SH_12_2999.tex'
  
    x, y = read_data(filename)

    # Plot the data
    plt.subplot(211)
    plt.plot(x, y, label=period.get(expt) )
    plt.title('SH Sea ice area')

plt.xlabel('year')
plt.ylabel('km2')
plt.legend()
plt.grid(True)


#do a 30 year running mean
for expt in exptnames:
    if expt == 'xpsij':
        filename = 'seaice_area_' + expt + '_SH_1991_2999.tex'
    else:
        filename = 'seaice_area_' + expt + '_SH_12_2999.tex'
  
    x, y = read_data(filename)
    window=30
    y_mean=running_mean(y,window)
    x_mean=x[window-1:]

    # Plot the data
    plt.subplot(212) 
    if expt == 'xpsie':
        plt.plot(x_mean, y_mean /1.0E6, label=period.get(expt),linestyle='--',
                 linewidth=3)
    else:
        plt.plot(x_mean, y_mean /1.0E6, label=period.get(expt) )
    plt.title('SH Sea ice area (30 year running mean)')

plt.xlabel('year',fontsize=14)
plt.ylabel('million km$^2$',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(fontsize=14)
#plt.xlim(1000,2000)
plt.grid(True)

plt.tight_layout()


# Save the plot to a file
plt.savefig('sea-ice.png')
plt.close()

# now do a 30 year running mean on its own for the paper
#do a 30 year running mean
plt.figure(figsize=(12, 4))

for expt in exptnames:
    if expt == 'xpsij':
        filename = 'seaice_area_' + expt + '_SH_1991_2999.tex'
    else:
        filename = 'seaice_area_' + expt + '_SH_12_2999.tex'
  
    x, y = read_data(filename)
    window=30
    y_mean=running_mean(y,window)
    x_mean=x[window-1:]

    # Plot the data
    if expt == 'xpsie':
        plt.plot(x_mean, y_mean /1.0E6, label=period.get(expt),linestyle='--',
                 linewidth=3)
    else:
        plt.plot(x_mean, y_mean /1.0E6, label=period.get(expt) )
    plt.title('SH Sea ice area (30 year running mean)')

#plt.xlabel('year',fontsize=16)
plt.ylabel('million km$^2$',fontsize=16)
#plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=20)
#plt.xlim(1000,2000)
plt.grid(True)
# remove x axis
ax = plt.gca()

# Hide bottom spine and all x ticks/labels
#ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

# Remove x label if present
ax.set_xlabel(None)

MARGINS = dict(left=0.14, right=0.98, bottom=0.16, top=0.96)  # tweak to suit labels
plt.subplots_adjust(**MARGINS)


# Save the plot to a file
plt.savefig('sea-ice_30yr.png')
