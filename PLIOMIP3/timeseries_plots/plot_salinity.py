#python program
#
# this program will plot the salinity SH from files created from SH_salinity.py
import matplotlib.pyplot as plt
import numpy as np
import sys

#exptnames = ['xpsid','xpsic','xpsig','xpsie']
exptnames = ['xpsid','xpsig','xpsie']

# Function to read data from a file, skipping the first line
def read_data(filename):
    x_vals = []
    y_vals = []
    with open(filename, 'r') as file:
        next(file)  # Skip the title line  
        for line in file:
            #print(line) # this should be the title
            #sys.exit(0)
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
    filestart = '/home/earjcti/um/' + expt + '/timeseries/salinity_' + expt
   
    if expt == 'xpsij':
        filename = filestart + '_60.0S-90S_1991_2999.tex'
    else:
        filename = filestart + '_60.0S-90S_12_2999.tex'
  
    x, y = read_data(filename)

    # Plot the data
    plt.subplot(211)
    plt.plot(x, y, label=period.get(expt) )
    plt.title('SSS 60S-90S')

plt.xlabel('year')
plt.ylabel('psu')
plt.legend()
plt.grid(True)


#do a 30 year running mean
for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/salinity_' + expt
    if expt == 'xpsij':
        filename = filestart + '_60.0S-90S_1991_2999.tex'
    else:
        filename = filestart + '_60.0S-90S_12_2999.tex'
  
    x, y = read_data(filename)
    window=30
    y_mean=running_mean(y,window)
    x_mean=x[window-1:]

    # Plot the data
    plt.subplot(212)
    if expt == 'xpsie':
        plt.plot(x_mean, y_mean, label=period.get(expt),
                 linewidth=3,linestyle='--')
    else:
        plt.plot(x_mean, y_mean, label=period.get(expt) )
    plt.title('SSS 60S-90S (30 year running mean)')

plt.xlabel('year',fontsize=14)
plt.ylabel('psu',fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
#plt.xlim(0,1000)
plt.legend(fontsize=14)
plt.grid(True)
plt.tight_layout()

# Save the plot to a file
plt.savefig('salinity.png')
print("Plot saved to 'combined_plot.png'.")
plt.close()


#do the 30 year running mean on its own for the paper
plt.figure(figsize=(12, 4))

for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/salinity_' + expt
    if expt == 'xpsij':
        filename = filestart + '_60.0S-90S_1991_2999.tex'
    else:
        filename = filestart + '_60.0S-90S_12_2999.tex'
  
    x, y = read_data(filename)
    window=30
    y_mean=running_mean(y,window)
    x_mean=x[window-1:]

    # Plot the data
    if expt == 'xpsie':
        plt.plot(x_mean, y_mean, label=period.get(expt),
                 linewidth=3,linestyle='--')
    else:
        plt.plot(x_mean, y_mean, label=period.get(expt) )
    plt.title('SSS 60S-90S (30 year running mean)')

#plt.xlabel('year',fontsize=16)
plt.ylabel('psu',fontsize=16)
#plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
#plt.xlim(0,1000)
#plt.legend(fontsize=14)
plt.grid(True)

# Hide bottom spine and all x ticks/labels
#ax.spines['bottom'].set_visible(False)
ax = plt.gca()
ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)

# Remove x label if present
ax.set_xlabel(None)

MARGINS = dict(left=0.14, right=0.98, bottom=0.16, top=0.96)  # tweak to suit labels
plt.subplots_adjust(**MARGINS)

# Save the plot to a file
plt.savefig('salinity_30yr.png')
print("Plot saved to 'combined_plot.png'.")
