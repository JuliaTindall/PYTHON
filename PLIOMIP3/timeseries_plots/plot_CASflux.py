#python program
#
# this program will plot the CAS flux created from flux_through_CAS.py
import matplotlib.pyplot as plt
import numpy as np
import sys

exptnames = ['xpsie','xpsig']
which_to_plot = 'heat'  # this is 'water' or 'heat'

# Function to read data from a file, skipping the first line
def read_data(filename):
    year_vals = []
    flow_vals = []
    heatflow_vals=[]
    with open(filename, 'r') as file:
        next(file)  # Skip the title line
        for line in file:
            try:
                x, y, z, a, b = map(float, line.strip().split(','))
                year_vals.append(x)
                flow_vals.append(y+z)
                heatflow_vals.append(a+b)
            except ValueError:
                print(f"Skipping invalid line in {filename}: {line.strip()}")
    return year_vals, flow_vals, heatflow_vals

# Function to compute running mean**
def running_mean(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')


##################################################################
plt.figure(figsize=(12, 9))
period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP400','xpsig':'EP'}

#plot data from all the files
for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/'

    if expt == 'xpsij':
        filename = filestart + 'CASflow_' + expt + '_1991_2999.tex'
    else:
        filename = filestart + 'CASflow_' + expt + '_12_2999.tex'
  
    year, waterflow, heatflow= read_data(filename)

    if which_to_plot == 'water':
        y=waterflow
        ylab='Sv'
    if which_to_plot == 'heat':
        y=heatflow
        ylab='PW'
    
    # Plot the data
    plt.subplot(211)
    plt.plot(year, y, label=expt + ': ' + period.get(expt,expt) )
    plt.title(which_to_plot + ' Flux throught CAS')

plt.xlabel('year')
plt.ylabel(ylab)
plt.legend()
plt.grid(True)


#do a 30 year running mean
for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/'

    if expt == 'xpsij':
        filename = filestart + 'CASflow_' + expt + '_1991_2999.tex'
    else:
        filename = filestart + 'CASflow_' + expt + '_12_2999.tex'

    year, waterflow, heatflow= read_data(filename)

    if which_to_plot == 'water':
        y=waterflow
        ylab='Sv'
    if which_to_plot == 'heat':
        y=heatflow
        ylab='PW'

    window=30
    y_mean=running_mean(y,window)
    x_mean=year[window-1:]

    # Plot the data
    plt.subplot(212)
    plt.plot(x_mean, y_mean, label=expt + ': ' + period.get(expt,expt) )
    plt.title(which_to_plot + ' flow through CAS (30 year running mean)')

plt.xlabel('year')
plt.ylabel(ylab)
plt.legend()
plt.grid(True)
#plt.show()
#sys.exit(0)

# Save the plot to a file
plt.savefig(which_to_plot + '_total_flow_through_CAS.png')
print("Plot saved to 'combined_plot.png'.")
