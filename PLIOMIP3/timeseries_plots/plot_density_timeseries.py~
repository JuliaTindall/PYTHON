#python program
#
# this program will plot the AABW.  This is defined as the minimum GMOC south of
# 60S
import matplotlib.pyplot as plt
import numpy as np
import sys

#exptnames = ['xpsid','xpsic','xpsig','xpsie']
exptnames = ['xpsid','xpsig','xpsie']

# Function to read data from a file, skipping the first line
def read_data(filename):
    print(filename)
    x_vals = []
    y_vals = []
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
period = {'xpsid':'LP','xpsij':'LP490','xpsic':'PI',
          'xpsie':'EP400','xpsig':'EP'}

#plot data from all the files
for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/GMOC_' + expt
    if expt == 'xpsij':
        filename = filestart + '#1991_2999_AABW.tex'
    else:
        filename = filestart + '#12_2999_AABW.tex'
  
    x, y = read_data(filename)

    # Plot the data
    plt.subplot(211)
    plt.plot(x, y, label=expt + ': ' + period.get(expt) )
    plt.title('GMOC min 60S-90S')

plt.xlabel('year')
plt.ylabel('Sv')
plt.legend()
plt.grid(True)


#do a 30 year running mean
for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/GMOC_' + expt
    if expt == 'xpsij':
        filename = filestart + '#1991_2999_AABW.tex'
    else:
        filename = filestart + '#12_2999_AABW.tex'
  
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
    plt.title('GMOC min 60S-90S (30 year running mean)')

plt.xlabel('year',fontsize=14)
plt.ylabel('Sv',fontsize=14)
plt.legend(fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(True)
#plt.xlim(0,1000)
plt.tight_layout()

# Save the plot to a file
plt.savefig('AABW.png')
print("Plot saved to 'combined_plot.png'.")
plt.close()


#now do a 30 year running mean all by itself for the paper
plt.figure(figsize=(12, 4))

for expt in exptnames:
    filestart = '/home/earjcti/um/' + expt + '/timeseries/GMOC_' + expt
    if expt == 'xpsij':
        filename = filestart + '#1991_2999_AABW.tex'
    else:
        filename = filestart + '#12_2999_AABW.tex'
  
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
    plt.title('GMOC min 60S-90S (30 year running mean)')

plt.xlabel('year',fontsize=16)
plt.ylabel('Sv',fontsize=16)
#plt.legend(fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.grid(True)
#plt.xlim(0,1000)
#plt.tight_layout()
MARGINS = dict(left=0.14, right=0.98, bottom=0.16, top=0.96)  # tweak to suit labels
plt.subplots_adjust(**MARGINS)

# Save the plot to a file
plt.savefig('AABW_30yr.png')
print("Plot saved to 'combined_plot.png'.")
plt.close()
