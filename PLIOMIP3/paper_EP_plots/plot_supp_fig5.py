#python program
#
# this program will plot sea ice in volume for supplementary figure 5
#
# 
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import numpy as np
from matplotlib.gridspec import GridSpec
import sys

#exptnames = ['xpsic','xpsid','xpsig','xpsie']
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


##################################################################
period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP$_{400}$','xpsig':'EP',
          'xpsic':'PI','xqbwg':'EP','xqbwj':'LP490'}

fig = plt.figure(figsize=[10.5,6])
#spec = GridSpec(ncols=2, nrows=2,
#                         width_ratios=[1,1], wspace=0.1,
#                         hspace=0.5, height_ratios=[1, 1])


# plot sea ice area as timeseries

ax=plt.subplot(111)

for expt in exptnames:
    if expt == 'xpsij':
        filename = '../timeseries_plots/seaice/seaice_area_' + expt + '_SH_1991_2999.tex'
    else:
        filename = '../timeseries_plots/seaice/seaice_area_' + expt + '_SH_2_100.tex'
  
    x, y = read_data(filename)
    if expt == 'xpsie':
        plt.plot(x,y,label=period.get(expt),linestyle='--')
    else:
        plt.plot(x,y,label=period.get(expt))
             
plt.title('SH Sea ice area',fontsize=14)
ax.tick_params(axis='both',which='major',labelsize=14)
plt.xlabel('year',fontsize=14)
plt.ylabel('km$^2$',fontsize=14)
plt.legend()
plt.xlim(0,50)
plt.grid(True)

ax.grid(True)

ax.legend(fontsize=14)

#plt.show()

    


#plt.tight_layout() 
#plt.show()
#sys.exit()
plt.savefig('supp_figure5.png')
plt.close()
sys.exit(0)

