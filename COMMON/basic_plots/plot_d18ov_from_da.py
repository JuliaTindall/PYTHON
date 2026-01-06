#NAME
#    plot_d18op_from_pd.py
#PURPOSE 
#
#  This program will plot d18ov from pd files

# Import necessary libraries

import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import sys
#from netCDF4 import Dataset, MFDataset
#from mpl_toolkits.basemap import Basemap,maskoceans, shiftgrid

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")



#dumpfile = '/nfs/hera1/earjcti/um/xqetd/restart/xqetd.astart.nc'
dumpfile = '/nfs/hera1/earjcti/um/xqetd/restart/xqetda#da000005901c1+.nc'

plot_vapour = 'n'
plot_snow = 'y'

if plot_vapour == 'y':

    cube_16ov = iris.load_cube(dumpfile,'QT IN THE EXTERNAL DUMP')
    cube_18ov = iris.load_cube(dumpfile,'Stash code = 321')
    
    d18o_data = ((cube_18ov.data / cube_16ov.data)-2005.2E-6)/2005.2E-9
    cube_d18ov = iris.util.squeeze(cube_16ov.copy(data=d18o_data))
    
    plt.figure(figsize=[12,12])
    vals=np.arange(-75,75,15)
    for lev in range(0,9):
        plt.subplot(3,3,lev+1)
        qplt.contourf(cube_d18ov[lev,:,:],levels=vals,extend='both')
        plt.title=lev

    plt.show()
