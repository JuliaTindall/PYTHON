import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

filename = ('/uolstore/Research/b/b0164/Data/PlioMIP3/AWI/AWI-ESM2/' +
            'PI/AWI-ESM2_PI_totalprecip_monthly_timeseries_native_645001-694912.nc')

cube = iris.load_cube(filename)
avgcube = cube.collapsed('time',iris.analysis.MEAN)

vals=np.arange(0,30,1)
qplt.pcolormesh(avgcube, vmin=0, vmax=6, cmap='viridis')
#qplt.contourf(avgcube,levels=vals,extend='both')
plt.gca().coastlines()
plt.show()

