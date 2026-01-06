# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a program that will create an orography netcdf file for input to HadCM3.
The orography will be pliomip2 everywhere except over n.america where it is p1
"""

import numpy as np
import iris
from mpl_toolkits.basemap import Basemap, shiftgrid, maskoceans
import matplotlib as mp
import matplotlib.pyplot as plt
import iris.quickplot as qplt
import iris.plot as iplt
import cartopy.crs as ccrs
import sys





def main():
    """
    1. read in the files
    2. Overwrite N America orography from P2 with P1
    3. write it out
    """
    cubenames = ['OROGRAPHY (/STRAT LOWER BC)',
                 'STANDARD DEVIATION OF OROGRAPHY',
                 'OROGRAPHIC GRADIENT XX COMPONENT',
                 'OROGRAPHIC GRADIENT XY COMPONENT',
                 'OROGRAPHIC GRADIENT YY COMPONENT',
                 'SILHOUETTE OROGRAPHIC ROUGHNESS',
                 'HALF OF  (PEAK TO TROUGH HT OF OROG)']

    # read in files and set up cubelist for output
    p1cubelist = iris.load(FILEINP1)
    print(p1cubelist)
    p2cubelist = iris.load(FILEINP2)
    print(p2cubelist)
    newcubelist = iris.cube.CubeList([])

    # loop through cubelist and replace p2 orography over N. America with
    # p1 orography over N. America.
    for name in cubenames:
        print(name)
        p1cube = iris.load_cube(FILEINP1,name)
        p2cube = iris.load_cube(FILEINP2,name)
        p1cubedata = p1cube.data
        p2cubedata = p2cube.data
        newcubedata = np.copy(p2cubedata)
 
        for j, lat in enumerate(p2cube.coord('latitude').points):
            if 10 < lat < 70:
                for i, lon in enumerate(p2cube.coord('longitude').points):
                    if 210 < lon < 300:
                        if (newcubedata[0,0,j,i] > 0 
                            and p1cubedata[0,0,j,i] > 0):
                            newcubedata[0, 0, j, i] = p1cubedata[0, 0, j, i]
                            
            
        newcube = p2cube.copy(data=newcubedata)
        newcubelist.append(newcube)
      # to check we have done it right
      #  if newcube.long_name == cubenames[0]:
      #      print('found it')
      #      V = np.arange(-400,440,40)
      #      ax1=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
      #      ax1.set_extent([-150, -60, 10, 70])
      #      anomcube = iris.util.squeeze(newcube - p1cube)
      #      cs = iplt.contourf(anomcube, V, extend='both')
      #      cbar = plt.colorbar(cs, orientation='horizontal')
      #      cbar.ax.tick_params(labelsize=7, rotation=90)
      #      plt.gca().coastlines()
      #      plt.title('new - P1 orog (m)',fontsize=10)
      #      plt.show()
      #  else:
      #      print(newcube.name,cubenames[0])
    iris.save(newcubelist, FILEOUT)     
       
    
  
   


FILEINP2 = '/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_qrparm.orog.nc'
FILEINP1 = '/nfs/hera1/earjcti/ancil/PRISM3/qrparm.orog.nc'
FILEOUT = '/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/Collaborators/N_America/P4_orog_P3_N_America.nc'



main()
