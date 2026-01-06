#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created 16 May 2023
#
#@author: earjcti
#"""
#
# Niels noted that for his site the summer warming was more than the
# winter warming.
# He wondered whether this was a global signal or just a local signal.
# This program is to see what the models show
#
#
#


import cartopy as cart
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib as mpl
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import iris.plot as iplt
import iris.coord_categorisation
import numpy as np
import sys

       

def get_data(exptname,fieldname):
    """
    gets the data if it is a single model
    gets the multimodel mean if it is multimodels
    """

    if MODEL == 'MMMsubset':
        for i, model in enumerate(MMMss_mods):
            cube = iris.load_cube(filestart + model + '/' + exptname + 
                                  '.' + fieldname + '_850.0.mean_month.nc')
            iris.util.squeeze(cube)
            if i == 0:
                print(np.shape(cube.data))
                nt, ny,nx = np.shape(cube.data)
                dataarr = np.zeros((len(MMMss_mods), nt, ny, nx))
                
            dataarr[i,:,:,:] = cube.data
      
        dataarr2 = np.where(dataarr > 1.0E10, np.nan, dataarr)
        meandata = np.mean(dataarr2,axis=0)
        meancube = cube.copy(data=meandata)

    return meancube

############
# MAIN PROGRAM
################


# read in multimodel mean monthly winds
MMMss_mods = ['CCSM4-UoT','COSMOS','EC-Earth3.3','GISS2.1G',
              'HadCM3','IPSLCM6A','MIROC4m','NorESM1-F']
#MMMss_mods = ['GISS2.1G']
MODEL = 'MMMsubset'

filestart = '/nfs/hera1/earjcti/regridded/'


EOI400_ua = get_data('EOI400','ua')
EOI400_va = get_data('EOI400','va')
E280_ua = get_data('E280','ua')
E280_va = get_data('E280','va')


EOI400_jan_ua = EOI400_ua.extract(iris.Constraint(month=1))
EOI400_jul_ua = EOI400_ua.extract(iris.Constraint(month=7))
EOI400_jan_va = EOI400_va.extract(iris.Constraint(month=1))
EOI400_jul_va = EOI400_va.extract(iris.Constraint(month=7))

E280_jan_ua = E280_ua.extract(iris.Constraint(month=1))
E280_jul_ua = E280_ua.extract(iris.Constraint(month=7))
E280_jan_va = E280_va.extract(iris.Constraint(month=1))
E280_jul_va = E280_va.extract(iris.Constraint(month=7))

x = E280_jan_ua.coord('longitude').points
y = E280_jan_ua.coord('latitude').points
ulon = E280_jan_ua.coord('longitude')

fig = plt.figure(figsize=[12,12])

ax = plt.subplot(3,2,1,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-80,35,20,90],crs=ccrs.PlateCarree())
#transform = ulon.coord_system.as_cartopy_projection()
Q=plt.quiver(x[::5],y[::5],
           E280_jan_ua.data[::5,::5], E280_jan_va.data[::5,::5], 
           pivot="middle",scale=50.00, scale_units='inches') 
qk=ax.quiverkey(Q, 0.9, 1.05, 10, r'$10\frac{m}{s}$',  #need to sort out quiver key
                labelpos='E', coordinates='axes')
          # transform=transform)
plt.title('Jan PI')
ax.coastlines()



ax = plt.subplot(3,2,3,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-80,35,20,90],crs=ccrs.PlateCarree())
#transform = ulon.coord_system.as_cartopy_projection()
Q=plt.quiver(x[::5],y[::5],
           EOI400_jan_ua.data[::5,::5], EOI400_jan_va.data[::5,::5], 
           pivot="middle",scale=50.00, scale_units='inches') 
qk=ax.quiverkey(Q, 0.9, 1.05, 10, r'$10\frac{m}{s}$',
                labelpos='E', coordinates='axes')
          # transform=transform)
plt.title('Jan MP')
ax.coastlines()

ax = plt.subplot(3,2,5,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-80,35,20,90],crs=ccrs.PlateCarree())
#transform = ulon.coord_system.as_cartopy_projection()
Q=plt.quiver(x[::5],y[::5],
           EOI400_jan_ua.data[::5,::5] - E280_jan_ua.data[::5,::5],
           EOI400_jan_va.data[::5,::5] - E280_jan_va.data[::5,::5], 
           pivot="middle",scale=20.00, scale_units='inches') 
qk=ax.quiverkey(Q, 0.9, 1.05, 10, r'$10\frac{m}{s}$',
                labelpos='E', coordinates='axes')
          # transform=transform)
plt.title('Jan MP - PI')
ax.coastlines()


ax = plt.subplot(3,2,2,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-80,35,20,90],crs=ccrs.PlateCarree())
#transform = ulon.coord_system.as_cartopy_projection()
Q=plt.quiver(x[::5],y[::5],
           E280_jul_ua.data[::5,::5], E280_jul_va.data[::5,::5], 
           pivot="middle",scale=50.00, scale_units='inches') 
qk=ax.quiverkey(Q, 0.9, 1.05, 10, r'$10\frac{m}{s}$',
                labelpos='E', coordinates='axes')
          # transform=transform)
plt.title('Jul PI')
ax.coastlines()



ax = plt.subplot(3,2,4,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-80,35,20,90],crs=ccrs.PlateCarree())
#transform = ulon.coord_system.as_cartopy_projection()
Q=plt.quiver(x[::5],y[::5],
           EOI400_jul_ua.data[::5,::5], EOI400_jul_va.data[::5,::5], 
           pivot="middle",scale=50.00, scale_units='inches') 
qk=ax.quiverkey(Q, 0.9, 1.05, 10, r'$10\frac{m}{s}$',
                labelpos='E', coordinates='axes')
          # transform=transform)
plt.title('Jul MP')
ax.coastlines()


ax = plt.subplot(3,2,6,projection=ccrs.PlateCarree(central_longitude=0.0))
ax.set_extent([-80,35,20,90],crs=ccrs.PlateCarree())
#transform = ulon.coord_system.as_cartopy_projection()
Q=plt.quiver(x[::5],y[::5],
           EOI400_jul_ua.data[::5,::5] - E280_jul_ua.data[::5,::5], 
           EOI400_jul_va.data[::5,::5] - E280_jul_va.data[::5,::5], 
           pivot="middle",scale=20.00, scale_units='inches') 
qk=ax.quiverkey(Q, 0.9, 1.05, 10, r'$10\frac{m}{s}$',
                labelpos='E', coordinates='axes')
          # transform=transform)
plt.title('Jul MP - PI')
ax.coastlines()




plt.savefig('summer_and_winter_winds_850mb' + MODEL + '.png')
plt.savefig('summer_and_winter_winds_850mb' + MODEL + '.eps')
