#!/usr/bin/env python3

# created by Julia Tindall on 19/01/2024

# This program will get the january and july data from the file containing
# all the months
#  

import matplotlib.pyplot as plt
import numpy as np
import iris
import iris.quickplot as qplt
import sys

exptname = ['E280','E400','E560']
fieldname = ['NearSurfaceTemperature','TotalPrecipitation']


for expt in exptname:
    for field in fieldname:
        filename = ('/nfs/hera1/earjcti/temporary/Maisie/' + expt + '.' + field + '.mean_month.nc')
        cube = iris.load_cube(filename)
        jan_cube = cube.extract(iris.Constraint(month=1))
        jul_cube = cube.extract(iris.Constraint(month=7))
        dec_cube = cube.extract(iris.Constraint(month=12))
      
        print('saving January', expt, field)
        print(jan_cube)
        iris.save(jan_cube,expt + '.' + field +'_January.nc')
        print('saving July', expt, field)
        print(jul_cube)
        iris.save(jul_cube,expt + '.' + field +'_July.nc')
        
