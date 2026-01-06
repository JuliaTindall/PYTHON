#!/usr/bin/env python2.7
#NAME
#    PLOT_DIFFERENCES
#PURPOSE
#    This program will plot the differences between two files
# search for 'main program' to find end of functions
# Julia 11/1/2018



import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np




################################
# main program

# annual mean

#cntrl_expt='xkvje'
#plio_expt='xkvjf'
#anom_expt='xkvjg'
#extra='n'
#HadCM3='n'

cntrl_expt='xoorb'
#plio_expt='xibot'
anom_expt='xoorf'
HadCM3='y'
extra=' '
field = 'TotalCloud'


fieldtitle=field
if field =='field207':
    fieldtitle = 'lw_cs_toa'

cntl_file = ('/nfs/hera1/earjcti/um/' + cntrl_expt + '/database_averages/' + 
             cntrl_expt + '_Annual_Average_a@pd_' + field + '.nc')
anom_file = ('/nfs/hera1/earjcti/um/' + anom_expt + '/database_averages/' + 
             anom_expt + '_Annual_Average_a@pd_' + field + '.nc')

cntl_cube = iris.util.squeeze(iris.load_cube(cntl_file))
anom_cube = iris.util.squeeze(iris.load_cube(anom_file))

print(cntl_cube)
print(anom_cube)
diff_cube = anom_cube - cntl_cube
V=np.arange(-10.0,11.0,1.0)
if field == 'TotalCloud':
   V = V / 100.0
qplt.contourf(diff_cube,levels=V,extend='both',cmap='bwr')
plt.title(fieldtitle)
plt.gca().coastlines()
plt.show()
