#!/usr/bin/env python2.7
#NAME
#    PLOT_D18O
#PURPOSE
#    This program will plot some d18o data over the Middle east for modern
#     
#    The programs used to create the netcdf files are:
#    1.  IDLPRGS/AMAZON/ GNIP_SPATIAL_COMPARISON
#    2.  IDLPRGS/8.2ka/JONATHAN/average_field_from_timeseries.pro
#    3.                         average_d18o_from_timeseries.pro
# Julia 2/11/2016



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import MFDataset, Dataset


# functions start here
def plotdata(plotdata,fileno,lons,lats,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(2,2,fileno+1)

    from mpl_toolkits.basemap import Basemap
    map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary
    x, y = map(lons, lats)
    map.drawcoastlines()
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal",extend='max')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-256,vmax=256),cmap='RdBu')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='max')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                cs = map.contourf(x,y,plotdata,V, extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
    #cbar.set_label('test label',verticalalignment='bottom')
    
 #   if fileno ==0:
 #       cax=plt.axes([0.1,0.55,0.4,0.02])
 #   if fileno ==1:
 #       cax=plt.axes([0.55,0.55,0.4,0.02])
 #   if fileno ==2:
 #       cax=plt.axes([0.1,0.05,0.4,0.02])
 #   if fileno ==3:
 #       cax=plt.axes([0.55,0.05,0.4,0.02])
   # plt.colorbar(cax=cax,orientation="horizontal")


################################
# main program

# read in the data from Kanhu's modern simulation

os.chdir("/nfs/hera1/earjcti/um/netcdf/xiboi_netcdf/")
filename="xiboia@pdy11*.nc"
f=MFDataset(filename)
f.dimensions
f.variables

# retrieving a variable
stashcode338 = f.variables['QCL'][:]
lat = f.variables['latitude'][:]
lon = f.variables['longitude'][:]

f.close()
print(np.shape(stashcode338))

mean338 = np.mean(stashcode338,axis=0)
total16o = mean338[0,:,:]+mean338[3,:,:]+mean338[6,:,:]+mean338[9,:,:]
total18o = mean338[1,:,:]+mean338[4,:,:]+mean338[7,:,:]+mean338[10,:,:]
delta18o = (((total18o/total16o)-2005.2e-6))/2005.2e-9
print (delta18o)
#print (np.shape(mean338))
permille=u'\u2030'
plotdata(delta18o,0,lon,lat,'delta18o',-15,15,1.0,0.0,'n',permille)
plt.show()
sys.exit()

# plot modern data

degC=u'\N{DEGREE SIGN}'+'C'
permille=u'\u2030'


plt.figure(0)
plotdata(avg_d18o_mod,0,lon,lat,'avg d18O modern',-15,0,1.0,0.0,'n',permille)
plotdata(stdev_d18o_mod,1,lon,lat,'stdev d18O modern',0,7,1.0,0.0,'n',permille)

V=[0,2,4,8,16,32,64,128,256]
plotdata(avg_precip_mod,2,lon,lat,'avg precip modern',0,100,10.0,V,'y','mm/month')
V=[0,2,4,8,16,32,64,128,256,512]
plotdata(stdev_precip_mod,3,lon,lat,'stdev precip modern',0,100,10,V,'y','mm/month')

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/0ka_d18o_p.png' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/0ka_d18o_p.eps' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()
plt.figure(1)


plotdata(avg_temp_mod,0,lon,lat,'avg temp modern',0,40,5.0,0,'n',degC)
plotdata(stdev_temp_mod,1,lon,lat,'stdev temp modern',0,1.5,0.1,0,'n',degC)

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/0ka_temperature.png' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/0ka_temperature.eps' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()

#=========================================
# read in the data from the 9ka simulation

os.chdir("/nfs/hera1/earjcti/IDLPLOTS/8.2ka/")

filename="tcnkx_d18o_300_839.nc"
f=Dataset(filename, mode='r')
f.dimensions
f.variables
avg_d18o_9k = f.variables['d18o - annual'][:]
stdev_d18o_9k = f.variables['d18o - stdev'][:]
f.close()

filename="tcnkx_tempsurf_300_839.nc"
f=Dataset(filename, mode='r')
f.dimensions
f.variables
avg_temp_9k = f.variables['tempsurf mean'][:]
stdev_temp_9k = f.variables['tempsurf stdev'][:]
f.close()

filename="tcnkx_precip_300_839.nc"
f=Dataset(filename, mode='r')
f.dimensions
f.variables
avg_precip_9k = f.variables['precip mean mmmonth'][:]
stdev_precip_9k = f.variables['precip stdev mmyear)'][:]
f.close()

# plot 9ka data

plt.figure(0)
plotdata(avg_d18o_9k,0,lon,lat,'avg d18O 9ka',-15,0,1.0,0.0,'n',permille)
plotdata(stdev_d18o_9k,1,lon,lat,'stdev d18O 9ka',0,7,1.0,0.0,'n',permille)



V=[0,2,4,8,16,32,64,128,256]
avg_precip_9k=avg_precip_9k * 60. * 60. * 24. * 30.
plotdata(avg_precip_9k,2,lon,lat,'avg precip 9ka',0,100,10.0,V,'y','mm/month')
V=[0,2,4,8,16,32,64,128,256,512]
plotdata(stdev_precip_9k,3,lon,lat,'stdev precip 9ka',0,100,10,V,'y','mm/month')

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka_d18o_p.png' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka_d18o_p.eps' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()
plt.figure(1)


plotdata(avg_temp_9k-273.15,0,lon,lat,'avg temp 9ka',0,40,5.0,0,'n',degC)
plotdata(stdev_temp_9k,1,lon,lat,'stdev temp 9ka',0,1.5,0.1,0,'n',degC)

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka_temperature.png' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka_temperature.eps' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()



# plot anomaly data data


plt.figure(0)
plotdata(avg_d18o_9k-avg_d18O_mod,0,lon,lat,'avg d18O 9ka - modern',-5,5.5,0.5,0.0,'a',permille)
plotdata(stdev_d18o_9k-stdev_d18o_mod,1,lon,lat,'stdev d18O 9ka - modern',-3,3.5,0.5,0.0,'a',permille)

V=[-128,-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64,128]
plotdata(avg_precip_9k-avg_precip_mod,2,lon,lat,'avg precip 9ka-modern',-50,60,10.0,V,'la','mm/month')
V=[-128,-64,-32,-16,-8,-4,-2,0,2,4,8,16,32,64,128]
plotdata(stdev_precip_9k-stdev_precip_mod,3,lon,lat,'stdev precip 9ka - modern',-100,120,20,V,'la','mm/month')

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka-mod_d18o_p.png' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka-mod_d18o_p.eps' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()
plt.figure(1)


plotdata(avg_temp_9k-273.15-avg_temp_mod,0,lon,lat,'avg temp 9ka - modern',-10,12,2.0,0,'a',degC)
plotdata(stdev_temp_9k-stdev_temp_mod,1,lon,lat,'stdev temp 9ka - modern',-1.0,1.2,0.2,0,'a',degC)

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka-mod_temperature.png' 
plt.savefig(fileout, bbox_inches='tight')  

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/MIDDLE_EAST/plot_anomalies/9ka-mod_temperature.eps' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()







import sys
sys.exit()

for plotno in range(0,6):
    if plotno % 2 == 0:
        monthname='jn'
    else:
        monthname='dc'
   
    #fileno=(plotno // 2) + 6 for xmrda
    if (plotno // 2) == 0:
        fileno='07'
    if (plotno // 2) == 1:
        fileno='13'
    if (plotno // 2) == 2:
        fileno='19'


    print(fileno,monthname)

    filename="xmrdaa@pd6%s%s.nc" % (fileno,monthname)
 
    print(filename)



    lons, lats = np.meshgrid(lon,lat)

    #X=np.zeros(lat.shape,lon.shape)
    X=np.zeros((lev.shape[0],lat.shape[0],lon.shape[0]))

    X[:,:,:] = np.squeeze(sm_data[0,:,:,:])

    # level is 0

    levuse=1
    sm=X[levuse,:,:] # single level only
    sm[sm>1e10]=0.0 # set missing values to zero

    plt.subplot(3,2,plotno+1)
    from mpl_toolkits.basemap import Basemap
    map=Basemap(projection='robin',resolution='l',lat_0=0,lon_0=180)
    map.drawmapboundary(fill_color='aqua')
    x, y = map(lons, lats)
    map.drawcoastlines()
    #if levuse == 0:
    #   V=np.arange(0,50,5)
    #if levuse == 3:
    V=np.arange(np.amin(sm),np.amax(sm),1)
       
    cs = map.contourf(x,y,sm,V)
    plt.title(filename)
    

    

plt.subplots_adjust(bottom=0.2)
cax=plt.axes([0.1,0.1,0.8,0.045])
plt.colorbar(cax=cax,orientation='horizontal')
plt.title('kg/m2')


fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/xmrda_smc_%s.eps' %levuse
plt.savefig(fileout, bbox_inches='tight')  

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/xmrda_smc_%s.png' %levuse
plt.savefig(fileout, bbox_inches='tight')  




