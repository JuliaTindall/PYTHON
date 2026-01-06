#!/usr/bin/env python2.7
#NAME
#    enso_anom.py
#PURPOSE
#    This program will plot the enso anomaly for the Pliocene and the 
# preindustrial for figure 1 of the iso_enso paper
#
#  It will use the files produced by the IDL program
#    IDLPRGS/ISOTOPE/elnino/plot_pliocene_patterns_seasonal.pro
#
# which are IDLPLOTS/ISOTOPE/elnino/PATTERNS
#     a) xiboi_0_299_elnino_graph_ONI_ann.nc
#     b) xibol_0_299_elnino_graph_ONI_ann.nc
#
# Note that the program is very specific to the task and is not intended
# for generic application
#
#
# search for 'main program' to find end of functions
# Julia 21/4/2017



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid


#functions are:
#  def plotdata

# functions start here
def plotdata(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)

   # this is good for a tropical region
    map=Basemap(llcrnrlon=90.0,urcrnrlon=300.0,llcrnrlat=-45.0,urcrnrlat=45.0,projection='cyl',resolution='c')
   # this is good for the globe
   # map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
    map.drawmapboundary


    x, y = map(lons, lats)
    map.drawcoastlines()

    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='both')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='ar':
                    cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    print(np.shape(plotdata))
                    cs = map.contourf(x,y,plotdata,V)
                    cbar = plt.colorbar(cs,orientation="horizontal")




#   overplot some syms
    map.scatter(190,0,c='k',marker='o',s=300)
    map.scatter(190,0,c='g',marker='o',s=250)

    map.scatter(279,-7,c='k',marker='s',s=300)
    map.scatter(279,-7,c='g',marker='s',s=250)

    map.scatter(141,-3,c='k',marker='^',s=300)
    map.scatter(141,-3,c='g',marker='^',s=250)
   
    map.scatter(175,-16,c='k',marker='D',s=300)
    map.scatter(175,-16,c='g',marker='D',s=250)
  
    map.scatter(124,14.375,c='k',marker='*',s=400)
    map.scatter(124,14.375,c='g',marker='*',s=340)
   

    plt.title(titlename,loc='left',fontsize=25)
    cbar.ax.tick_params(labelsize=15)
    #cbar.set_label(cbarname,fontsize=15)
    cbar.ax.set_title(cbarname,fontsize=20)
#end def plotdata



################################
# main program

# open and read in preindstrial data

datasetname='/nfs/hera1/earjcti/IDLPLOTS/ISOTOPE/elnino/PATTERNS/xiboi_0_299_elnino_graph_ONI_ann.nc'
fpi=Dataset(datasetname)
latpd = fpi.variables['latitudepd'][:]
lonpd = fpi.variables['longitudepd'][:]
latpf = fpi.variables['latitudepf'][:]
lonpf = fpi.variables['longitudepf'][:]
pitemp=fpi.variables['a_Preind_temp'][:]
piprecip=fpi.variables['d_preind_precip'][:]
pid18op=fpi.variables['g_preind_d18op'][:]
pid18osw=fpi.variables['j_preind_d18osw'][:]
fpi.close()


# plot a) preindustrial temperature
degC=u'\N{DEGREE SIGN}'+'C'
plotdata(pitemp,0,lonpf,latpf,'a) Preindustrial EN-NT SST',-2.0,2.2,0.2,0.0,'ar',degC)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1a.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1a.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# plot d) preindustrial precipitation
plotdata(piprecip,0,lonpd,latpd,'d) Preindustrial EN-NT precipitation',-100.0,120.0,20.,0.0,'a','mm/month')
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1d.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1d.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# plot g) preindustrial d18op
permil=u'\u2030'
plotdata(pid18op,0,lonpd,latpd,'g) Preindustrial EN-NT d18O_p',-3.0,3.6,0.6,0.0,'ar',permil)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1g.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1g.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
 

# plot j) preindustrial d18osw
permil=u'\u2030'
plotdata(pid18osw,0,lonpf,latpf,'j) Preindustrial EN-NT d18O_sw',-0.15,0.18,0.03,0.0,'ar',permil)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1j.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1j.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
 

####################################
# open and read in Pliocene data

datasetname='/nfs/hera1/earjcti/IDLPLOTS/ISOTOPE/elnino/PATTERNS/xibol_0_299_elnino_graph_ONI_ann.nc'
fpi=Dataset(datasetname)
latpd = fpi.variables['latitudepd'][:]
lonpd = fpi.variables['longitudepd'][:]
latpf = fpi.variables['latitudepf'][:]
lonpf = fpi.variables['longitudepf'][:]
pltemp=fpi.variables['b_Plio_temp'][:]
plprecip=fpi.variables['e_plio_precip'][:]
pld18op=fpi.variables['h_plio_d18op'][:]
pld18osw=fpi.variables['k_plio_d18osw'][:]
anomtemp=fpi.variables['c_temp_anom'][:]
anomprecip=fpi.variables['f_anom_precip'][:]
anomd18op=fpi.variables['i_anom_d18op'][:]
anomd18osw=fpi.variables['l_anom_d18osw'][:]
fpi.close()


# plot b) mPWP temperature
degC=u'\N{DEGREE SIGN}'+'C'
plotdata(pltemp,0,lonpf,latpf,'b) mPWP EN-NT SST',-2.0,2.2,0.2,0.0,'ar',degC)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1b.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1b.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# plot e) mPWP precipitation
plotdata(plprecip,0,lonpd,latpd,'e) mPWP EN-NT precipitation',-100.0,120.0,20.,0.0,'a','mm/month')
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1e.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1e.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# plot h) mPWP d18op
permil=u'\u2030'
plotdata(pld18op,0,lonpd,latpd,'h) mPWP EN-NT d18O_p',-3.0,3.6,0.6,0.0,'ar',permil)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1h.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1h.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
 

# plot k) mPWP d18osw
permil=u'\u2030'
plotdata(pld18osw,0,lonpf,latpf,'k) mPWP EN-NT d18O_sw',-0.15,0.18,0.03,0.0,'ar',permil)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1k.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1k.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
 

##### plot anomalies

# plot c) mPWP temperature
degC=u'\N{DEGREE SIGN}'+'C'
plotdata(anomtemp,0,lonpf,latpf,'c) mPWP-Preindustrial [EN-NT] SST',-1.0,1.2,0.2,0.0,'ar',degC)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1c.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1c.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# plot f) mPWP-Preindustrialprecipitation
plotdata(anomprecip,0,lonpd,latpd,'f) mPWP-Preind [EN-NT] precip',-100.0,120.0,20.,0.0,'a','mm/month')
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1f.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1f.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# plot i) mPWP-Preindustriald18op
permil=u'\u2030'
plotdata(anomd18op,0,lonpd,latpd,'i) mPWP-Preind [EN-NT] d18O_p',-3.0,3.6,0.6,0.0,'ar',permil)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1i.eps' 
plt.savefig(fileout, bbox_inches='tight') 
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1i.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
 

# plot l) mPWP-Preindustrial d18osw
permil=u'\u2030'
plotdata(anomd18osw,0,lonpf,latpf,'l) mPWP-Preind [EN-NT] d18O_sw',-0.15,0.18,0.03,0.0,'ar',permil)
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1l.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/enso_anom/fig_1l.png' 
plt.savefig(fileout, bbox_inches='tight')  
plt.close()
 









sys.exit(0)

####

