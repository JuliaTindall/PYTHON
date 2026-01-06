 #!/usr/bin/env python2.7
#NAME
#    PLOT_OCEAN HEAT TRANSPORT
#PURPOSE
#    This program will do all the plots to do with ocean heat transport
#   
#    We want to see how much heat the THC transported in the pliocene and if 
#    if this is different from the modern
# search for 'main program' to find end of functions
# Julia 19/1/2017



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
#  def annmean
#  def seasmean

# functions start here
def plotmap(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(2,2,fileno+1)


   # this is good for a tropical region
   # map=Basemap(llcrnrlon=10.0,urcrnrlon=70.0,llcrnrlat=10.0,urcrnrlat=55.0,projection='cyl',resolution='c')
   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
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
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu_r')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='max')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                cs = map.contourf(x,y,plotdata,V)
                cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
 
#end def plotdata
#####################################
def plot_lat_dep(plotdata,fileno,lat,dep,titlename,minval,maxval,valinc,cbarname):
    lats, deps  = np.meshgrid(lat,dep)

    V=np.arange(minval,maxval,valinc)
    
    cs = plt.contourf(lats,deps,plotdata,V,extend="both",cmap='RdBu_r')
    plt.gca().invert_yaxis()

    cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
 
#end def plotdata

#  to check if a character is numeric
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#============================
def Atlantic_heat_transport(expt_name):



#  get data from files

    f=MFDataset('/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pg2/'+expt_name+'o@pg*.nc')
    lat = f.variables['latitude_1'][:]
    lon = f.variables['longitude_1'][:]
    dep = f.variables['depth_1'][:]
    Vvel=f.variables['field704'][:]
    PotT=f.variables['temp'][:]
    Vvel=np.squeeze(Vvel)
    PotT=np.squeeze(PotT)
    f.close()

    ntimes,ndep,nlat,nlon=np.shape(Vvel)
    
    # calculate T(j)-T(j+1)  which is notthwards anomaly
    
    PotT_diff=np.zeros((ntimes,ndep,nlat,nlon))
    for j in range(0,nlat):
        PotT_diff[:,:,j,:]=PotT[:,:,j,:]-PotT[:,:,j+1,:] 

   
    vT=Vvel * PotT_diff
    mean_VT=np.mean(vT,axis=0)
    meanPotT_diff=np.mean(PotT_diff,axis=0)

  

   

#  get mask of Oceans
    filemask='/nfs/see-fs-02_users/earjcti/MOC/merid/basin_hadgom_216'
    f1=open(filemask,'r')
    # discard 3 title lines
    textline=f1.readline()
    textline=f1.readline()
    textline=f1.readline()
    # next line is details from basin
    nrows=nlat+1  # because it's on T grid not V grid
    Indian_d=np.zeros((nrows,4))
    Pacific_d=np.zeros((nrows,4))
    Atlantic_d=np.zeros((nrows,4))
    Combined_d=np.zeros((nrows,4))
    for line in f1:
        linesplit=line.split()  # split line by space
        if is_number(linesplit[0]):
            rowno=int(linesplit[0])
            Indian_d[rowno-1,:]=int(linesplit[1]),int(linesplit[2]),\
                int(linesplit[3]),int(linesplit[4])
            Pacific_d[rowno-1,:]=int(linesplit[5]),int(linesplit[6]),\
                int(linesplit[7]),int(linesplit[8])
            Atlantic_d[rowno-1,:]=int(linesplit[9]),int(linesplit[10]),\
                int(linesplit[11]),int(linesplit[12])
            Combined_d[rowno-1,:]=int(linesplit[13]),int(linesplit[14]),\
                int(linesplit[15]),int(linesplit[16])

    f1.close()

    # make some attempt to convert to V grid
    Atlantic_V=np.zeros((nlat,4))
    for j in range(0,nlat):
        Atlantic_V[j,:]=(Atlantic_d[j,:]+Atlantic_d[j+1,:])/2.0


    Atlantic_mask=np.ones((nlat,nlon),dtype=bool)
    for j in range(0,nlat):
        Atlantic_mask[j,Atlantic_V[j,0]:Atlantic_V[j,1]]=0
        Atlantic_mask[j,Atlantic_V[j,2]:Atlantic_V[j,3]]=0



    # get mean heat transport over Atlantic
    Atlantic_VT=np.ma.masked_array(mean_VT,mask=np.tile(Atlantic_mask,(mean_VT.shape[0],1)))


    Atlantic_VT=Atlantic_VT / 100. # convert from K cm/s to K m/s


    lontemp=lon
    PotT_diff,lon = shiftgrid(180.,PotT_diff,lon,start=False)
    lon=lontemp
    mean_VT,lon = shiftgrid(180.,mean_VT,lon,start=False)
    lon=lontemp
    Atlantic_VT, lon = shiftgrid(180.,Atlantic_VT,lon,start=False)
    lon=lontemp
    meanPotT_diff, lon = shiftgrid(180.,meanPotT_diff,lon,start=False)


    #plotmap(meanPotT_diff[0,:,:],0,lon,lat,'theta diff',-2,2,0.1,0.0,'n','degC')
    #plotmap(mean_VT[0,:,:],1,lon,lat,'meanvt',-2,2,0.1,0.0,'n','degC')
    #plt.show()
    
 
    #plotmap(Atlantic_VT[0,:,:],2,lon,lat,'atl_vt0',-2,2,0.1,0.0,'n','degC')
    #plotmap(Atlantic_VT[30,:,:],3,lon,lat,'atlvt30',-2,2,0.1,0.0,'n','degC')
    #plt.show()


    # we have extracted the data we want.  Now have a go at
    # calculating surfact heat transport.
    # (This is Cp * integral [(dens * VT)]dx


    dens=1029 # in reality this will change depending on temperature and 
              # saltiness  but mean value is 1029kg/m3
    cp = 3985      # This is seawater heat capacity 3985 J /kg / K

    Atl_qsurf=np.zeros((ndep,nlat))

    for k in range(0,ndep):
        for j in range(0,nlat):
            dx=(lon[1]-lon[0]) * 111000. * np.cos(np.radians(lat[j]))
            Atl_weight=Atlantic_VT[k,j,:]*cp*dens*dx
            Atl_qsurf[k,j]=np.sum(Atl_weight)
            if k == 0 and lat[j] < 5. and lat[j]> -5.:
                print(Atlantic_VT[k,j,:])
                print('pott',meanPotT_diff[k,j,:])
                print('vvel',np.mean(Vvel[:,k,j,:],axis=0))
                print(dx)
                sys.exit()

    # units of Atl_qsurf[j] are currently J/s/m.  This is because we have 
    # not accounted for depth.  To find the total heat transport we would 
    # multiply by the depth of the layer which would give the results in Watts

    dep_lay=np.zeros(ndep)
    dep_lay[0]=dep[0] * 2.0
    dep_lay[ndep-1]=dep[ndep-1] -dep[ndep-2]
    for k in range(1,ndep-1):
        dep_lay[k]=(dep[k+1]-dep[k-1])/2.0
   
   
    Atl_qsurf=Atl_qsurf / 10**15  # give results in petawatts
    Atl_qsurf=Atl_qsurf * dep_lay[:,None]
    print('shape',np.shape(Atl_qsurf))


    titlename='Atlantic Heat transport: '+expt_name
    plot_lat_dep(Atl_qsurf,0,lat,dep,titlename,-1.5,1.5,0.1,'Petawatts')


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_ocean_heat_transport/'+expt_name+'_Atl_lat_depth.eps'
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

    return(lat,dep,Atl_qsurf)
    

#end def Atlantic_salinity_depth_plot
#end def annmean

#========================================
def lat_heat_transport(expt_name1,expt_OHT1,expt_name2,expt_OHT2,lat,dep):
# this subroutine will do a plot of ocean heat transport integrated over all the ocean 

    print('in lat heat transport')
    allOHT1=np.nansum(expt_OHT1,axis=0)
    plt.plot(lat,allOHT1)
    plt.plot(lat,expt_OHT1[0,:],'b.')

    if expt_name2 != 'none':
        allOHT2=np.nansum(expt_OHT2,axis=0)
        plt.plot(lat,allOHT2,'r-')
        plt.plot(lat,expt_OHT2[0,:],'r.')




    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_ocean_heat_transport/'+expt_name1+'_'+expt_name2+'_Atl_lat.eps'
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()
    print('end lat heat transport')



#end def lat_heat_transport

################################
# main program

# annual mean salinity by depth

retdata=Atlantic_heat_transport('xkvjg')
lat=retdata[0]
dep=retdata[1]
xkvjg_OHT=retdata[2]

retdata=Atlantic_heat_transport('xkvjf')
xkvjf_OHT=retdata[2]

retdata=Atlantic_heat_transport('xkvje')
xkvje_OHT=retdata[2]

# do an anomaly plot

anomplot=xkvjg_OHT-xkvje_OHT
titlename='Atlantic Heat transport: xkvjg-xkvje'
plot_lat_dep(anomplot,0,lat,dep,titlename,-0.5,0.5,0.05,'Petawatts')
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_ocean_heat_transport/xkvjg-xkvje_Atl_lat_depth.eps'
plt.savefig(fileout, bbox_inches='tight')  
plt.close()


anomplot=xkvjf_OHT-xkvje_OHT
titlename='Atlantic Heat transport: xkvjf-xkvje'
plot_lat_dep(anomplot,0,lat,dep,titlename,-0.5,0.5,0.05,'Petawatts')
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_ocean_heat_transport/xkvjf-xkvje_Atl_lat_depth.eps'
plt.savefig(fileout, bbox_inches='tight')  
plt.close()

# integrate over ocean depth and then plot
print('calling lat heat transport')

lat_heat_transport('xkvjg',xkvjg_OHT,'xkvje',xkvje_OHT,lat,dep)


sys.exit()



####

