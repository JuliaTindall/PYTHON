#!/usr/bin/env python2.7
#NAME
#    PLOT_FLOW_FRAM_STRAIT
#PURPOSE
#   This program will find the volume of water flowing through the fram strait
#
# search for 'main program' to find end of functions
# Julia 20/02/2019


import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid
import pdb

#functions are:
#  def plotdata
#  def annmean_throughflow

# functions start here
def plotdata(region,plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,lonmin,lonmax,latreq):

    if region == 'NP':
        proj='npstere'
        latbb=45.
    if region == 'SP':
        proj='spstere'
        latbb=-45.
    #proj='stere'

    lons, lats = np.meshgrid(lon,lat)
    if fileno != 99:
        plt.subplot(2,2,fileno+1)

  
    map=Basemap(projection=proj,resolution='c',lon_0=0,boundinglat=latbb,round=True,lat_0=90)
  
    #map.drawmapboundary(fill_color='green')
    map.drawmapboundary

    x, y = map(lons, lats)

    map.drawcoastlines()
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal")
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu_r')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='max')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='i': #increasing
                    print(V)
                    cs = map.contourf(x,y,plotdata,V,norm=mp.colors.LogNorm(vmin=0,vmax=32),cmap='Reds')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    mycmap=mp.cm.get_cmap('Reds',len(V+2))
                    newcolors=mycmap(np.linspace(0,1,len(V+2)))
                    white=([1,1,1,1])
                    newcolors[0:2,:]=white
                    mycmap=ListedColormap(newcolors)
                    cs = map.contourf(x,y,plotdata,V,cmap=mycmap)
                   # map.drawparallels(np.arange(-80.,81.,20.))
                   # map.drawmeridians(np.arange(-180.,181.,20.))
                    map.fillcontinents()
                  
                    cbar = plt.colorbar(cs,orientation="horizontal")
                    # overplot in contours where control sea ice was over 50% (or 0.5)
  


    if fileno != 99:
        plt.title(titlename)
        cbar.set_label(cbarname,labelpad=-40)
    else:
        cbar.set_label(cbarname,labelpad=-70,size=15)
        cbar.ax.tick_params(labelsize=15)
        plt.title(titlename,loc='left',fontsize=15)
   
    lon_transect=np.linspace(lonmin,lonmax)
    lat_transect=np.linspace(latreq,latreq)
    x1,y1=map(lon_transect,lat_transect)
    map.plot(x1,y1,linewidth=1.5,color='yellow')
        


#end def plotdata

def annmean_throughflow(expt_list,extra_list,endlist,names_list,latreq,lonmin,lonmax):


    nexpts=len(expt_list)
 

    for expt in range(0,nexpts):
        f=MFDataset('/nfs/hera1/earjcti/um/'+expt_list[expt]+'/netcdf/'+expt_list[expt]+'o@pg'+extra_list[expt]+'[7-9]*'+endlist[expt]+'.nc')
        # get velocity
        lat = f.variables['latitude_1'][:]
        lon = f.variables['longitude_1'][:]
        dep = f.variables['depth_1'][:]
        atemp=f.variables['field704'][:] # velocity in cm/s
        atemp=np.squeeze(atemp)
        ntimes,nz,ny,nx=np.shape(atemp)

        # get potential temperature
        lat_theta = f.variables['latitude'][:]
        lon_theta = f.variables['longitude'][:]
        btemp=f.variables['temp'][:] # degC


        f.close()

        velocity_arr=np.mean(atemp,axis=0)
        theta_arr=np.mean(btemp,axis=0)+273.15  #(theta in kelvin)


        velocity_arr,lon = shiftgrid(180.,velocity_arr,lon,start=False)

            

        #plotdata('NP',velocity_arr[0,:,:],99,lon,lat,'velocity_map (v)',-5.0,5.5,0.5,0,'a','cm/s',lonmin,lonmax,latreq)
        #plt.show()
        #plt.close()

        # setup array for longitude height plot through transect
        if expt ==0:
            lon_req=np.where((lonmin<=lon) & (lonmax>=lon)) # note that lonreq are subscripts
                                                      # so lon[lonreq] will give the longitudes
           
            vel_trans=np.zeros((nexpts,nz,len(lon_req[0])))
            theta_trans=np.zeros((nexpts,nz,len(lon_req[0])))

        # put the data into the transect array
        lat_ix=np.where(lat==latreq)
        lonmin_ix=lon_req[0][0]
        lonmax_ix=lon_req[0][len(lon_req[0])-1]
     
        vel_trans[expt,:,:]=np.squeeze(velocity_arr[:,lat_ix[0],lonmin_ix:lonmax_ix+1])
        lon_trans=lon[lonmin_ix:lonmax_ix+1]

        # find out where the theta grid of latitude crosses the
        # velocity grid
        lat_theta_ix1=lat_ix[0]
        lat_theta_ix2=lat_ix[0]+1

        # check velocity lat is in middle of theta latitudes
        if lat_ix[0] > lat_theta_ix2[0] or lat_ix[0] < lat_theta_ix1[0]:
            print('your grid spacing is wrong')
            sys.exit()


        theta_diff=theta_arr[:,lat_theta_ix1,:]-theta_arr[:,lat_theta_ix2,:]
        theta_diff=np.squeeze(theta_diff)
        theta_diff,lon_theta = shiftgrid(180.,theta_diff,lon_theta,start=False)

        ncount=np.zeros(len(lon_trans))
        for i in range(0,len(lon_trans)):
            for i2 in range(0,len(lon_theta)):
                if np.abs(lon_theta[i2]-lon_trans[i]) < 1.0:
                    theta_trans[:,:,i]=theta_trans[:,:,i]+theta_diff[:,i2]
                    ncount[i]=ncount[i]+1

        theta_trans=theta_trans/ncount
    
   


       
    # plot velocity data

    for expt in range(0,nexpts):
        
        plt.subplot(2,3,expt+1)
        V=np.arange(-5.,5.5,0.5)
        cs=plt.contourf(lon_trans,dep,vel_trans[expt,:,:],V,extend='both',cmap='RdBu_r')
        plt.gca().invert_yaxis()
        plt.title(names_list[expt])
        plt.ylim(500,0)
        plt.ylabel('depth (m)')

        plt.tick_params(axis='both',labelsize=8)
    
        locs=np.arange(-20,30,10)
        labs2=[]
        for i in range(0,len(locs)):
            if locs[i]<0:
                labs2.append(np.str(locs[i]*(-1))+'W')
            else:
                labs2.append(np.str(locs[i])+'E')
        
   
        plt.xticks(locs,labs2)
    
    
        if expt==1:
            plt.subplot(2,3,6)
            plt.gca().set_visible(False)
            cbar = plt.colorbar(cs,orientation="vertical",fraction=0.7)
            cbar.set_label('cm/s',rotation=0,fontsize=10,labelpad=20)


    plt.tight_layout()
    
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/timeslices/plot_flow_fram_strait/flow_fram_strait.eps' 
    plt.savefig(fileout, bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/timeslices/plot_flow_fram_strait/flow_fram_strait.png' 
    plt.savefig(fileout, dpi=300)
  
    plt.close()

    # plot anomalies from pi

    for expt in range(1,nexpts):
        vel_tran_anom=vel_trans[expt,:,:]-vel_trans[0,:,:]
        plt.subplot(2,2,expt)
        V=np.arange(-2.,2.2,0.2)
        cs=plt.contourf(lon_trans,dep,vel_tran_anom,V,extend='both',cmap='RdBu_r')
        plt.gca().invert_yaxis()
        plt.title(names_list[expt]+'minus'+names_list[0])

    plt.subplots_adjust(bottom=0.25)
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    plt.colorbar(cax=cax,orientation="horizontal")


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/timeslices/plot_flow_fram_strait/flow_fram_strait_anom_pi.eps' 
    plt.savefig(fileout, bbox_inches='tight')
  
    plt.close()


    # plot anomalies from km5c

    for expt in range(2,nexpts):
        vel_tran_anom=vel_trans[expt,:,:]-vel_trans[1,:,:]
        plt.subplot(2,2,expt)
        V=np.arange(-1.,1.1,0.1)
        cs=plt.contourf(lon_trans,dep,vel_tran_anom,V,extend='both',cmap='RdBu_r')
        plt.gca().invert_yaxis()
        plt.title(names_list[expt]+'minus'+names_list[1])

    plt.subplots_adjust(bottom=0.25)
    cax = plt.axes([0.1, 0.1, 0.8, 0.05])
    plt.colorbar(cax=cax,orientation="horizontal")


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/timeslices/plot_flow_fram_strait/flow_fram_strait_anom_km5c.eps' 
    plt.savefig(fileout, bbox_inches='tight')
  
    plt.close()


    # #############################################################
    # calculate total transport through fram strait in top 500m only
    
    
    # calculate size of depths
    boundary_depth=np.zeros(len(dep)+1)
    boundary_depth[0]=0.
    boundary_depth[len(dep)]=dep[len(dep)-1]
    for z in range(1,len(dep)):
        boundary_depth[z]=(dep[z]+dep[z-1])/2.
    
    layer_size=boundary_depth[1:]-boundary_depth[:-1]
   
    # calculate transport for top 500m
    
    flux=np.zeros((nexpts,len(lon_trans)))
    north_flux=np.zeros(nexpts)
    south_flux=np.zeros(nexpts)
   
    onegridbox = (lon[1]-lon[0])*111000.*np.cos(np.deg2rad(latreq))
    area=layer_size * onegridbox
    print(onegridbox,lon[1],lon[0])
    print(area)
    
    # plot flux within 50m of surface
    plt.subplot(2,1,1)
    for expt in range(0,nexpts):
        for z in range(0,len(dep)): 
            if dep[z] < 50.:
                waterflux=vel_trans[expt,z,:] * area[z] / 100.
                flux[expt,:]=flux[expt,:]+waterflux
                for i in range(0,len(lon_trans)):
                    if waterflux[i] > 0:
                        north_flux[expt]=north_flux[expt]+waterflux[i]
                    else:
                        south_flux[expt]=south_flux[expt]+waterflux[i]
                
      
        #legname=(names_list[expt]+':'+np.str(np.round(north_flux[expt]/1.0E6,2))+
        #         'Sv N :'+np.str(np.round(south_flux[expt]/1.0E6,2))+'Sv S')
        legname=(names_list[expt]+':'+np.str(np.round(np.sum(flux[expt,:])/1.0E6,2))+'Sv')
        plt.plot(lon_trans,flux[expt,:]/1.0E6,label=legname)
    
    plt.plot([-100,180],[0,0],color='black')
    print(lon_req[0])
 
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylabel('Sv')
    #plt.xlabel('longitude')
    locs=np.arange(-20,30,10)
    labs2=[]
    for i in range(0,len(locs)):
        if locs[i]<0:
            labs2.append(np.str((locs[i]*(-1)))+'W')
        else:
            labs2.append(np.str(locs[i])+'E')
        print(np.str(locs[i]))
    
   
    plt.xticks(locs,labs2)
    
    plt.title('flux in top 50m')
    plt.xlim(np.min(lon_trans),np.max(lon_trans))
    
    
    # plot flux within 500m of surface
    flux=np.zeros((nexpts,len(lon_trans))) 
    north_flux=np.zeros(nexpts)
    south_flux=np.zeros(nexpts)
    plt.subplot(2,1,2)
    for expt in range(0,nexpts):
        for z in range(0,len(dep)): 
            if dep[z] < 500.:
                waterflux=vel_trans[expt,z,:] * area[z] / 100.
                flux[expt,:]=flux[expt,:]+waterflux
                for i in range(0,len(lon_trans)):
                    if waterflux[i] > 0:    
                        north_flux[expt]=north_flux[expt]+waterflux[i]
                    else:
                        south_flux[expt]=south_flux[expt]+waterflux[i]
        print(names_list[expt],np.sum(flux[expt,:]/1.0E6))
        #legname=(names_list[expt]+':'+np.str(np.round(north_flux[expt]/1.0E6,2))+
        #         'Sv N :'+np.str(np.round(south_flux[expt]/1.0E6,2))+'Sv S')
        legname=(names_list[expt]+':'+np.str(np.round(np.sum(flux[expt,:])/1.0E6,2))+'Sv')
        plt.plot(lon_trans,flux[expt,:]/1.0E6,label=legname)
    
    plt.plot([-100,180],[0,0],color='black')
    print(lon_req[0])
 
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.ylabel('Sv')
    #plt.xlabel('longitude')
    
    locs=np.arange(-20,30,10)
    labs2=[]
    for i in range(0,len(locs)):
        if locs[i]<0:
            labs2.append(np.str((locs[i]*(-1)))+'W')
        else:
            labs2.append(np.str(locs[i])+'E')
        print(np.str(locs[i]))
    
   
    plt.xticks(locs,labs2)
    
    plt.title('flux in top 500m')
    plt.xlim(np.min(lon_trans),np.max(lon_trans))
    plt.tight_layout()
  
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/timeslices/plot_flow_fram_strait/waterflux.eps' 
    plt.savefig(fileout, bbox_inches='tight')
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/timeslices/plot_flow_fram_strait/waterflux.png' 
    plt.savefig(fileout, dpi=300)
    plt.close()

    
    
   
    

#

################################
# main program

# annual mean
figureno=0

# timeslices are xiboi=preindustrial, xibol=3205 - km5c', xjplc=3205-km5c, xjpld=3060 (K1), xjple=2950 (G17), xjplf=3155 (KM3)


expt_list=['xiboi','xjplc','xjpld','xjple','xjplf']
extra_list=['y','6','6','6','6']
end_list=['c1','11','11','11','11']
names_list=['PreInd','KM5C','K1','G17','KM3']

#djf mean
season=['dc','ja','fb']
seasonname='djf'

# seasonmean throughflow will do a longitude depth plot of velocity at a 
# given latitude

#Fraim strait
latreq=78.75   # fraim strait is between 77N and 81N
lonmin=338.75 - 360. # this is Greenland
lonmax=30.       # Svalbard is about 15E but I don't think it is in the UM

#latreq=30.0
#lonmin=-180.
#lonmax=180.
annmean_throughflow(expt_list,extra_list,end_list,names_list,latreq,lonmin,lonmax)

#

####

