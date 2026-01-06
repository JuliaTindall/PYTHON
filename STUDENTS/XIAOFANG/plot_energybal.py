#!/usr/bin/env python2.7
#NAME
#    PLOT_ENERGYBAL
#PURPOSE
#    This program will plot the energy balance for Xiaofangs simulation
# xoorb is the standard mid pliocene experiment, xoorf is the one with
# the reduced height of the antarctic ice sheet
#
# search for 'main program' to find end of functions
# Julia 11/1/2018



import os
import numpy as np
import scipy as sp
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
#from mpl_toolkits.basemap import Basemap, shiftgrid


#functions are:
#  def global_enbal
#  def seasmean


def global_enbal(exptname,HadCM3):

    #==============
    # read in data from  average temperature files produced for Dan

    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/'+exptname+'/database_averages/'+exptname
    else:
        filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/database_averages/'+exptname

    print(filestart)

    f=Dataset(filestart+'_Annual_Average_a@pd_Temperature.nc')
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['temp'][:]
    atemp=np.squeeze(atemp)
    ny,nx=np.shape(atemp)

    # get upward and downward sw radiation
    # incoming sw
    f=Dataset(filestart+'_Annual_Average_a@pd_field200.nc')
    sw_down=f.variables['field200'][:]
    sw_down=np.squeeze(sw_down)

    # outgoing sw (clear sky flux)
    f=Dataset(filestart+'_Annual_Average_a@pd_field207.nc')
    sw_up=f.variables['field207'][:]
    sw_up=np.squeeze(sw_up)
 

   # outgoing lw  (toa)
    filename='olr'
    f=Dataset(filestart+'_Annual_Average_a@pd_olr.nc')
    lw_toa=f.variables['olr'][:]
    lw_toa=np.squeeze(lw_toa)
  
    # outgoing lw  (toa clear sky)
    f=Dataset(filestart+'_Annual_Average_a@pd_csolr.nc')
    lwcs_toa=f.variables['field207'][:]
    lwcs_toa=np.squeeze(lwcs_toa)
    

    # net downward longwave surface
    filename='longwave'
    f=Dataset(filestart+'_Annual_Average_a@pd_'+filename+'.nc')
    net_lwdown_surf=f.variables['longwave'][:]
    net_lwdown_surf=np.squeeze(net_lwdown_surf)

 
    # outgoing lw  (surface)
    filename='ilr'
    f=Dataset(filestart+'_Annual_Average_a@pd_'+filename+'.nc')
    lw_surf=f.variables['ilr'][:]
    lw_surf=np.squeeze(lw_surf)

    
    # JULIA NOTE  THIS IS MUCH WORSE WHEN SUBTRACTING NET LWDOWN SURF
    lw_surf=lw_surf-net_lwdown_surf # upward lw rad is downlw - net downlw
    
    # ====================================
    # get the global average fields

    grid_alpha=sw_up/sw_down
    grid_epsilon=lw_toa/lw_surf

    # create weighting array
    weightarr=np.zeros(np.shape(atemp))
    for i in range(0,len(lon)):
        weightarr[:,i]=np.cos(np.deg2rad(lat))

    meantemp=np.average(atemp,weights=weightarr)
    mean_sw_up=np.average(sw_up,weights=weightarr)
    mean_sw_down=np.average(sw_down,weights=weightarr)
    mean_lw_toa=np.average(lw_toa,weights=weightarr)
    mean_lwcs_toa=np.average(lwcs_toa,weights=weightarr)
    mean_lw_surf=np.average(lw_surf,weights=weightarr)
    mean_alpha=np.average(grid_alpha, weights=weightarr)
    mean_epsilon=np.average(grid_epsilon, weights=weightarr)

    #============================================
    # calculate terms in equation
    So=1367
    sigma=5.67 * (10.0 ** (-8.))

    alpha=mean_sw_up / mean_sw_down
    mean_lw_surf_up=sigma * (meantemp ** 4)

    epsilon=mean_lw_toa / mean_lw_surf # julia using mean lwcs_toa is worse


    print('alphas',alpha,mean_alpha)
    print('epsilon',epsilon,mean_epsilon)

    t4=So / 4.0 * (1.0-alpha) / (epsilon * sigma)
    t=t4 ** (1./4.)

    t4_mean=So / 4.0 * (1.0-mean_alpha) / (mean_epsilon * sigma)
    t_mean=t4_mean ** (1./4.)


    print('globvals',So/4.0,mean_alpha,mean_epsilon,sigma)



    print('    ')
    print(exptname)
    print('=========')
    print('t from formula=',t,' K   ',t-273.15,' degC')
    print('t obs=',meantemp,' K   ',meantemp-273.15,' degC')
    print('t mean=',t_mean,' K   ',t_mean-273.15,' degC')

    # julia note the mean energy balance is not quite right.  I am ignoring it because the latitudinal energy balance seems fine.

#end def global_enbal

################################################
def dh_zonal_enbal(exptname,HadCM3):

  
    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/'+exptname+'/database_averages/'+exptname
    else:
        filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/database_averages/'+exptname


    #==============
    # read in data from  average temperature files produced for Dan

    f=Dataset(filestart+'_Annual_Average_a@pd_Temperature.nc')
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['temp'][:]
    atemp=np.squeeze(atemp)
    ny,nx=np.shape(atemp)

    # get upward and downward sw radiation
    # incoming sw
    f=Dataset(filestart+'_Annual_Average_a@pd_field200.nc')
    sw_down=f.variables['field200'][:]
    sw_down=np.squeeze(sw_down)

    # outgoing sw 
    f=Dataset(filestart+'_Annual_Average_a@pd_field201.nc')
    sw_up=f.variables['field201'][:]
    sw_up=np.squeeze(sw_up)
 
    # outgoing sw (clear sky flux)
    f=Dataset(filestart+'_Annual_Average_a@pd_field207.nc')
    swcs_up=f.variables['field207'][:]
    swcs_up=np.squeeze(swcs_up)
 

   # outgoing lw  (toa)
    f=Dataset(filestart+'_Annual_Average_a@pd_olr.nc')
    lw_toa_up=f.variables['olr'][:]
    lw_toa_up=np.squeeze(lw_toa_up)
  
    # outgoing lw  (toa clear sky)
    f=Dataset(filestart+'_Annual_Average_a@pd_csolr.nc')
    lwcs_toa=f.variables['field207'][:]
    lwcs_toa=np.squeeze(lwcs_toa)
    

    # net downward longwave surface
    f=Dataset(filestart+'_Annual_Average_a@pd_longwave.nc')
    net_lwdown_surf=f.variables['longwave'][:]
    net_lwdown_surf=np.squeeze(net_lwdown_surf)
 
 
    # incoming lw  (surface)
    f=Dataset(filestart+'_Annual_Average_a@pd_ilr.nc')
    lw_surf_down=f.variables['ilr'][:]
    lw_surf_down=np.squeeze(lw_surf_down)


    lw_surf_up=lw_surf_down-net_lwdown_surf # upward lw rad is downlw - net downlw


    print(net_lwdown_surf)    
    
    # ====================================
    # get the zonal average fields


    grid_alpha=sw_up/sw_down
    grid_alpha_cs=swcs_up/sw_down
    grid_epsilon=lw_toa_up/(lw_surf_up)
    grid_epsilon_cs=lwcs_toa / lw_surf_up

   
    meantemp=np.average(atemp,axis=1)
    mean_sw_up=np.average(sw_up,axis=1)
    mean_sw_down=np.average(sw_down,axis=1)
    mean_lw_toa_up=np.average(lw_toa_up,axis=1)
    mean_lwcs_toa_up=np.average(lwcs_toa,axis=1)
    mean_lw_surf_down=np.average(lw_surf_down,axis=1)
    mean_net_lwdown_surf=np.average(net_lwdown_surf,axis=1)
    mean_alpha=np.average(grid_alpha,axis=1)
    mean_alpha_cs=np.average(grid_alpha_cs,axis=1)
    mean_epsilon=np.average(grid_epsilon,axis=1)
    mean_epsilon_cs=np.average(grid_epsilon_cs,axis=1)
    mean_lw_surf_up=np.average(lw_surf_up,axis=1)

    plt.subplot(2,1,1)
    plt.plot(lat,mean_alpha,label='mean_alpha')
    plt.plot(lat,mean_alpha_cs,label='mean_alpha_cs')
    plt.legend()
    plt.title('different alphas')



    sigma=5.67 * (10.0 ** (-8.))

    #============================================
    # calculate terms in equation
    H=(-1.0) * ((mean_sw_down - mean_sw_up) -(mean_lw_toa_up))
   
    
    print('means',np.average(mean_sw_down),np.average(mean_sw_up), np.average(mean_net_lwdown_surf),np.average(mean_lw_toa_up))

    t4=((mean_sw_down * (1.0-mean_alpha)) + H) / (mean_epsilon * sigma)
    t=t4 ** (1./4.)

    plt.subplot(2,1,2)
    plt.plot(lat,meantemp-t,label='temperature')
    plt.legend()
    plt.title('temperature difference')
   
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/dh_temp_anom'+anom_expt+'-'+cntrl_expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()


    
    components=[mean_epsilon,mean_epsilon_cs,mean_alpha,mean_alpha_cs,H,t,mean_sw_down,lat]
    return(components)

#end def dh_zonal_enbal

def split_gge(deltaT_gge):
    """
    splits up the gge and topography for the DH energy balance
    """
    orog_fname_cntl ='/nfs/hera1/earjcti/Xiaofang/P4_enh_qrparm.orog.nc'
    orog_fname_expt = '/nfs/hera1/earjcti/Xiaofang/P4_enh_qrparm.orog_no_antarctica.nc'

    f=Dataset(orog_fname_cntl)
    cntl_orog=f.variables['ht'][:]
    cntl_orog=np.squeeze(cntl_orog)
    cntl_orog_lat=np.average(cntl_orog,axis=1)
   

    f=Dataset(orog_fname_expt)
    expt_orog=f.variables['ht'][:]
    expt_orog=np.squeeze(expt_orog)
    expt_orog_lat=np.average(expt_orog,axis=1)
 
    deltaT_topo=(expt_orog_lat - cntl_orog_lat) * (-5.5) / 1000.
    deltaT_gge_only = deltaT_gge - deltaT_topo

    return deltaT_topo, deltaT_gge_only

def rf_zonal_enbal(exptname,HadCM3):

    if HadCM3 == 'y':
        filestart='/nfs/hera1/earjcti/um/'+exptname+'/database_averages/'+exptname
    else:
        filestart='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/database_averages/'+exptname


    #==============
    # read in data from  average temperature files produced for Dan

    f=Dataset(filestart+'_Annual_Average_a@pd_Temperature.nc')
    lat = f.variables['latitude'][:]
    lon = f.variables['longitude'][:]
    atemp=f.variables['temp'][:]
    atemp=np.squeeze(atemp)
    ny,nx=np.shape(atemp)

    # get upward and downward sw radiation
    # incoming sw
    f=Dataset(filestart+'_Annual_Average_a@pd_field200.nc')
    sw_down=f.variables['field200'][:]
    sw_down=np.squeeze(sw_down)

    # outgoing sw 
    f=Dataset(filestart+'_Annual_Average_a@pd_field201.nc')
    sw_up=f.variables['field201'][:]
    sw_up=np.squeeze(sw_up)
 
    # surface upwards sw (clear sky flux)
    f=Dataset(filestart+'_Annual_Average_a@pd_field207_1.nc')
    sw_surf_cs_up=f.variables['field207_1'][:]
    sw_surf_cs_up=np.squeeze(sw_surf_cs_up)
 
    # surface downwards sw (clear sky flux)
    f=Dataset(filestart+'_Annual_Average_a@pd_field208.nc')
    sw_surf_cs_down=f.variables['field208'][:]
    sw_surf_cs_down=np.squeeze(sw_surf_cs_down)
 
  # outgoing sw (clear sky flux)
    f=Dataset(filestart+'_Annual_Average_a@pd_field207.nc')
    swcs_up=f.variables['field207'][:]
    swcs_up=np.squeeze(swcs_up)
 

   # outgoing lw  (toa)
    f=Dataset(filestart+'_Annual_Average_a@pd_olr.nc')
    lw_toa_up=f.variables['olr'][:]
    lw_toa_up=np.squeeze(lw_toa_up)
  
    # outgoing lw  (toa clear sky)
    f=Dataset(filestart+'_Annual_Average_a@pd_csolr.nc')
    lwcs_toa=f.variables['field207'][:]
    lwcs_toa=np.squeeze(lwcs_toa)
    

    # net downward longwave surface
    f=Dataset(filestart+'_Annual_Average_a@pd_longwave.nc')
    net_lwdown_surf=f.variables['longwave'][:]
    net_lwdown_surf=np.squeeze(net_lwdown_surf)

 
    # incoming lw  (surface)
    f=Dataset(filestart+'_Annual_Average_a@pd_ilr.nc')
    lw_surf_down=f.variables['ilr'][:]
    lw_surf_down=np.squeeze(lw_surf_down)


    # net downward shortwave flux surface
    f=Dataset(filestart+'_Annual_Average_a@pd_solar.nc')
    net_swdown_surf=f.variables['solar'][:]
    net_swdown_surf=np.squeeze(net_swdown_surf)

    # total downward shortwave flux surface
    f=Dataset(filestart+'_Annual_Average_a@pd_field203.nc')
    sw_surf_down=f.variables['field203'][:]
    sw_surf_down=np.squeeze(sw_surf_down)


    # total cloud amount random overlap
    f=Dataset(filestart+'_Annual_Average_a@pd_TotalCloud.nc')
    cloud_frac=f.variables['field30'][:]
    cloud_frac=np.squeeze(cloud_frac)


    lw_surf_up=lw_surf_down-net_lwdown_surf # upward lw rad is downlw - net downlw

    sw_surf_up=sw_surf_down - net_swdown_surf # upwards sw rad at surf

 
    # topography
    if exptname=='xkvje':
        orog_fname='/nfs/hera1/earjcti/um/HadGEM_ancils/qrparm.orog_new.nc'
    if exptname=='xkvjg':
        orog_fname='/nfs/hera2/apps/metadata/ancil/PRISM3_ALT/HadGEM2/HadGEM_pliocene_orog.nc'
    if exptname=='xibos':
        orog_fname='/nfs/hera2/apps/metadata/ancil/cntrl2/qrparm.orog.nc'
    if exptname=='xibot':
        orog_fname='/nfs/hera2/apps/metadata/ancil/PRISM3_ALT/qrparm.orog.nc'
    f=Dataset(orog_fname)
    orog=f.variables['ht'][:]
    orog=np.squeeze(orog)


    
    # ====================================
    # get gridded epsilon and alpha



    grid_epsilon=lw_toa_up/(lw_surf_up)
    grid_epsilon_cs=lwcs_toa / lw_surf_up


    # from Taylor 2007 
    # alpha is now surface albedo.  We denote planetary albedo:  A.
    # alpha=(1-c)*alpha_clr + (c*alpha_oc)
    # oc=overcast, clr=clear sky, c=cloud fraction

    grid_A=sw_up/sw_down
    grid_A_clr=swcs_up/sw_down
    grid_alpha_clr=sw_surf_cs_up / sw_surf_cs_down
    grid_alpha=sw_surf_up / sw_surf_down
    grid_alpha_oc=(grid_alpha - ((1-cloud_frac) * grid_alpha_clr)) / cloud_frac


    #===============================
    # get amount scattered (gamma) and 
    # amount not absorbed (mu)  (absorbtion =1-mu)
    # see Taylor 2007 equation 9 for the equations

    grid_mu=(sw_surf_down / sw_down) * (1.0 - grid_alpha) + grid_A
    grid_mu_clr=(sw_surf_cs_down / sw_down) * (1.0 - grid_alpha_clr)+grid_A_clr
    grid_mu_oc=(grid_mu - ((1-cloud_frac) * grid_mu_clr)) / cloud_frac
    grid_mu_cld=grid_mu_oc / grid_mu_clr

    grid_gamma=((grid_mu - (sw_surf_down / sw_down)) /
                (grid_mu - (grid_alpha * (sw_surf_down / sw_down))))
    grid_gamma_clr=((grid_mu_clr - (sw_surf_cs_down / sw_down)) /
                (grid_mu_clr - (grid_alpha_clr * (sw_surf_cs_down / sw_down))))
    grid_gamma_oc=(grid_gamma - ((1-cloud_frac) * grid_gamma_clr)) / cloud_frac
    grid_gamma_cld=1.0 - ((1.0 - grid_gamma_oc)/(1.0 - grid_gamma_clr))
                 

    # note these values of gamma, alpha and mu have been checked against
    # equations 7 and 8 of Taylor 2007.  Therefore if A is correct and
    # qs_hat_down is correct then gamma alpha and mu are also correct.  

  
    #===============================
    # get fields for use in energy balance calculation
 
    meantemp=np.average(atemp,axis=1)
    mean_sw_up=np.average(sw_up,axis=1)
    mean_sw_down=np.average(sw_down,axis=1)
    mean_sw_surf_cs_down=np.average(sw_surf_cs_down,axis=1)
    mean_sw_surf_down=np.average(sw_surf_down,axis=1)
    mean_lw_toa_up=np.average(lw_toa_up,axis=1)
    mean_lwcs_toa_up=np.average(lwcs_toa,axis=1)
    mean_lw_surf_down=np.average(lw_surf_down,axis=1)
    mean_net_lwdown_surf=np.average(net_lwdown_surf,axis=1)
    mean_alpha=np.average(grid_alpha,axis=1)
    mean_alpha_clr=np.average(grid_alpha_clr,axis=1)
    mean_alpha_oc=np.average(grid_alpha_oc,axis=1)
    mean_alpha=np.average(grid_alpha,axis=1)
    mean_epsilon=np.average(grid_epsilon,axis=1)
    mean_epsilon_cs=np.average(grid_epsilon_cs,axis=1)
    mean_lw_surf_up=np.average(lw_surf_up,axis=1)
    mean_cloud=np.average(cloud_frac,axis=1)
    mean_mu=np.average(grid_mu,axis=1)
    mean_mu_clr=np.average(grid_mu_clr,axis=1)
    mean_mu_oc=np.average(grid_mu_oc,axis=1)
    mean_mu_cld=np.average(grid_mu_cld,axis=1)
    mean_gamma=np.average(grid_gamma,axis=1)
    mean_gamma_clr=np.average(grid_gamma_clr,axis=1)
    mean_gamma_oc=np.average(grid_gamma_oc,axis=1)
    mean_gamma_cld=np.average(grid_gamma_cld,axis=1)
    mean_A=np.average(grid_A,axis=1)
    mean_orog=np.average(orog,axis=1)

    #plt.subplot(3,1,1)
    #plt.plot(lat,mean_alpha-0.1,label='mean_alpha-0.1')
    #plt.plot(lat,mean_alpha_clr,label='mean_alpha_clr')
    #plt.plot(lat,mean_alpha_oc,label='mean_alpha_oc')
    #plt.legend()
    #plt.title('different values')

    #plt.subplot(3,1,2)
    #plt.plot(lat,mean_mu,label='mean_mu')
    #plt.plot(lat,mean_mu_clr,label='mean_mu_clr')
    #plt.plot(lat,mean_mu_oc,label='mean_mu_oc')
    #plt.legend()


    #plt.subplot(3,1,3)
    #plt.plot(lat,mean_gamma,label='mean_gamma')
    #plt.plot(lat,mean_gamma_clr,label='mean_gamma_clr')
    #plt.plot(lat,mean_gamma_oc,label='mean_gamma_oc')  
    #plt.legend()

    #plt.show()
   


    sigma=5.67 * (10.0 ** (-8.))

    #============================================
    # calculate terms in equation
    H=(-1.0) * ((mean_sw_down - mean_sw_up) -(mean_lw_toa_up))
   
    

    t4=((mean_sw_down * (1.0-mean_alpha)) + H) / (mean_epsilon * sigma)
    t=t4 ** (1./4.)


    components=[mean_epsilon,mean_epsilon_cs,mean_alpha_clr,mean_alpha_oc,mean_mu_cld,mean_gamma_cld,mean_mu_clr,mean_gamma_clr,mean_cloud,H,meantemp,mean_sw_down,lat,mean_alpha,mean_mu,mean_gamma,mean_mu_oc,mean_gamma_oc,mean_A,mean_orog]
    return(components)

#end def rf_zonal_enbal


##########################
def main_dh_zonal(cntrl_expt,anom_expt,extra,HadCM3):

    icesheetred = {'xoorf' : '0%', 'xoorg' : '25%', 'xoorh': '50%', 'xoori': '75%'}

    components=dh_zonal_enbal(cntrl_expt,HadCM3)
    emis_pi=components[0]
    emis_pi_cs=components[1]
    alpha_pi=components[2]
    alpha_pi_cs=components[3]
    H_pi=components[4]
    temp_pi=components[5]
    sw_down_pi=components[6]
    lat=components[7]
    
    components=dh_zonal_enbal(anom_expt,HadCM3)
    emis_plio=components[0]
    emis_plio_cs=components[1]
    alpha_plio=components[2]
    alpha_plio_cs=components[3]
    H_plio=components[4]
    temp_plio=components[5]
    sw_down_plio=components[6]
    lat=components[7]

    sigma=5.67 * (10.0**-8)

  

    # we are now decomposing but will also find the average value from 55N-90N
    # create weighting array
    weightarr=np.cos(np.deg2rad(lat))
    for i in range(0,len(lat)):
        if lat[i] < 55:
            weightarr[i]=0.

    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_plio_cs * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_pi_cs * sigma)
    t_2=t4 ** (1./4.)

    deltaT_gge=t_1-t_2

    deltaT_height, deltaT_gge_only = split_gge(deltaT_gge)
    #plt.subplot(2,1,1)
    plt.plot(lat,deltaT_gge,label='Topography + GHG',color="cyan")
    #plt.plot(lat,deltaT_height,label='topography',color="blue")
    #print(deltaT_height)
    #plt.plot(lat,deltaT_gge_only,label='GHG',color="pink")
    print('Arctic T gge',np.average(deltaT_gge,weights=weightarr))


    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_plio * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_plio_cs * sigma)
    t_2=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_pi * sigma)
    t_3=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_pi_cs * sigma)
    t_4=t4 ** (1./4.)

    deltaT_ce=(t_1-t_2) - (t_3 - t_4)
    print('Arctic T cloud emisivity',np.average(deltaT_ce,weights=weightarr))
    plt.plot(lat,deltaT_ce, label='cloud emissivity',color="orange")

    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_plio * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_plio_cs)) + H_plio) / (emis_plio * sigma)
    t_2=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_pi)) + H_plio) / (emis_plio * sigma)
    t_3=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_pi_cs)) + H_plio) / (emis_plio * sigma)
    t_4=t4 ** (1./4.)

    deltaT_ca=(t_1-t_2) - (t_3 - t_4)
    plt.plot(lat,deltaT_ca, label='cloud albedo',color="purple")
    print('Arctic T cloud albedo',np.average(deltaT_ca,weights=weightarr))



    t4=((sw_down_plio * (1.0-alpha_plio_cs)) + H_plio) / (emis_plio * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_pi_cs)) + H_plio) / (emis_plio * sigma)
    t_2=t4 ** (1./4.)

    deltaT_csa=t_1-t_2
    plt.plot(lat,deltaT_csa,label='clear sky albedo',color="chocolate",linestyle="dashdot")
    print('Arctic T clear sky albedo',np.average(deltaT_csa,weights=weightarr))



    t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_plio * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio * (1.0-alpha_plio)) + H_pi) / (emis_plio * sigma)
    t_2=t4 ** (1./4.)

    deltaT_H=t_1-t_2
    plt.plot(lat,deltaT_H,label='heat transport',color="red")
    print('Arctic T Heat transport',np.average(deltaT_H,weights=weightarr))
    plt.legend()

    deltaT=temp_plio-temp_pi
    plt.plot(lat,deltaT,label='Total temperature change')
    print('Arctic T total change',np.average(deltaT,weights=weightarr))
    plt.legend()


    plt.ylim(-2.0,22.0)
    mp.pyplot.axhline(y=0,xmin=-90,xmax=90,color='black')
#    if HadCM3 == 'y':
#        plt.title('Energy balance HadCM3 - Dan Hill')
#    else:
#        plt.title('Energy balance HadGEM2- Dan Hill')
        
    plt.xlabel('latitude')
    #plt.ylabel('Warming ' + anom_expt + ' - ' + cntrl_expt)
    ylabel = 'Warming:' + icesheetred.get(anom_expt) + 'EAIS - MPControl, degC'
    plt.ylabel(ylabel)
    
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/DH_energybal'+anom_expt+'-'+cntrl_expt
    plt.savefig(fileout + '.eps', bbox_inches='tight')  
    plt.savefig(fileout + '.png', bbox_inches='tight')  

    plt.close()

    f = open(fileout + '.txt','w+')
    line = ('lat, Topo+GHG, cloud emissivity, cloud albedo, clear sky albedo,Total temperature change \n')
    f.write(line)

    for i, latuse in enumerate(lat):
        print(latuse)
        line = (str(latuse) + ',' + str(np.round(deltaT_gge[i],4)) + ',' +  str(np.round(deltaT_ce[i],4))+ ',' +  str(np.round(deltaT_ca[i],4))+ ',' +  str(np.round(deltaT_csa[i],4))+ ',' +  str(np.round(deltaT_H[i],4))+ ',' +  str(np.round(deltaT[i],4))+ '\n')
        f.write(line)
    f.close()

        

# end of main part of Dan Hills energy balance



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
anom_expt='xoorf' # xoorf - Antarctica 0%
                  # xoorg - Antarctica 25%
                  # xoorh - Antarctica 50%
                  # xoori - Antarctica 75%
HadCM3='y'
extra=' '


global_enbal(cntrl_expt,HadCM3)
global_enbal(anom_expt,HadCM3)


# Dan Hills energy balance
main_dh_zonal(cntrl_expt,anom_expt,extra,HadCM3)


# Ran Fengss energy balance


components=rf_zonal_enbal(cntrl_expt,HadCM3)
emis_pi=components[0]
emis_pi_cs=components[1]
alpha_pi_clr=components[2]
alpha_pi_oc=components[3]
mu_pi_cld=components[4]
gamma_pi_cld=components[5]
mu_pi_clr=components[6]
gamma_pi_clr=components[7]
cloud_pi=components[8]
H_pi=components[9]
temp_pi=components[10]
sw_down_pi=components[11]
lat=components[12]
alpha_pi=components[13]
mu_pi=components[14]
gamma_pi=components[15]
mu_pi_oc=components[16]
gamma_pi_oc=components[17]
pi_A=components[18]
topo_pi=components[19]


components=rf_zonal_enbal(anom_expt,HadCM3)
emis_plio=components[0]
emis_plio_cs=components[1]
alpha_plio_clr=components[2]
alpha_plio_oc=components[3]
mu_plio_cld=components[4]
gamma_plio_cld=components[5]
mu_plio_clr=components[6]
gamma_plio_clr=components[7]
cloud_plio=components[8]
H_plio=components[9]
temp_plio=components[10]
sw_down_plio=components[11]
lat=components[12]
alpha_plio=components[13]
mu_plio=components[14]
gamma_plio=components[15]
mu_plio_oc=components[16]
gamma_plio_oc=components[17]
plio_A=components[18]
topo_plio=components[19]


delta_T_topo=(topo_plio-topo_pi) * (-5.5) / 1000.
sigma=5.67 * (10.0**-8)

# print j corresponding to 75deg
for j in range(0,len(lat)):
    if lat[j]==75:
        lat75=j
        print(j,lat[j])

#taylor 2007 equation 16 a-c and 7
# note the equation 16a is misleading.  I think it should be
#delta_A_alpha=(1-c)delta_A_alpha_clr + c*delta_A_alpha_oc 

# get the change in albedo due to alpha


# new use mean values
mu=(mu_plio+mu_pi)/2.0
gamma=(gamma_plio+gamma_pi)/2.0
alpha=(alpha_plio+alpha_pi)/2.0
cloud=(cloud_plio+cloud_pi)/2.0

A_alpha_pi_clr=((mu * gamma) + ((mu * alpha_pi_clr * (1.0-gamma)**2.0)/ (1.0 - alpha_pi_clr * gamma)))
A_alpha_plio_clr=((mu * gamma) + ((mu * alpha_plio_clr * (1.0-gamma)**2.0)/ (1.0 - alpha_plio_clr * gamma)))
A_alpha_pi_oc=((mu * gamma) + ((mu * alpha_pi_oc * (1.0-gamma)**2.0)/ (1.0 - alpha_pi_oc * gamma)))
A_alpha_plio_oc=((mu * gamma) + ((mu * alpha_plio_oc * (1.0-gamma)**2.0)/ (1.0 - alpha_plio_oc * gamma)))
A_alpha_diff=(1.0-cloud)*(A_alpha_plio_clr - A_alpha_pi_clr)+cloud*(A_alpha_plio_oc - A_alpha_pi_oc)
A_alpha_plio=(1.0-cloud)*A_alpha_plio_clr + cloud*A_alpha_plio_oc
A_alpha_pi=(1.0-cloud)*A_alpha_pi_clr + cloud*A_alpha_pi_oc
# julia print out at 75N
print('surface alpha at 75',A_alpha_diff[lat75],A_alpha_plio[lat75],A_alpha_pi[lat75])

# get the change in albedo due to clouds eqn 16b


A_cld_pi_gamma=((mu * gamma_pi_cld) + ((mu * alpha * (1.0-gamma_pi_cld)**2.0)/ (1.0 - alpha * gamma_pi_cld)))
A_cld_plio_gamma=((mu * gamma_plio_cld) + ((mu * alpha * (1.0-gamma_plio_cld)**2.0)/ (1.0 - alpha * gamma_plio_cld)))
A_cld_pi_mu=((mu_pi_cld * gamma) + ((mu_pi_cld * alpha * (1.0-gamma)**2.0)/ (1.0 - alpha * gamma)))
A_cld_plio_mu=((mu_plio_cld * gamma) + ((mu_plio_cld * alpha * (1.0-gamma)**2.0)/ (1.0 - alpha * gamma)))

# equation 15 A=(a-c)Aclr + c Aoc

mu_clr=(mu_plio_clr + mu_pi_clr)/2.0
gamma_clr=(gamma_plio_clr + gamma_pi_clr)/2.0
alpha_clr=(alpha_plio_clr + alpha_pi_clr)/2.0
mu_oc=(mu_plio_oc + mu_pi_oc)/2.0
gamma_oc=(gamma_plio_oc + gamma_pi_oc)/2.0
alpha_oc=(alpha_plio_oc + alpha_pi_oc)/2.0

clravg=((mu_clr * gamma_clr) + ((mu_clr * alpha_clr * (1.0-gamma_clr)**2.0)/ (1.0 - alpha_clr * gamma_clr)))
ocavg=((mu_oc * gamma_oc) + ((mu_oc * alpha_oc * (1.0-gamma_oc)**2.0)/ (1.0 - alpha_oc * gamma_oc)))
A_deltaC=((cloud_pi-cloud_plio)*clravg)+((cloud_plio-cloud_pi)*ocavg)


A_cld_diff=A_cld_plio_mu-A_cld_pi_mu + A_cld_plio_gamma-A_cld_pi_gamma + A_deltaC

# julia print out at 75N
print('cloud A at 75',A_cld_diff[lat75])



# get the change in albedo due to clear skies


A_clr_pi_mu=((mu_pi_clr * gamma) + ((mu_pi_clr * alpha * (1.0-gamma)**2.0)/ (1.0 - alpha * gamma)))
A_clr_plio_mu=((mu_plio_clr * gamma) + ((mu_plio_clr * alpha * (1.0-gamma)**2.0)/ (1.0 - alpha * gamma)))
A_clr_pi_gamma=((mu * gamma_pi_clr) + ((mu * alpha * (1.0-gamma_pi_clr)**2.0)/ (1.0 - alpha * gamma_pi_clr)))
A_clr_plio_gamma=((mu * gamma_plio_clr) + ((mu * alpha * (1.0-gamma_plio_clr)**2.0)/ (1.0 - alpha * gamma_plio_clr)))
A_clr_diff=(A_clr_plio_mu - A_clr_pi_mu)+(A_clr_plio_gamma - A_clr_pi_gamma)
print('clr A at 75',A_clr_diff[lat75],A_clr_plio_mu[lat75],A_clr_pi_mu[lat75],A_clr_plio_gamma[lat75],A_clr_pi_gamma[lat75],cloud_pi[lat75],cloud_plio[lat75])


# get change in planetary albedo to check budgets

A_diff=plio_A-pi_A


#plt.subplot(2,1,1)
plt.plot(lat,A_alpha_diff,label='A_alpha')
plt.plot(lat,A_cld_diff,label='A_cld')
plt.plot(lat,A_clr_diff,label='A clr diff')
plt.plot(lat,A_diff,label='total planetary Albedo diff')
plt.plot(lat,A_alpha_diff+A_cld_diff+A_clr_diff,label='sum')
plt.legend()
plt.title('plio-pi A components')

fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_albedo'+anom_expt+'-'+cntrl_expt+'.eps' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()


#print values at 75N in plot
print('A_alpha_diff',A_alpha_diff[lat75])
print('A_cld_diff',A_cld_diff[lat75])
print('A_clr_diff',A_clr_diff[lat75])
print('A_diff',A_diff[lat75])
print(' ')



# we are now decomposing but will also find the average value from 55N-90N
# create weighting array
weightarr_Arctic=np.cos(np.deg2rad(lat))
for i in range(0,len(lat)):
    if lat[i] < 55:
        weightarr_Arctic[i]=0.
weightarr=np.cos(np.deg2rad(lat))

t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_plio_cs * sigma)
t_1=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-alpha_plio)) + H_plio) / (emis_pi_cs * sigma)
t_2=t4 ** (1./4.)

deltaT_gge=t_1-t_2 - delta_T_topo


#plt.subplot(2,1,2)
fig=plt.figure()
ax=plt.subplot(111)
ax.plot(lat,deltaT_gge,label='Greenhouse gas emissivity',color="blue")
#plt.plot(lat,delta_T_topo,label='topography')
print('Arctic T gge',np.average(deltaT_gge,weights=weightarr_Arctic))
print('Arctic T topo',np.average(delta_T_topo,weights=weightarr_Arctic))


t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_plio * sigma)
t_1=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_plio_cs * sigma)
t_2=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_pi * sigma)
t_3=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_pi_cs * sigma)
t_4=t4 ** (1./4.)

deltaT_ce=(t_1-t_2) - (t_3 - t_4)
print('Arctic T cloud emisivity',np.average(deltaT_ce,weights=weightarr_Arctic))
ax.plot(lat,deltaT_ce, label='cloud emissivity',color="orange")

# surface albedo

t4=((sw_down_plio * (1.0-A_alpha_plio)) + H_plio) / (emis_plio * sigma)
t_1=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-A_alpha_pi)) + H_plio) / (emis_plio * sigma)
t_2=t4 ** (1./4.)

deltaT_surfalpha=(t_1-t_2) 
ax.plot(lat,deltaT_surfalpha, label='surface albedo',color="chocolate",linestyle='dotted')
print('Arctic T surface albedo',np.average(deltaT_surfalpha,weights=weightarr_Arctic))

# check alternate value for surface (this seems to work fine)
# we did this to check how we had done clear sky and cloud
t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_plio * sigma)
t_1=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-(plio_A-A_alpha_diff))) + H_plio) / (emis_plio * sigma)
t_2=t4 ** (1./4.)

deltaT_altsurf=t_1-t_2
#plt.plot(lat,deltaT_altsurf,label='alt surface albedo')
print('Arctic T altsurf',np.average(deltaT_altsurf,weights=weightarr_Arctic))


# clear sky albedo

t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_plio * sigma)
t_1=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-(plio_A-A_clr_diff))) + H_plio) / (emis_plio * sigma)
t_2=t4 ** (1./4.)

deltaT_csalbedo=t_1-t_2
ax.plot(lat,deltaT_csalbedo,label='clear sky albedo',color="chocolate",linestyle="dashed")
print('Arctic T csalbedo',np.average(deltaT_csalbedo,weights=weightarr_Arctic))

#cloud albedo

t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_plio * sigma)
t_1=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-(plio_A-A_cld_diff))) + H_plio) / (emis_plio * sigma)
t_2=t4 ** (1./4.)

deltaT_cldalbedo=t_1-t_2
ax.plot(lat,deltaT_cldalbedo,label='cloud albedo',color="purple")
print('Arctic T cloud albedo',np.average(deltaT_cldalbedo,weights=weightarr_Arctic))



t4=((sw_down_plio * (1.0-plio_A)) + H_plio) / (emis_plio * sigma)
t_1=t4 ** (1./4.)
t4=((sw_down_plio * (1.0-plio_A)) + H_pi) / (emis_plio * sigma)
t_2=t4 ** (1./4.)

deltaT_H=t_1-t_2
ax.plot(lat,deltaT_H,label='heat transport',color="red")
print('Arctic T Heat transport',np.average(deltaT_H,weights=weightarr_Arctic))

deltaT=temp_plio-temp_pi
#plt.plot(lat,deltaT,label='actual temperature change')
print('Arctic T total change',np.average(deltaT,weights=weightarr_Arctic))


total_components=deltaT_gge+deltaT_ce+deltaT_surfalpha+deltaT_csalbedo+deltaT_cldalbedo+deltaT_H
#plt.plot(lat,total_components,label='total accountable')


deltaT_cloud=deltaT_ce+deltaT_cldalbedo

print('Mean T gge',np.average(deltaT_gge,weights=weightarr))
print('Mean T topo',np.average(delta_T_topo,weights=weightarr))
print('Mean T cloud emisivity',np.average(deltaT_ce,weights=weightarr))
print('Mean T surface albedo',np.average(deltaT_surfalpha,weights=weightarr))
print('Mean T csalbedo',np.average(deltaT_csalbedo,weights=weightarr))
print('Mean T cloud albedo',np.average(deltaT_cldalbedo,weights=weightarr))
print('Mean T Heat transport',np.average(deltaT_H,weights=weightarr))
print('Mean T total change',np.average(deltaT,weights=weightarr))






print('Mean cloud changes',np.average(deltaT_cloud,weights=weightarr))


plt.ylim(-4.0,8.0)
mp.pyplot.axhline(y=0,xmin=-90,xmax=90,color='black')
if HadCM3 == 'y':
    plt.title('d) Energy Balance: HadCM3',loc='left',fontsize=18)
else:
    plt.title('c) Energy Balance: HadGEM2',loc='left',fontsize=18)
plt.xlabel('Latitude',fontsize=15)
degC=u'\N{DEGREE SIGN}'+'C'
plt.ylabel('Warming ('+degC+')',fontsize=15)


plt.legend(fontsize=15)
plt.tick_params(axis='both',labelsize=15)
box=ax.get_position()
ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
ax.legend(loc='lower center',bbox_to_anchor=(0.5,-0.4),ncol=3)


fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_energybal'+anom_expt+'-'+cntrl_expt+'.eps' 
plt.savefig(fileout, bbox_inches='tight')  
fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_energybal'+anom_expt+'-'+cntrl_expt+'.png' 
plt.savefig(fileout, bbox_inches='tight')  

plt.close()

filetext='/home/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_energybal'+anom_expt+'-'+cntrl_expt+'.tex'
f= open(filetext,"w+")
f.write('latitude,Greenhouse gas emissivity,Cloud emissivity,surfacealbedo,'+
        'clear sky albedo,cloud albedo,heat_transport \n')
for i in range(0,len(lat)):
    f.write(np.str(lat[i])+','+np.str(deltaT_gge[i])+','+
            np.str(deltaT_ce[i])+','+
            np.str(deltaT_surfalpha[i])+','+
            np.str(deltaT_csalbedo[i])+','+
            np.str(deltaT_cldalbedo[i])+','+np.str(deltaT_H[i])
            +'\n')
f.close()


####

