#!/usr/bin/env python2.7
#NAME
#    PLOT_ENERGYBAL
#PURPOSE
#    This program will plot the energy balance for the pliocene simulations
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
import iris
import iris.plot as iplt
import iris.quickplot as qplt
from mpl_toolkits.basemap import Basemap, shiftgrid


#functions are:
#  def global_enbal
#  def seasmean


def global_enbal(exptid):
    """
    this looks to be a check that the temperature contribution from all the different sources is equal to the modelled temperature
    """

    def mean_data(filename):
        """
        gets the mean data from the file"
        """
        
        f = open(filename,"r")
        lines = f.readlines()
        avgs = lines[2]
        mean, sd = avgs.split(',')
       
        return mean

    #===================================================================
    # read in data from  average temperature files produced for Dan

    T_mean = np.float(mean_data(FILESTART + '/' + exptid + '.NearSurfaceTemperature.data.txt'))

    # get upward and downward sw radiation at the top of the atmosphere
    # remember alpha = sw_up_toa / sw_down_toa (rsut / rsdt)
    # incoming sw

    rsut_mean = np.float(mean_data(FILESTART + '/' + exptid + '.rsut.data.txt'))
    rsdt_mean = np.float(mean_data(FILESTART + '/' + exptid + '.rsdt.data.txt'))
    rsut_cube = iris.load_cube(FILESTART + '/' + exptid + '.rsut.allmean.nc')
    rsdt_cube = iris.load_cube(FILESTART + '/' + exptid + '.rsdt.allmean.nc')


    # get terms for effctive longwave emissivity e = lw_up_toa / lw_up_surf
    # = rlut / rlus
    
    rlut_mean = np.float(mean_data(FILESTART + '/' + exptid + '.rlut.data.txt'))
    rlut_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlut.allmean.nc')

    if MODELNAME == 'HadCM3' or (MODELNAME == 'CESM2'):
        # lw upward surface = lw_down_surf - netdown_surf I think.  
        rlds_mean = np.float(mean_data(FILESTART + '/' + exptid + '.rlds.data.txt'))
        rlds_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlds.allmean.nc')
        flns_mean = np.float(mean_data(FILESTART + '/' + exptid + '.flns.data.txt'))
        flns_cube = iris.load_cube(FILESTART + '/' + exptid + '.flns.allmean.nc')
        rlus_mean = rlds_mean - flns_mean
        rlus_cube = rlds_cube - flns_cube
    else:
        rlus_mean = np.float(mean_data(FILESTART + '/' + exptid + '.rlus.data.txt'))
        rlus_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlus.allmean.nc')

    # ====================================
    # get alpha and epsilon

    mean_alpha=rsut_mean/rsdt_mean
    alpha_cube = rsut_cube / rsdt_cube

    mean_epsilon=rlut_mean/rlus_mean
    epsilon_cube = rlut_cube / rlus_cube

  
    #============================================
    # calculate terms in equation
    So=1367 # solar constant
    sigma=5.67 * (10.0 ** (-8.))

    print('alphas',mean_alpha)
    print('epsilon',mean_epsilon)

    #t4=So / 4.0 * (1.0-alpha) / (epsilon * sigma)
    #t=t4 ** (1./4.)

    # calculate t_mean using average values of globe
    t4_mean=(So / 4.0) * (1.0-mean_alpha) / (mean_epsilon * sigma)
    t_mean_formula=t4_mean ** (1./4.)

    # calculate t_mean at every point of the globe and average afterwards
    t4_data=(So / 4.0) * (1.0-alpha_cube.data) / (epsilon_cube.data * sigma)
    t_mean_data=t4_data ** (1./4.)
    t_mean_cube = alpha_cube.copy(data=t_mean_data)

    t_mean_cube.coord('latitude').guess_bounds()
    t_mean_cube.coord('longitude').guess_bounds()
   
    grid_areas = iris.analysis.cartography.area_weights(t_mean_cube)
    t_mean_avg_cube = t_mean_cube.collapsed(['latitude', 'longitude'],
                              iris.analysis.MEAN, weights=grid_areas)
    
    
    print('globvals',So/4.0,mean_alpha,mean_epsilon,sigma)



    print('    ')
    print(exptid)
    print('=========')
    print('t obs=',T_mean,' degC')
    print('t mean formula=',t_mean_formula-273.15,' degC')
    print('t mean from average=',t_mean_avg_cube.data - 273.15)

   
#end def global_enbal

################################################
def dh_zonal_enbal(exptid):

   
    #==============
    # read in data from  average temperature files produced for Dan

    T_cube = iris.load_cube(FILESTART + '/' + exptid + '.NearSurfaceTemperature.allmean.nc')
 
    # get upward and downward sw radiation at top of atmos
    rsut_cube = iris.load_cube(FILESTART + '/' + exptid + '.rsut.allmean.nc')
    rsdt_cube = iris.load_cube(FILESTART + '/' + exptid + '.rsdt.allmean.nc')
  
    # get upward sw radiation at toa clear sky (downwards same as rsdt)
    if (MODELNAME == 'NorESM1-F' or MODELNAME == 'NorESM-L' or 
        (MODELNAME == 'CESM2' and exptid == 'E400')):
        rsntcs_cube = iris.load_cube(FILESTART + '/' + exptid + '.rsntcs.allmean.nc')
        rsnt_cube = iris.load_cube(FILESTART + '/' + exptid + '.rsnt.allmean.nc')
        rsutcs_cube = rsdt_cube - rsntcs_cube
       
    else:
        rsutcs_cube = iris.load_cube(FILESTART + '/' + exptid + '.rsutcs.allmean.nc')
        
  
   # outgoing lw  (toa) rlut and clear sky rlutcs
    rlutcs_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlutcs.allmean.nc') 
    rlut_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlut.allmean.nc')


    # outgoing lw  (surface) rlus and clear sky rluscs
    if MODELNAME == 'HadCM3':
        # lw upward surface = lw_down_surf - netdown_surf I think.  
        rlds_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlds.allmean.nc')
        rlns_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlns.allmean.nc')
        rlus_cube = rlds_cube - rlns_cube
    else:
        rlus_cube = iris.load_cube(FILESTART + '/' + exptid + '.rlus.allmean.nc')
    if MODELNAME == 'COSMOS':
        rlus_cube = rlus_cube * -1.0
        rsut_cube = rsut_cube * -1.0

    # rluscs = tested using rlus
  
    rluscs_cube = rlus_cube
    
    # ====================================
    # get the zonal average fields


    alpha_cube = rsut_cube / rsdt_cube
    alpha_cube_cs = rsutcs_cube / rsdt_cube
    epsilon_cube = rlut_cube / rlus_cube
    epsilon_cube_cs = rlutcs_cube / rluscs_cube


    sigma=5.67 * (10.0 ** (-8.))

    #============================================
    # calculate terms in equation
    H_cube=(-1.0) * ((rsdt_cube - rsut_cube) -(rlut_cube))
   
    components=[epsilon_cube, epsilon_cube_cs, alpha_cube, alpha_cube_cs, 
                H_cube, T_cube, rsdt_cube, rlutcs_cube, rluscs_cube]
    return(components)

#end def dh_zonal_enbal


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
        orog_fname='/nfs/hera2/apps/metadata/ancil/preind2/qrparm.orog.nc'
    if exptname=='xibot' or exptname=='xoorb' or exptname=='xoorf':
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
def main_dh_zonal(preind_expt,plio_expt):
    """
    performs the energy baland calculation as detailed in dan hills paper
    """
    # get all the fields needed for the energy balance
    components=dh_zonal_enbal(preind_expt)
    emis_pi_cube=components[0]
    emis_pi_cs_cube=components[1]
    alpha_pi_cube=components[2]
    alpha_pi_cs_cube=components[3]
    H_pi_cube=components[4]
    temp_pi_cube=components[5]
    sw_down_pi_cube=components[6]
    rlutcs_pi_cube = components[7]
    rluscs_pi_cube = components[8]
    
    components=dh_zonal_enbal(plio_expt)
    emis_plio_cube = components[0]
    emis_plio_cs_cube = components[1]
    alpha_plio_cube = components[2]
    alpha_plio_cs_cube = components[3]
    H_plio_cube = components[4]
    temp_plio_cube = components[5]
    sw_down_plio_cube = components[6]
    rlutcs_plio_cube = components[7]
    rluscs_plio_cube = components[8]
  
 
    sigma=5.67 * (10.0**-8)

    # do energy balance.  

    # greenhouse gas and topography
    t4 = ((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_plio_cs_cube.data * sigma)
    t_1=t4 ** (1./4.)
    t4 = ((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_pi_cs_cube.data * sigma)
    t_2 = t4 ** (1./4.)

    deltaT_gge_cube = alpha_plio_cube.copy(data = t_1 - t_2)
    deltaT_gge_lat_cube= deltaT_gge_cube.collapsed(['longitude'],
                              iris.analysis.MEAN)

    emis_plio_cs_cube.long_name = 'pliocene cs emissivity'
    emis_pi_cs_cube.long_name = 'pi cs emissivity'
    emis_diff = emis_plio_cs_cube - emis_pi_cs_cube
    emis_diff.long_name = 'cs emissivity difference'
    emis_diff_lat = emis_diff.collapsed('longitude',iris.analysis.MEAN)
    rlutcs_diff = rlutcs_plio_cube - rlutcs_pi_cube
    rlutcs_diff_lat = rlutcs_diff.collapsed('longitude', iris.analysis.MEAN)
    rlutcs_diff_lat.long_name = 'rlutcs diff'
    rluscs_diff = rluscs_plio_cube - rluscs_pi_cube
    rluscs_diff_lat = rluscs_diff.collapsed('longitude', iris.analysis.MEAN)
    rluscs_diff_lat.long_name = 'rluscs diff'
    iris.save([emis_plio_cs_cube, emis_pi_cs_cube, emis_diff, emis_diff_lat, rlutcs_diff_lat, rluscs_diff_lat],MODELNAME + '_emis.nc')
   
    lat = deltaT_gge_lat_cube.coord('latitude').points
    plt.plot(lat,deltaT_gge_lat_cube.data,label='GHG+ topography',color="blue")
  
    # cloud emisivity
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_plio_cs_cube.data * sigma)
    t_2=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_pi_cube.data * sigma)
    t_3=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_pi_cs_cube.data * sigma)
    t_4=t4 ** (1./4.)

    deltaT_ce_cube = alpha_plio_cube.copy(data = (t_1-t_2) - (t_3 - t_4))
    deltaT_ce_lat_cube= deltaT_ce_cube.collapsed(['longitude'],
                              iris.analysis.MEAN)
   
    plt.plot(lat,deltaT_ce_lat_cube.data, 
             label='cloud emissivity',color="orange")

  
    # cloud albedo
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cs_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_2=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_pi_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_3=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_pi_cs_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_4=t4 ** (1./4.)

    deltaT_ca_cube = alpha_plio_cube.copy(data = (t_1-t_2) - (t_3 - t_4))
    deltaT_ca_lat_cube= deltaT_ca_cube.collapsed(['longitude'],
                              iris.analysis.MEAN)
   
    plt.plot(lat,deltaT_ca_lat_cube.data, label='cloud albedo',color="purple")
   
   
    # clear sky albedo
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cs_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_pi_cs_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_2=t4 ** (1./4.)

    print('t1 components',sw_down_plio_cube.data[178,0],alpha_plio_cs_cube.data[178,0],H_plio_cube.data[178,0],emis_plio_cube.data[178,0])
    print('t2 components',sw_down_plio_cube.data[178,0],alpha_pi_cs_cube.data[178,0],H_plio_cube.data[178,0],emis_plio_cube.data[178,0])

    deltaT_csa_cube = alpha_plio_cube.copy(data = t_1-t_2)
    deltaT_csa_lat_cube= deltaT_csa_cube.collapsed(['longitude'],
                              iris.analysis.MEAN)

    for i, latval in enumerate(deltaT_csa_lat_cube.coord('latitude').points):
        if latval == 88.5:
            print('j2',i,latval,deltaT_csa_lat_cube.data[i],t_1[i,0],t_2[i,0],t_1[i,0] - t_2[i,0])
    
  
    plt.plot(lat,deltaT_csa_lat_cube.data,label='clear sky albedo',color="chocolate",linestyle="dashdot")
   
    # heat transport

    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_plio_cube.data) / (emis_plio_cube.data * sigma)
    t_1=t4 ** (1./4.)
    t4=((sw_down_plio_cube.data * (1.0-alpha_plio_cube.data)) + H_pi_cube.data) / (emis_plio_cube.data * sigma)
    t_2=t4 ** (1./4.)

    deltaT_H_cube = alpha_plio_cube.copy(data=t_1-t_2)
    deltaT_H_lat_cube= deltaT_H_cube.collapsed(['longitude'],
                              iris.analysis.MEAN)
  
    plt.plot(lat,deltaT_H_lat_cube.data,label='heat transport',color="red")
   
    # compare with total temperature change
    deltaT_cube=temp_plio_cube-temp_pi_cube
    deltaT_lat_cube= deltaT_cube.collapsed(['longitude'],
                              iris.analysis.MEAN)
  
    plt.plot(lat,deltaT_lat_cube.data,label='actual temperature change')


    # sum of all energy balance terms
    deltaT_sum = (deltaT_gge_lat_cube.data +  deltaT_ce_lat_cube.data + 
                  deltaT_ca_lat_cube.data + deltaT_csa_lat_cube.data + 
                  deltaT_H_lat_cube.data)
    plt.plot(lat,deltaT_sum, color='black',linestyle='dotted',
             label='sum EB terms')
    #print(deltaT_gge_lat_cube.data)
    


    for i,latind in enumerate(lat):
        if latind == 75.5:
            print('found totals',deltaT_sum[i], deltaT_lat_cube.data[i])
            print('comps',deltaT_gge_lat_cube.data[i],deltaT_ce_lat_cube.data[i],
                  deltaT_ca_lat_cube.data[i], deltaT_csa_lat_cube.data[i],
                  deltaT_H_lat_cube.data[i])
  #          sys.exit(0)
        
   

    plt.legend()
    plt.ylim(-8.0,15.0)
    mp.pyplot.axhline(y=0,xmin=-90,xmax=90,color='black')
    plt.title('Energy balance: '+ MODELNAME+ ' - Dan Hill')
    plt.xlabel('latitude')
    plt.ylabel('pliocene warming')
    
    filestart='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/PLIOMIP2/energybal/DH_energybal_'+MODELNAME + '_' + EXPT + '_' + CNTL
    plt.savefig(filestart + '.eps', bbox_inches='tight')  
    plt.savefig(filestart + '.png', bbox_inches='tight')  

    plt.close()
    sys.exit(0)

# end of main part of Dan Hills energy balance

def main_rf_energybal(modelname):
    components=rf_zonal_enbal(preind_expt,HadCM3)
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


    components=rf_zonal_enbal(pliop2_expt,HadCM3)
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

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_albedo'+pliop2_expt+'-'+preind_expt+'.eps' 
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
    plt.ylabel('Pliocene Warming ('+degC+')',fontsize=15)


    plt.legend(fontsize=15)
    plt.tick_params(axis='both',labelsize=15)
    box=ax.get_position()
    ax.set_position([box.x0,box.y0+box.height*0.1,box.width,box.height*0.9])
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,-0.4),ncol=3)


    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_energybal'+pliop2_expt+'-'+preind_expt+'.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_energybal'+pliop2_expt+'-'+preind_expt+'.png' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()

    filetext='/home/earjcti/PYTHON/PLOTS/HadGEM2/plot_energybal/RF_energybal'+pliop2_expt+'-'+preind_expt+'.tex'
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






################################
# main program

# annual mean

MODELNAME = 'NorESM1-F'
FILESTART = '/nfs/hera1/earjcti/regridded/' + MODELNAME
CNTL = 'E280'
EXPT = 'EOI400'

# check to see if Dan Hills global energy balance equation is correct
#global_enbal(CNTL)
#global_enbal(EXPT)
#sys.exit(0)


# Dan Hills energy balance
main_dh_zonal(CNTL,EXPT)


# Ran Fengss energy balance

main_rf_energybal(modelname)
