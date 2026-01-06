
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 14.07.2022 by Julia

We are looking at ETCCDI Climate change indicies.  This program will deal with indices 23-24 these are dry spell and wet spell duration index

23 CDD: Maximum length of dry spell, maximum number of consecutive days with RR< 1mm.
Let RRij be the daily predipitation amount on day i in preiod j.  Count the largest number of consecutive days where RRij < 1mm

(ETCCDI did not specify a period.  However I will use a period of 1 year.  However if the 31st december is a dry spell then we will just keep going.  We will not let a dry spell last for more than 180 days).  

23 WDD: Maximum length of wet spell, maximum number of consecutive days with RR> 1mm.
Let RRij be the daily predipitation amount on day i in preiod j.  Count the largest number of consecutive days where RRij > 1mm
This is also known as consecutive wet days.


"""
import numpy as np
import iris
import iris.quickplot as qplt
from iris.coords import DimCoord
import matplotlib.pyplot as plt
import sys
from os.path import exists


def just_testing():
    """
    just check if the sum of the daily precipitation matches the monthly precip
    """
    monthname = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
    daily_cubes = iris.cube.CubeList([])
    monthly_cubes = iris.cube.CubeList([])
    for month in range(0,12):
        pointsstart =  month *30.
        points=np.arange(pointsstart, pointsstart+30, 1.0)
        day = iris.coords.DimCoord(points,
                             long_name = 't', attributes = None,
                             bounds=None, circular=False)
     
        cubes = iris.load('/nfs/hera1/earjcti/um/'+ EXPTNAME + '/pb/' + EXPTNAME + 'a@pb' + EXTRA + '02'+ monthname[month] + '.nc')
        if EXPTNAME == 'xozzm':
            varreq = 'precip'
        else:
            varreq = 'precip_2'
        for cube in cubes:
            if cube.var_name == varreq:
                try:
                    cube.coord('t_1'.rename('t'))
                except:
                    pass
                cube = cube * 24.* 60. * 60.
                cube.units='mm/day'
                cube.long_name = 'precip from daily data'
                iris.util.demote_dim_coord_to_aux_coord(cube, 't')
                cube.remove_coord('t')
                cube.add_dim_coord(day,0)
                cube = iris.util.squeeze(cube)
                daily_cubes.append(cube)
                
                
        cube = iris.load_cube('/nfs/hera1/earjcti/um/'+EXPTNAME+'/pd/' + EXPTNAME + 'a@pd'+ EXTRA + '02' + monthname[month] + '.nc','TOTAL PRECIPITATION RATE     KG/M2/S')
        cube = cube * 24.* 60. * 60.
        cube.units='mm/day'
        cube = iris.util.squeeze(cube)
        monthly_cubes.append(cube)

    iris.util.unify_time_units(daily_cubes)
    iris.util.equalise_attributes(daily_cubes)
    daily_cube = daily_cubes.concatenate_cube()

    print(monthly_cubes)
    iris.util.unify_time_units(monthly_cubes)
    iris.util.equalise_attributes(monthly_cubes)
    monthly_cube = monthly_cubes.merge_cube()
              
    daily_mean_cube = daily_cube.collapsed('t',iris.analysis.MEAN)
    daily_max_cube = daily_cube.collapsed('t',iris.analysis.MAX)
    daily_mean_cube = iris.util.squeeze(daily_mean_cube)
    daily_max_cube = iris.util.squeeze(daily_max_cube)
    monthly_mean_cube = monthly_cube.collapsed('t',iris.analysis.MEAN)
    monthly_mean_cube = iris.util.squeeze(monthly_mean_cube)
    monthly_mean_cube.long_name = 'precip from monthly data'
   
    fig = plt.figure(figsize=[11.0,8.0])
    plt.subplot(221)
    qplt.contourf(daily_mean_cube,levels=np.arange(0,20,2))
    plt.subplot(222)
    qplt.contourf(monthly_mean_cube,levels=np.arange(0,20,2))
    plt.subplot(223)
    qplt.contourf(monthly_mean_cube - daily_mean_cube,levels=np.arange(-1.0,2.0,0.2))
    plt.title('difference between daily and monthly')
    plt.subplot(224)
    qplt.contourf(daily_max_cube,levels=np.arange(0,110,10),extend='max')
    plt.title('maximum precip from daily data')

    plt.savefig('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag23_24/test_' + EXPTNAME + '.eps')
    sys.exit(0)

def get_ann_precip(year):
    """
    gets the precipitation for a year
    """
    cubelist_precip = iris.cube.CubeList([])
      
    filestart = ('/nfs/hera1/earjcti/um/' + EXPTNAME + '/pb/' + 
                 EXPTNAME +  'a@pb' + EXTRA + np.str(year).zfill(2))
    for month in range(0,12):
        filename = filestart + MONTHNAMES.get(month) + '.nc'
        cubes = iris.load(filename)

        for cube in cubes:
            if cube.var_name == 'precip':
                precipcube = cube
                precipcube = iris.util.squeeze(precipcube)
                maxp = np.mean(precipcube.data)
            if cube.var_name == 'precip_1':
                precip1cube = cube
                precip1cube = iris.util.squeeze(precip1cube)
                maxp1 = np.mean(precip1cube.data)
            if cube.var_name == 'precip_2':
                precip2cube = cube
                precip2cube = iris.util.squeeze(precip2cube)
                maxp2 = np.mean(precip2cube.data)
        
       # it didn't like the time coordinate.  Set a new one
        try:
            tcoord = precip1cube.coord('t') 
            tcoordname = 't'
        except: 
            tcoord = precip1cube.coord('t_1')
            tcoordname = 't_1'

        pointsstart = (year * 360.) + (month *30.)
        points=np.arange(pointsstart, pointsstart+30, 1.0)
        day = iris.coords.DimCoord(points,
                             long_name = 't', 
                             attributes = tcoord.attributes,
                             bounds=None, circular=False)
       
        iris.util.demote_dim_coord_to_aux_coord(precip1cube, tcoordname)
        precip1cube.remove_coord(tcoordname)
        precip1cube.add_dim_coord(day,0)
        
        try:
            iris.util.demote_dim_coord_to_aux_coord(precipcube,'t')
            precipcube.remove_coord('t')
        except:
            iris.util.demote_dim_coord_to_aux_coord(precipcube,'t_1')
            precipcube.remove_coord('t_1')
        precipcube.add_dim_coord(day,0)

        try:
            iris.util.demote_dim_coord_to_aux_coord(precip2cube,'t')
            precip2cube.remove_coord('t')
        except:
            iris.util.demote_dim_coord_to_aux_coord(precip2cube,'t_1')
            precip2cube.remove_coord('t_1')
        precip2cube.add_dim_coord(day,0)
    

        maxps = np.max([maxp, maxp1, maxp2])
        minps = np.min([maxp, maxp1, maxp2])
        found = 'n'
        #print(maxp,maxp1,maxp2)
        #print(maxps,minps,maxp)
        
        if minps < maxp < maxps:
            meanp = maxp
            meanprecipcube = precipcube
            print('using precipcube')
            found = 'y'

        if minps < maxp2 < maxps:
            if found == 'y':
                print('error')
                sys.exit(0)
            else:
                meanp = maxp2
                meanprecipcube = precip2cube
                print('using precipcube2')
                found = 'y'

        if minps < maxp1 < maxps:
            if found == 'y':
                print('error')
                sys.exit(0)
            else:
                meanp = maxp1
                meanprecipcube = precip1cube
                print('using precip1cube')
                found = 'y'

        meanprecipcube = meanprecipcube * 60. * 60. * 24.
        print(month,meanprecipcube)
        meanprecipcube.units = 'mm/day'
        cubelist_precip.append(meanprecipcube)
        meanprecipcube=0
      

    iris.util.equalise_attributes(cubelist_precip)
    iris.util.unify_time_units(cubelist_precip)
    ann_precip_cube = cubelist_precip.concatenate_cube()
    ann_precip_cube = iris.util.squeeze(ann_precip_cube)
    return ann_precip_cube

def get_wet_spell(this_year_precip, next_year_precip):
    """
    this will calculate the wet spell as the number of days when
    the precip was > 1mm day
    """
    nt,ny,nx = np.shape(this_year_precip)
    nwet_array = np.zeros((ny,nx))
    ndry_array = np.zeros((ny,nx))
    for j in range(0,ny): 
        for i in range(0,nx):
            nwetdays=0
            count_wet=0
            ndrydays=0
            count_dry=0
            for t in range(0,nt):
                # get wet/dry info
                #print(this_year_precip[t,j,i])
                if this_year_precip[t,j,i] > 1.0: # move to max dry spell if req
                    count_wet = count_wet+1
                    
                    if count_dry > ndrydays:
                        ndrydays=count_dry
                    count_dry=0
                else: # move to maximum wetspell if appropriate and reset
                    nwetdays = np.max([count_wet, nwetdays])
                    count_wet = 0

                    count_dry=count_dry+1
                #print(t,count_wet,count_dry,this_year_precip[t,j,i])
            nwet_array[j,i] = nwetdays
            ndry_array[j,i] = ndrydays
            # if we have got to the end of the year and we are in
            # the middle of the wet spell add some from the next year
            
            if count_dry ==360:
                ndry_array[j,i]=360.
            if count_wet == 360:
                nwet_array[j,i]=360.
            if count_wet > 0 and nwetdays < 360:
                for t in range(0,179):
                    if next_year_precip[t,j,i] > 1.0:
                        count_wet = count_wet+1
                    else:
                        nwetdays = np.max([count_wet,nwetdays])
                        nwet_array[j,i] = np.min([360,nwetdays])
                        break
            if nwet_array[j,i] > 200:
                print('wet > 200 at',j,i,nwet_array[j,i])
            
            # if we have got to the end of the year and we are in
            # the middle of the wet spell add some from the next year
            
            if count_dry > 0 and ndrydays <360:
                for t in range(0,179):
                    if next_year_precip[t,j,i] < 1.0:
                        count_dry = count_dry+1
                    else:
                        ndrydays = np.max([count_dry,ndrydays])
                        ndry_array[j,i] = np.min([360,ndrydays])
                        break
    return nwet_array, ndry_array

def create_cube(datareqd,cube,long_name,year):
    """
    we are going to move the data to a cube.
    but then we are going to add a new time axis for concatenation
    """
    datacube = cube.copy(data=datareqd)
    datacube.long_name = long_name
    datacube.units='days'
    datacube.coord('surface').rename('year')
    datacube = iris.util.new_axis(datacube,scalar_coord='year')
    datacube.coord('year').points = [year]

    return datacube

def calculate_extreme_spell():
    """
    this function will calculate the number of days which belong to an
    extreme spell for each gridbox in the 100 year period
    """
    startyear=1
    endyear=99

    # get mean precipitation data in mm/day 
    # (collapse the time dimension to create a cube to store data)
    this_year_cube = get_ann_precip(startyear)
    lat_lon_cube = this_year_cube.collapsed('t',iris.analysis.MEAN)

    wet_spell_allcubes = iris.cube.CubeList([])
    dry_spell_allcubes = iris.cube.CubeList([])
   
    for year in range(startyear,endyear):
        print(year)
        next_year_cube = get_ann_precip(year+1)
        print('1')
        (wet_spell_data,
         dry_spell_data)= get_wet_spell(this_year_cube.data,next_year_cube.data)
        print('2')

        wet_spell_cube = create_cube(wet_spell_data,lat_lon_cube,'wet spell duration',year)
        dry_spell_cube = create_cube(dry_spell_data,lat_lon_cube,'dry spell duration',year)
        
        wet_spell_allcubes.append(wet_spell_cube)
        dry_spell_allcubes.append(dry_spell_cube)
        
    wet_spell_allyears_cube = wet_spell_allcubes.concatenate_cube()
    dry_spell_allyears_cube = dry_spell_allcubes.concatenate_cube()
    outfile = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
               EXPTNAME + '/diag_23_24/'+ EXTRA+'_wet_spell_duration.nc')
    iris.save([wet_spell_allyears_cube, dry_spell_allyears_cube],
               outfile, netcdf_format="NETCDF3_CLASSIC")
   
def average_and_plot_wet_spell():
    """
    will average the wet spell over all the years in the file and plot
    """

    # read in cubes difference and remove lsm
    filein = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
               EXPTNAME + '/diag_23_24/'+ EXTRA+'_wet_spell_duration.nc')
    wetcube = iris.load_cube(filein,'wet spell duration')
    drycube = iris.load_cube(filein,'dry spell duration')

    wet_avg_expt_cube = wetcube.collapsed('year',iris.analysis.MEAN)
    dry_avg_expt_cube = drycube.collapsed('year',iris.analysis.MEAN)
    

    filein = ('/nfs/hera1/earjcti/PLIOMIP2/LEEDS/HadCM3/ETCCDI/' + 
               CNTLNAME + '/diag_23_24/'+ CNTLEXTRA+'_wet_spell_duration.nc')
    wetcube = iris.load_cube(filein,'wet spell duration')
    drycube = iris.load_cube(filein,'dry spell duration')

    wet_avg_cntl_cube = wetcube.collapsed('year',iris.analysis.MEAN)
    dry_avg_cntl_cube = drycube.collapsed('year',iris.analysis.MEAN)

    wet_anom_cube = wet_avg_expt_cube - wet_avg_cntl_cube
    dry_anom_cube = dry_avg_expt_cube - dry_avg_cntl_cube

    maskcube=iris.load_cube('/nfs/b0164/Data/LEEDS/HadCM3/eoi400/P4_enh_qrparm.mask.nc','LAND MASK (LOGICAL: LAND=TRUE)')
    vals = np.arange(-30,35,5)
    wet_anom_cube.data.mask = (maskcube.data - 1.0) * (-1.0)
    dry_anom_cube.data.mask = (maskcube.data - 1.0) * (-1.0)

    cs = qplt.contourf(wet_anom_cube, levels=vals, cmap='BrBG', extend='both')
    plt.gca().coastlines()
    plt.title('wet spell duration Plio-Pi ')
   

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag23_24/wet_spell_' + EXPTNAME + '-' + CNTLNAME)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

    vals = np.arange(-50,60,10)
   
    cs = qplt.contourf(dry_anom_cube, levels=vals, cmap='BrBG_r', extend='both')
    plt.gca().coastlines()
    plt.title('dry spell duration Plio-Pi ')
   

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag23_24/dry_spell_' + EXPTNAME + '-' + CNTLNAME)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

    # plot pi values
    wet_avg_cntl_cube.data.mask = (maskcube.data - 1.0) * (-1.0)
    dry_avg_cntl_cube.data.mask = (maskcube.data - 1.0) * (-1.0)

    vals = np.arange(0,110,10)
   
    cs = qplt.contourf(dry_avg_cntl_cube, levels=vals, cmap='copper_r', extend='max')
    plt.gca().coastlines()
    plt.title('dry spell duration PI ')
   

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag23_24/dry_spell_' + CNTLNAME)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

    cs = qplt.contourf(wet_avg_cntl_cube, levels=vals, cmap='BuGn', extend='max')
    plt.gca().coastlines()
    plt.title('wet spell duration PI ')
   

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag23_24/wet_spell_' + CNTLNAME)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()


    # plot plio values
    wet_avg_expt_cube.data.mask = (maskcube.data - 1.0) * (-1.0)
    dry_avg_expt_cube.data.mask = (maskcube.data - 1.0) * (-1.0)

    vals = np.arange(0,110,10)
   
    cs = qplt.contourf(dry_avg_expt_cube, levels=vals, cmap='copper_r', extend='max')
    plt.gca().coastlines()
    plt.title('dry spell duration mPWP ')
   

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag23_24/dry_spell_' + EXPTNAME)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

    cs = qplt.contourf(wet_avg_expt_cube, levels=vals, cmap='BuGn', extend='max')
    plt.gca().coastlines()
    plt.title('wet spell duration mPWP ')
   

    fileout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/extremes/diag23_24/wet_spell_' + EXPTNAME)
    plt.savefig(fileout + '.eps')
    plt.savefig(fileout + '.png')
    plt.close()

    
    #wet_leeds = wet_avg_cube.extract(iris.Constraint(longitude=0,latitude=52.5))
    #dry_leeds = dry_avg_cube.extract(iris.Constraint(longitude=0,latitude=52.5))

    #print('wet',wet_leeds.data)
    #print('dry',dry_leeds.data)

    

##########################################################
# main program
MODELNAME = 'HadCM3'  # 'CESM2', 'IPSLCM6A', 'COSMOS', 'EC-Earth3.3', 
                      # 'CESM1.2', 'IPSLCM5A', 'MIROC4m', 'IPSLCM5A2',
                      # 'HadCM3', 'GISS2.1G', 'CCSM4',  'CCSM4-Utr', 
                      # 'CCSM4-UoT','NorESM-L', 'MRI2.3', 'NorESM1-F'
MIN_MAX = 'max' # max is warm spell duration index
                # min is cold spell duration index
print('start of program')
MONTHNAMES = {0:'ja',1:'fb',2:'mr',3:'ar',4:'my',5:'jn',6:'jl',
                  7:'ag',8:'sp',9:'ot',10:'nv',11:'dc'}
  
EXPTNAME = 'tenvj'
EXTRA='o'

# this will check if the average of the daily precipitation matches the
# monthly precipitation
#just_testing()

#nextreme_spell = calculate_extreme_spell()

CNTLNAME = 'xozza'
CNTLEXTRA = 'p'
average_and_plot_wet_spell()
