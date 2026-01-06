#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
#NAME
#    orbit_forcing
#PURPOSE
#    This program will assess the orbital forcing between two different
# experiments.  We are actually trying to replicate the difference in incoming
# shortwave radiation
#
#  It actually also does the calendar correction for some fields.
#
#
# search for 'main program' to find end of functions
# Julia January 2019
# ammended in September 2019 to add new fields including d18o
"""


import numpy as np
import iris
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import sys
import cf_units

#import basemap
#from mpl_toolkits.basemap import Basemap, shiftgrid


def get_incom_sw_read(filename):
    """
    called by get_incom_sw.  Gets the incoming sw radiation from 'filename'
    """
    f = Dataset(filename)
    lat = f.variables['latitude'][:]
    days = f.variables['t'][:]
    atemp = f.variables['field200'][:] # get incoming sw radiation
    atemp = np.squeeze(atemp)
    f.close()

    return lat, days, atemp

def get_incom_sw_calcshift(lat, d1, d2):
    """
    called by get_incom_sw to calculate the shift required
    d1 and d2 are the data arrays containing the incoming
              shortwave radiation for the experiments
    """

 # get maximum from first experiment
    max_data1_np = np.argmax(d1[0, :]) # find max at northpole
    max_data1_sp = np.argmax(d1[len(lat)-1, :]) # find max at southpole

    print('np', max_data1_np, d1[0, max_data1_np])
    print('sp', max_data1_sp, d1[0, max_data1_sp])


    # get maximum from second experiment
    max_data2_np = np.argmax(d2[0, :]) # find max at northpole
    max_data2_sp = np.argmax(d2[len(lat)-1, :]) # find max at southpole

    print('np', max_data2_np, d2[0, max_data2_np])
    print('sp', max_data2_sp, d2[0, max_data2_sp])

    nh_shift = max_data2_np - max_data1_np
    sh_shift = max_data2_sp - max_data1_sp

    # if the shift calculated from the NH and the SH are within two days
    # use the mean

    if np.abs(nh_shift-sh_shift) > 300:
        if sh_shift < 0:
            sh_shift = sh_shift + 360.
        if nh_shift < 0:
            nh_shift = nh_shift + 360.

    if np.abs(nh_shift - sh_shift) <= 5:
        shiftreq = np.int(np.fix((nh_shift + sh_shift)/2.0))
    else:
        print('we are shifting differently in NH and SH')
        print('nh shift', max_data2_np - max_data1_np)
        print('sh shift', max_data2_sp - max_data1_sp)
        sys.exit(0)
        #shiftreq = nh_shift

    print('nh shift', nh_shift)
    print('sh shift', sh_shift)
    print('shift', shiftreq)
    #sys.exit(0)
    return shiftreq

def get_incom_sw_plotdata(lat_, timesteps_, d1,
                          d2_raw, d2_shifted, expt1_, expt2_, shiftreq_, daily_):

    """
    called by get_incom_sw
    plots the incoming short wave radition and shows
    what the anomaly is like after the data has been shifted in order to
    check that the shift is correct
    """

    anom = d2_shifted - d1

    V = np.arange(0, 600, 50)

    plt.subplot(2, 2, 1)


    plt.contourf(np.arange(0, timesteps_), lat_, d1, V, extend='max')
    plt.title('PI')

    plt.colorbar()

    plt.subplot(2, 2, 2)
    plt.contourf(np.arange(0, timesteps_), lat_, d2_raw, V, extend='max')
    plt.title('new orbit ' + expt2_ + ' - raw')
    plt.colorbar()

    plt.subplot(2, 2, 3)
    plt.contourf(np.arange(0, timesteps_), lat_, d2_shifted, V, extend='max')
    plt.title('new orb (month+' + np.str(shiftreq_) + ')')

    plt.colorbar()

    plt.subplot(2, 2, 4)
    V = np.arange(-100, 105, 5)
    plt.contourf(np.arange(0, timesteps_), lat_, anom, V, cmap='RdBu_r', extend='both')
    plt.title('new orbit - PI')
    if daily_ == 'y':
        plt.plot([151, 151], [-90, 90])
        plt.plot([270, 270], [-90, 90])
    else:
        plt.plot([5, 5], [-90, 90])
        plt.plot([9, 9], [-90, 90])

    plt.colorbar()

    plt.savefig(plotout, bbox_inches='tight')
    plt.close()


def get_incom_sw(filestart1, filestart2, expt1, expt2, daily):
    """
    the top level for getting the incoming shortwave radiation from
    two experiments and seeing how much we have to shift by to get the
    solstaces to line up

    inputs:
    filestart1/2, expt1/2 are the experiment name and the file path
                             for the control/anomaly experiment
    daily is a y/n indicator saying whether the input data is daily data
          (ie from a pa file) or monthly data (ie from a pd file)

    returns:
        shiftreq:  the number of days/months to shift the second experiment
                   by to get the two experiments to line up
    """

    print('j2')
    for i in range(0, len(monthnames)):
        # read in data from expt1 files
        filename1 = filestart1 + monthnames[i] + '.nc'
        lat1, days1, atemp1 = get_incom_sw_read(filename1)

        if i == 0:
            data1 = np.zeros((len(lat1), len(monthnames)*len(days1)))
            data2 = np.zeros((len(lat1), len(monthnames)*len(days1)))
            timesteps = len(monthnames)*len(days1)

        daystart = i*len(days1)
        dayend = (i + 1)*len(days1)


        # read in data from expt2 files
        filename2 = filestart2 + monthnames[i] + '.nc'
        lat2, days2, atemp2 = get_incom_sw_read(filename2)


        if daily == 'y':
            data1[:, daystart:dayend] = np.swapaxes(atemp1[0:len(days1), :, 0], 0, 1)
            data2[:, daystart:dayend] = np.swapaxes(atemp2[0:len(days1), :, 0], 0, 1)
        else:
            data1[:, i] = atemp1[:, 0]
            data2[:, i] = atemp2[:, 0]

    print('j3')
    if daily == 'y':
        shiftreq = get_incom_sw_calcshift(lat1, data1, data2)
    else:
        shiftreq=0.0
    print('j4')

   # shift data

    data2_raw = data2
    data2_shifted = np.zeros(np.shape(data2))

    # shift data2 so that the maximum value in the northern hemisphere is in
    # the same place as it is in data1

    shiftreq = shiftreq * (-1)

    for j in range(0, len(lat1)):
        # shift values using np.roll
        data2_shifted[j, :] = np.roll(data2_raw[j, :], shiftreq, axis=0)




    get_incom_sw_plotdata(lat1, timesteps, data1, data2_raw,
                          data2_shifted, expt1, expt2, shiftreq, daily)



    return shiftreq


#end def get_incom_sw

def shiftexpt_init(mon_, year_, filestart__):

    daystart_ = mon_ * 30
    dayend_ = ((mon_ + 1) * 30)

    if year_ < 100:
        extra_ = 'o'
        yearuse_ = year_
    else:
        extra_ = '80'
        yearuse_ = year_ - 100

    if yearuse_ < 10:
        yearchar_ = '0' + np.str(yearuse_)
    else:
        yearchar_ = np.str(yearuse_)

    filename_ = (filestart__ + extra_ + yearchar_ + monthnames[mon_] + '.nc')

    return daystart_, dayend_, filename_

def shiftexpt_monchunk(cubeorig, mon_, dataarr_sh_, fieldname_):

    daystart = mon_ * 30
    dayend = ((mon_ + 1) * 30)

    # take a copy of the original cube but put in the new data
    newcube = cubeorig.copy(data=dataarr_sh_[daystart:dayend, :, :, :])

    plt.plot(newcube.data[:, 0, 20, 20])


    # rename if required

    if fieldname_ == 'Stash code = 8244':
        newcube.rename('SOIL MOISTURE CONTENT 18O')
    if fieldname_ == 'Stash code = 8246':
        newcube.rename('SOIL MOISTURE CONTENT IN LAYER 18O')

    return newcube

def shiftexpt_d18o(cube):

    cube.coord('level-1').rename('level')
    precip16o = (cube.extract(iris.Constraint(level=1.0)) +
                 cube.extract(iris.Constraint(level=4.0)) +
                 cube.extract(iris.Constraint(level=7.0)) +
                 cube.extract(iris.Constraint(level=10.0)))
    precip18o = (cube.extract(iris.Constraint(level=2.0)) +
                 cube.extract(iris.Constraint(level=5.0)) +
                 cube.extract(iris.Constraint(level=8.0)) +
                 cube.extract(iris.Constraint(level=11.0)))

    cubed18o = ((precip18o / precip16o) - 2005.2E-6)/2005.2E-9
    tempdata = cubed18o.data
    tempdata2 = np.where(np.isfinite(tempdata), tempdata, 1.0E20)
    cubed18o.data = tempdata2
    cubed18o.rename('d18o')
    precip18o.rename('precip 18o')
    precip16o.rename('precip 16o')
    #u = cf_units.Unit('permille')
    cubed18o.units = None
    precip18o.units = 'kg/m2/s'
    precip16o.units = 'kg/m2/s'

    return cubed18o, precip18o, precip16o



#######################################
def shiftexpt(shiftreq, filestart_, fieldname):


# read in 30 years of data for each month

    for mon in range(0, len(monthnames)):

        for year in range(startyear, endyear):
            # get initial data
            daystart, dayend, filename = shiftexpt_init(mon, year, filestart_)

            # read in the initial data (30 days worth into a cube)
            cubes = iris.load(filename, fieldname)
            cube = cubes[0]



            # put cube data from each month into a file for the
            # average year data
            if (year == startyear) & (mon == 0):
                # set up empty cube from
                days, ndepths, nlats, nlons = np.shape(cube)
                ndays = days * 12

                # create a masked array place to put data
                dataarr1d = np.zeros((ndepths, nlats, nlons))
                dataarr1d = np.ma.masked_where(np.ma.getmask(cube.data[0, :, :, :]), dataarr1d)
                dataarr = np.expand_dims(dataarr1d, axis=0)
                dataarr = np.repeat(dataarr, ndays, axis=0)

            dataarr[daystart:dayend, :, :, :] = (
                (dataarr[daystart:dayend, :, :, :] + cube.data))


    # data should all be loaded in now
    # shift it using np.roll


    dataarr = dataarr / (endyear - startyear)

    dataarr_sh = np.roll(dataarr, shiftreq, axis=0)

    plt.subplot(1, 2, 2)
    plt.plot(dataarr[:, 0, 20, 20])
    print(np.shape(dataarr_sh))
    plt.plot(dataarr_sh[:, 0, 20, 20])
    plt.subplot(2, 2, 2)


    # looks like we have put data in correct place so try and split into
    # monthly chunks

    cubemon = iris.cube.CubeList([])
    precip16omon = iris.cube.CubeList([])
    precip18omon = iris.cube.CubeList([])
    for mon in range(0, len(monthnames)):
        chunkcube = shiftexpt_monchunk(cube, mon, dataarr_sh, fieldname)

        if fieldname == 'Stash code = 338': # need to extract d18o
            d18ocube, precip18ocube, precip16ocube = shiftexpt_d18o(chunkcube)
            cubemon.append(d18ocube)
            precip18omon.append(precip18ocube)
            precip16omon.append(precip16ocube)
        else:
            cubemon.append(chunkcube)


    return [cubemon, precip18omon, precip16omon]


#enddef shiftexpt




################################
# main program

monthnames = ['ja', 'fb', 'mr', 'ar', 'my', 'jn', 'jl', 'ag', 'sp', 'ot', 'nv', 'dc']

# timeslices are xiboi=preindustrial, xibol=3205 - km5c',
#xoekc=3205-km5c, xoekd=3060 (K1), xoeke=2950 (G17), xoekf=3155 (KM3)
# others are xogzb, xogzc, xogzd, xogze, xogzf


#############################
# from monthly data         #
#############################
expt1='xogzb'
filestart1='/nfs/hera1/earjcti/um/xogzb/netcdf/xogzba@pdt39'

expt2='xogzc'
# doesnt matter what year we use
filestart2='/nfs/hera1/earjcti/um/xogzc/netcdf/xogzca@pdt39'
#filestart2='/nfs/hera1/earjcti/um/netcdf/xibol_netcdf/xibola@pdy99'
#filestart2='/nfs/hera1/earjcti/Xiaofang/xhgfk_netcdf/pdfiles/xhgfka@pdy99'
plotout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/orbit_forcing/'
               + expt2 + '-' + expt1 + 'monthorig.eps')


print('j1')
shiftreq=get_incom_sw(filestart1,filestart2,expt1,expt2,'n')
print('new shift required',shiftreq)
sys.exit(0)


######################################################
# from daily data                                    #
######################################################
linuxwin = 'l'
if linuxwin == 'w':
    filestart = 'C:\\Users\\julia\\OneDrive\\WORK\\DATA\\'
else:
    filestart = '/nfs/hera1/earjcti/um/'
##############################################
# get the shift required from the daily data
e1 = 'xozza'  # preindustrial = xogzs, xozza
f1 = filestart + 'xozza/netcdf/xozzaa@pam34'

e2 = 'xozzd'
f2 = filestart + e2 + '/netcdf/' + e2 + 'a@pao87'
filestart2='\\Users\\julia\\OneDrive\\DATA\\HadCM3_DATA\\'

if linuxwin == 'l':
    plotout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/orbit_forcing/'
               + e2 + '-' + e1 + 'dailyorig.eps')
else:
    plotout = filestart + e2 + '-' + e1 + 'dailyorig.png'

toshift = get_incom_sw(f1, f2, e1, e2, 'y')
print('old shift required', toshift)
 

########################################
#shift by appropriate amount
#fieldnames = ['Stash code = 338',
#              'INCOMING SW RAD FLUX (TOA): ALL TSS',
#              'TOTAL PRECIPITATION RATE     KG/M2/S',
#              'SOIL MOISTURE CONTENT',
#              'SOIL MOISTURE CONTENT IN A LAYER',
#              'Stash code = 8244',
#              'Stash code = 8246',
#              'U COMPNT OF WIND ON PRESSURE LEVELS',
#              'V COMPNT OF WIND ON PRESSURE LEVELS',
#              'GEOPOTENTIAL HEIGHT: PRESSURE LEVELS',
#              'SURFACE TEMPERATURE AFTER TIMESTEP']

fieldnames = ['INCOMING SW RAD FLUX (TOA): ALL TSS',
               'SURFACE TEMPERATURE AFTER TIMESTEP']



startyear = 51
#endyear=53
endyear = 99


cube_list_ja = iris.cube.CubeList([])
cube_list_fb = iris.cube.CubeList([])
cube_list_mr = iris.cube.CubeList([])
cube_list_ar = iris.cube.CubeList([])
cube_list_my = iris.cube.CubeList([])
cube_list_jn = iris.cube.CubeList([])
cube_list_jl = iris.cube.CubeList([])
cube_list_ag = iris.cube.CubeList([])
cube_list_sp = iris.cube.CubeList([])
cube_list_ot = iris.cube.CubeList([])
cube_list_nv = iris.cube.CubeList([])
cube_list_dc = iris.cube.CubeList([])


#filestarto = (filestart + e2 + '/netcdf/' + e2 + 'o@pf')
#cube1, dummy1, dummy2 = shiftexpt(toshift, filestarto,
#                                  'OCN TOP-LEVEL TEMPERATURE          K')
#cube_list_ja.append(cube1[0])
#cube_list_fb.append(cube1[1])
#cube_list_mr.append(cube1[2])
#cube_list_ar.append(cube1[3])
#cube_list_my.append(cube1[4])
#cube_list_jn.append(cube1[5])
#cube_list_jl.append(cube1[6])
#cube_list_ag.append(cube1[7])
#cube_list_sp.append(cube1[8])
#cube_list_ot.append(cube1[9])
#cube_list_nv.append(cube1[10])
#cube_list_dc.append(cube1[11])


f2 = (filestart + e2 + '/netcdf/' + e2 + 'a@pa')
print(f2)

for field in range(0, len(fieldnames)):
    cube1, cube2, cube3 = shiftexpt(toshift, f2, fieldnames[field])


    cube_list_ja.append(cube1[0])
    cube_list_fb.append(cube1[1])
    cube_list_mr.append(cube1[2])
    cube_list_ar.append(cube1[3])
    cube_list_my.append(cube1[4])
    cube_list_jn.append(cube1[5])
    cube_list_jl.append(cube1[6])
    cube_list_ag.append(cube1[7])
    cube_list_sp.append(cube1[8])
    cube_list_ot.append(cube1[9])
    cube_list_nv.append(cube1[10])
    cube_list_dc.append(cube1[11])

    if fieldnames[field] == 'Stash code = 338':
        cube_list_ja.append(cube2[0])
        cube_list_fb.append(cube2[1])
        cube_list_mr.append(cube2[2])
        cube_list_ar.append(cube2[3])
        cube_list_my.append(cube2[4])
        cube_list_jn.append(cube2[5])
        cube_list_jl.append(cube2[6])
        cube_list_ag.append(cube2[7])
        cube_list_sp.append(cube2[8])
        cube_list_ot.append(cube2[9])
        cube_list_nv.append(cube2[10])
        cube_list_dc.append(cube2[11])
        cube_list_ja.append(cube3[0])
        cube_list_fb.append(cube3[1])
        cube_list_mr.append(cube3[2])
        cube_list_ar.append(cube3[3])
        cube_list_my.append(cube3[4])
        cube_list_jn.append(cube3[5])
        cube_list_jl.append(cube3[6])
        cube_list_ag.append(cube3[7])
        cube_list_sp.append(cube3[8])
        cube_list_ot.append(cube3[9])
        cube_list_nv.append(cube3[10])
        cube_list_dc.append(cube3[11])



fileout = (f2 + '_avgja.nc')
print(f2, fileout)
iris.save(cube_list_ja, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)
print('saved first cube')

fileout = (f2 + '_avgfb.nc')
iris.save(cube_list_fb, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgmr.nc')
iris.save(cube_list_mr, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgar.nc')
iris.save(cube_list_ar, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgmy.nc')
iris.save(cube_list_my, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgjn.nc')
iris.save(cube_list_jn, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgjl.nc')
iris.save(cube_list_jl, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgag.nc')
iris.save(cube_list_ag, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgsp.nc')
iris.save(cube_list_sp, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgot.nc')
iris.save(cube_list_ot, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgnv.nc')
iris.save(cube_list_nv, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)

fileout = (f2 + '_avgdc.nc')
iris.save(cube_list_dc, fileout, netcdf_format='NETCDF3_CLASSIC', fill_value=1.0E20)
print('saved all cubes')


# check shift required is now zero with new files
if linuxwin == 'l':
    plotout = ('/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadCM3/orbit_forcing/'
               + e2 + '-' + e1 + 'dailynew.eps')
else:
    plotout = filestart + e2 + '-' + e1 + 'dailynew.png'
f2 = filestart + e2 + '/netcdf/' + e2 + 'a@pa_avg'

newshift = get_incom_sw(f1, f2, e1, e2, 'y')
print('new shift required', newshift)
