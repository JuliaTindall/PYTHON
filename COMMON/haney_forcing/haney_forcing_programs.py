# -*- coding: utf-8 -*-
"""
Spyder Editor

Here we are going to store all the classes and programs for doing our 
Permanent El nino experiments.  

We need to create the reference temperature for Haney forcing
and also the weighting (the relaxation to the forcing at each gridbox)

Julia 7/12/2019
"""
import iris
import iris.quickplot as qplt
import numpy as np
import itertools
import sys

class haney_nc():
    """
    this class contains all the stuff needed to make the netcdf file 
    of the reference temperatures which the program will relax to
    """
    def __init__(self):
        """
        sets up the filenames and variables
        (note we will use ple as the mask as this is on a 2d grid,
        with nan over land points)
        """
        
        self.ORIG_HANEY_FILE = FILESTART + 'qrclim.Haney_elnino.nc'
        self.MASK_FILE = FILESTART + 'xogzi.ostart.p50c1.nc'
        self.MASK_VARNAME = 'PLE:PRECIP-EVAP INTO OCEAN KG/M2/S A'
        self.NEW_HANEY_FILE = FILESTART + 'qrclim.Haney_elnino_new.nc'
        
    def create_nc(self):
        """
        creates an haney forcing sst with the same Lsm as is in the maskfile

        Returns
        -------
        None.

        """
        
        # get a lsm where sea =20. land is masked
        lsm_cube = self.get_mask()
        lon_lsm = lsm_cube.coord('longitude').points
        ydim, xdim = lsm_cube.shape
        lsm12_init = (np.broadcast_to(lsm_cube.data,(12, 1, ydim, xdim)))
        lsm12_alt = lsm12_init - 20.
        
        
        # get original cubes
        orig_cubes = iris.load(self.ORIG_HANEY_FILE)
        ncubes=len(orig_cubes)
        
         # if we want to vary the temperature by latitude
        if MASKTYPE == 'latvar':
            lsm12 = self.get_mask_latvar(lsm12_init)
        elif MASKTYPE == 'LANINA':
            lsm12 = self.get_mask_lanina(lsm12_init, lon_lsm)
        else: # set to constant
            lsm12 = lsm12_init
            
        
        
        # overwrite all lsm cubes with data from lsm cubes
        # note there are 12 months
        cubelist = iris.cube.CubeList([])
        for i in range(0, ncubes):
            cube = orig_cubes[i]
            
            if i == 0:
                cubelist.append(cube.copy(data=lsm12))
            else:
                cubelist.append(cube.copy(data=lsm12_alt))
                newcube = cube.copy(data=lsm12_alt)
                print(newcube.data)
        
        outfile = FILESTART + OUTFILEEND
        iris.save(cubelist, outfile, netcdf_format='NETCDF3_CLASSIC',fill_value=1.0E20)
       
        
    def get_mask(self):
        """
        get's a 2d cube of the land sea mask
        land is missing data indicator
        sea is set to 20degC

        Returns
        -------
        lsmcube

        """
        
        cube = iris.load_cube(self.MASK_FILE, self.MASK_VARNAME)
        
        # remove all dimensions except longitude and latitude
        for coord in cube.coords():        
            name=coord.var_name
            if name !='latitude' and name!='longitude':
                cubenew = cube.collapsed(name, iris.analysis.MEAN)
                cube = cubenew
               
              
        mask_data = (cube.data)*0. + 20. # set to 20 everywhere
        mask_data  = np.ma.masked_where(np.ma.getmask(cube.data), 
                                        mask_data)
        mask_data = np.ma.filled(mask_data, 1.0E20)
        
        mask_cube = cube.copy(data = mask_data)
       
        return mask_cube


    def get_mask_latvar(self, lsmin):
        """
         Parameters
        ----------
        lsmin : array x y z t (t=12)
            array containing the lsm with all non masked values set to 20degC
      

        Returns
        -------
        lsmout

        This function will use the lsmin 
        It will overwrite it with the zonally averaged temperatures from paul's 
        Haney forcing file
        """
        
        # get data and find zonal mean
        cube = iris.load_cube(self.ORIG_HANEY_FILE, 'REF. SEA SURF. TEMPERATURE  DEG.C  A')
        haney_zm_cube = cube.collapsed('longitude', iris.analysis.MEAN)
        haney_zm_data = haney_zm_cube.data
       
        
        # overwrite data with zonal means
        nt, nz, ny, nx = np.shape(lsmin)
        lsmout = np.ma.masked_all_like(lsmin)
        
        for t in range(0, nt):
            for k in range(0, nz):
                for j in range(0, ny):
                    for i in range(0,nx):
                        if lsmin[t, k, j, i] < 100.0:
                            lsmout[t, k, j, i] = haney_zm_data[t, k, j]
                      
            
       
        return lsmout
    
    def get_mask_lanina(self, lsmin, lonlsmin):
        """
        Parameters
        ----------
        lsmin : array x y z t (t=12)
            array containing the lsm with all non masked values set to 20degC
        lonlsmin ; array containing the longitudes of the lsmin

        Returns
        -------
        lsmout

        This function will use the lsmin as a template
        It will get then get the average gradient across the Pacific from each month
        and will multiply it by 1.5 to get a strong La Nina
        
        The lsmin will be overwritten with the average monthly temperature everywhere
        except across the Pacific where it will be overwritten with the strong la nina
        """
        monthnames = ['ja','fb','mr','ar','my','jn','jl','ag','sp','ot','nv','dc']
        
        # setup outfile to overwrite with lanina
        nt, nz, ny, nx = np.shape(lsmin)
        lsmout = np.ma.masked_all_like(lsmin)
        
        for monthno, month in enumerate(monthnames):
            # get average data
            cube = iris.load_cube('AVG_FIELDS/xoorb-sst-mean-' + month + '.nc', 'TEMPERATURE (OCEAN)  DEG.C')
            nlats = len(cube.coord('latitude').points)
            nlons = len(cube.coord('longitude').points)
            
            # get pacific basin (from the writenc class)
            if monthno == 0:
              obj = haney_txt()
              pacific_region = obj.get_pacific_basin(FILESTART+'/basin_hadcm3', nlats, nlons)
      
            #####################################
            # set up temperature in Pacific region
            # this will be the difference between the Pacific mean and the local mean multiplied by 1.5
        
            lanina = cube.data #currently a copy of the average data
           
            pac_temp = pacific_region * lanina # this is the temperatures in the pacific but zero elsewhere
            pac_temp[pac_temp == 0.] = np.nan #this is temperature in pacific but nan elsewhere
            for j in range(0,nlats):
                total_pac_lat = np.sum(pacific_region[j, :])
                if total_pac_lat > 0.:
                    avg_pac_lat = np.nanmean(pac_temp[j, :]) # average temperature at that latitude in pacific
            
                # if we are in Pacific replace temperature with 
                # 1.5 * the local difference from the zonalmean
                for i in range(0, nlons):
                    if pacific_region[j,i] == 1.0:
                        lanina[j, i] = avg_pac_lat + (1.5 * (lanina[j,i] - avg_pac_lat))
                
        
            
        
           
       
            lanina_lons = cube.coord('longitude').points
        
            for i in range(0,nx):
                # note the input array has more longitudes than the mean array (use modulus)
                i2 = np.where(lanina_lons == np.mod(lonlsmin[i],360.))
          
                for j in range(0, ny):
                
                    if lsmin[monthno, :, j, i] < 100.0:
                        lsmout[monthno, :, j, i] = lanina[j, i2]
                        
                      
            print('monthno is',monthno,lsmout[monthno,0,20,20])
            
       
        return lsmout
# end of class
  


class haney_txt():
    """
    this class contains all the stuff needed to make the txt file which 
    contains the strength of the relaxation for each gridbox
    """
    def __init__(self):
        """
        sets up the filenames and variables
        """
        
        self.ORIG_HANEY_NC = FILESTART + 'qrclim.Haney_elnino.nc'
        self.HANEY_FILE_TXT = FILESTART + 'haney_elnino.dat'
        self.MAX_NHANEY = 163.  # this is the haney forcing on the equator in W/m2
        
    def create_txt(self):
        """
        creates the text file and writes out

        Returns
        -------
        None.

        """
        #1. get latitudes
        cube = iris.load_cube(self.ORIG_HANEY_NC, 'REF. SEA SURF. TEMPERATURE  DEG.C  A')
        lats = cube.coord('latitude').points
        nlats = len(lats)
        nlons = len(cube.coord('longitude').points)
      
        # Haney restoring function 14 day timescale 40S-40N with
        # value set to that of the Western Pacific warm pool.  ie expand western Pacific warm pool across 
        # Pacific      
        #The coupling coefficient κ has units of watts per square meter per kelvin. 
        #A timescale can be derived as τ = ρcph/κ, where ρ is the density, 
        #cp is the specific heat of seawater, and h is a typical mixed layer depth. 
        # Weaver and Sarachik (1991) used a timescale of 25 days.
        
        # if timescale is 14 days what is k
        # k= density * cp * mixed layer depth / 14
        #  = 1022 (kg/m3) * 4011 (J / kg/ k) * 50 (m) / 1209600 secs
        #  = 170 (J)  / (m2 * K * seconds) or 170 W/m2/K
        # so assuming that the mixed layer depth is 50m the
        # restoring timescale of 14 days correspoinds to haney forcing
        # coefficient of 170W/m2/s
        
        outarr = np.zeros((nlats,nlons))
      
        for j, lat in enumerate(lats):
            if np.abs(lat) < 20.0:
                outarr[j, :] = 170.0
            elif np.abs(lat) < 40.0:
                #lat=40 forcing 0, lat=20 forcing 170
                # forcing changes linearly between
                outarr[j, :] = 170. * (1.0- ((np.abs(lat) - 20.)/20.))
            else:
                outarr[j, :] = 0.0
                
        # if mask region is PACIFIC set to zero in other regions
        if MASKREGION == 'Pacific':
            # note that this was adjusted by me so that the Indian ocean was not included
            pac_basin = self.get_pacific_basin(FILESTART+ '/basin_hadcm3', nlats, nlons)
            outarr = outarr * pac_basin # zeros away from pacific
          
        # check the mask
       
        #cubetoplot = cube.collapsed(['t','unspecified'], iris.analysis.MEAN)
        #cubeplot2 = cubetoplot.copy(data=cubetoplot.data * pac_basin)
        #qplt.contourf(cubeplot2)
        
            
        f=open(self.HANEY_FILE_TXT,'w+')
        #np.savetxt(f, np.transpose(lats[:]), fmt='%9.2f')
        juheader = ('Each column is a latitude from ' + np.str(lats[0])
                   + ' to ' + np.str(lats[-1]))
        np.savetxt(f, np.transpose(outarr), header=juheader, fmt='%9.1f') 
        f.close()

    def get_pacific_basin(self,filename, nlats, nlons):
        """

        Parameters
        ----------
        filename : the name of the basin indices file

        Returns
        -------
        a numpy array containing the basin indices (1-pacific, 0 = not pacific)

        """
        
        f1 = open(filename, "r")
      
        lines = f1.readlines()
        
        # lines 1-2 are title, line 3-146 are data, line 147 onwards metadata
        rowno = np.zeros(144, dtype=int)
        pacstart = np.zeros(144, dtype=int)
        pacend = np.zeros(144, dtype=int)
        for count, line in enumerate(lines):
            if count in range(3,147):
               linedata = line.split() # index0 is row
                                       # index5 is pacstart
                                       # index6 is pacend
               rowno[count-3] = np.int(linedata[0])
               pacstart[count-3] = np.int(linedata[5])
               pacend[count-3] = np.int(linedata[6])
               
        print(pacstart)
        
        region = np.zeros((nlats, nlons))
        
        for count, row in enumerate(rowno):
            colstart = pacstart[count]
            colend = pacend[count]
           
            for i in range(colstart, colend+1):
                if (row >= 40) and (row <= 104):
                    region[row-1, i] = 1.0
       
        return region
     

LINUX_WIN = 'l'
FILESTART_DICT = {
                   "w" : "C:\\Users\\julia\\OneDrive\\WORK\\DATA\\TEMPORARY\\",
                   "l" : "/nfs/hera1/earjcti/Xiaofang/HANEY/"
                 }

OUTFILE_DICT = {
                   "LANINA" : "Hadley_lanina_new.nc",
                   "latvar" : "Hadley_latconst_new.nc",
                   "fixed" : "Hadley_constant.nc"
                   
                 }

MASKTYPE = 'LANINA' # fixed; fixed at 20deg, latvar: varies with latitude only
                    # La Nina : we will increase the gradient across the Pacific (like a lanina)
MASKREGION = 'Pacific' # Pacific, None

FILESTART = FILESTART_DICT.get(LINUX_WIN)
OUTFILEEND = OUTFILE_DICT.get(MASKTYPE)


#set up and run the code to produce the Haney nc file
#obj = haney_nc()
#obj.create_nc()

#set up and run the code to produce the Haney forcing text file
#this file contains the relaxiation coefficient for each gridcell
obj = haney_txt()
obj.create_txt()
