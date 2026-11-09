#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#Created on June2019 2019

# adapted 2026 for PlioMIP3


#
# This very simple program will convert data that is on a tripolar grid onto a 
# rectilinear grid.  It will still need to be passed through regrid_ocn in order
# to calculate means etc.


import numpy as np
import cf
import iris
#import cfplot as cfp
#import matplotlib.pyplot as plt
from netCDF4 import Dataset, MFDataset
#import sys
#import os


def reformat_with_iris(exptnamein):
    """
    we are going to load the data as an iris cube and add some
    auxillary coordinates and then write out to a file called temporary.nc
    
    the data is currently in filename
    """
    
    print('reformatting ',filename)
    
    origcube = iris.load_cube(filename, 'Near-Surface Air Temperature')
    latcube = iris.load_cube(filename, 'array of t-grid latitudes')
    loncube = iris.load_cube(filename, 'array of t-grid longitudes')
    

    print(origcube)
   
    # promote the auxillary coordinates to dimension coordinates
    nt, ny, nx = origcube.shape
    origcube.coord('nlat').points=np.arange(0,ny,1)
    origcube.coord('nlat').rename('y')
    origcube.coord('y').var_name='y'
    origcube.coord('y').long_name=None
    origcube.coord('y').units=None
    origcube.coord('nlon').points=np.arange(0,nx,1)
    origcube.coord('nlon').rename('x')
    origcube.coord('x').var_name='x'
    origcube.coord('x').long_name=None
    origcube.coord('x').units=None
    print('j3',origcube)
    iris.util.promote_aux_coord_to_dim_coord(origcube, 'y')
    iris.util.promote_aux_coord_to_dim_coord(origcube, 'x')
    print(origcube.coord('x'))
    
 
    
    
    # add an auxillary coordinate for latitude and longitude these are 
    # 2d coordinates
    loncoord=iris.coords.AuxCoord(loncube.data,standard_name='longitude', 
                                  long_name='Longitude',var_name='nav_lon',
                                  units='degrees_east')
    latcoord=iris.coords.AuxCoord(latcube.data,standard_name='latitude', 
                                  long_name='Latitude',var_name='nav_lat',
                                  units='degrees_north')


    origcube.add_aux_coord(loncoord,[1,2])
    origcube.add_aux_coord(latcoord,[1,2])
    
    print('j6')
    print(origcube)
    print(origcube.coord('latitude'))
    #sys.exit(0)
    
    
    iris.save(origcube, exptnamein + '_temporary.nc', 
              fill_value=2.0E20)
    print('saved')

    


#####################################
def regrid_data(fieldnamein,fieldnameout,exptnamein,exptnameout,filename,modelname,linux_win,fielduse,filenameout):

   
    
    print('moodelname is',modelname)
    print('filename is',filename)
    print('fielduse is',fielduse)
      
    
    if ((modelname=='IPSLCM5A') and exptnamein=='Eoi400'):
       # there is a bit of an error in the file calendar so we will 
       # copy the data to a new file but without the error
       with Dataset(filename) as src, Dataset("temporary.nc", "w",format='NETCDF3_CLASSIC') as dst:
        # copy attributes
            for name in src.ncattrs():
                dst.setncattr(name, src.getncattr(name))
                #print('att',name)   
                # copy dimensions
            for name, dimension in src.dimensions.iteritems():
              
                if name != 'tbnds':   # don't copy across time counter bounds
                   # if fielduse=='SeasurfaceTemp_sst' and name=='y':
                   #     name='latitude'
                   # if fielduse=='SeasurfaceTemp_sst' and name=='x':
                   #     name='longitude'
                    dst.createDimension(name, (len(dimension)))
           
                     
            # copy all file data 
            for name, variable in src.variables.iteritems():
                print('name is',name,variable)
                if name !='time_counter_bnds' and name!='time_centered':
                    x = dst.createVariable(name, variable.datatype, 
                                               variable.dimensions)
                       
                    if name=='time_counter':
                    # convert from seconds to days and start at middle of month
                        dst.variables[name][:] = (src.variables[name][:] / (60.*60.*24))-(src.variables[name][0] / (60.*60.*24))+15.
                    else:
                        dst.variables[name][:] = src.variables[name][:]
                    # copy attributes for this variable
                    for ncattr in src.variables[name].ncattrs():
                        attribute=src.variables[name].getncattr(ncattr)
                        #print('j2',name,ncattr,attribute)
                           
                        if ncattr=='calendar' and exptnamein=='Eoi400':
                            dst.variables[name].setncattr(ncattr,'360_day')
                        else:
                            if (ncattr=='units' and name=='time_counter'):
                            # change units from seconds to days
                                dst.variables[name].setncattr(ncattr,attribute.replace('seconds','days'))
                            else:
                                dst.variables[name].setncattr(ncattr,attribute) 
               
      
    
        
       origf=cf.read_field('temporary.nc')
       print('read copied dataset')
    else:
       if (modelname == 'IPSLCM6A'):
           print (filename)
           origf = cf.read_field(filename, select='sea_surface_temperature')
       elif (modelname == 'CESM1.0.5'):
           #origf = cf.read(filename)
           #print(origf)
           #sys.exit(0)
           reformat_with_iris(exptnamein) # the data is not in a very good format
                                # we will reformat and write to temporary.nc
           
           print('about to read temporaray.nc')
           
           origf=cf.read_field(exptnamein + '_temporary.nc', )
           # iterate over auxillary coordinates and see if they are
           # longitud eor latitude
           for aux in origf.auxs():
               print (origf.auxs(aux))
           a=origf.aux('latitude')
           a.units=('degrees_north')
           print(a)
           
           b=origf.aux('longitude')
           b.units=('degrees_east')
           print(b)
           print(origf)
           #sys.exit(0)  
               
           #for key, aux in self.auxs(ndim=2).iteritems():
           #     if aux.Units.islongitude:
           #sys.exit(0)
           print (origf)
           print (origf.data_axes())
           print (origf.coords())
           
         
           
       else:
            origf=cf.read_field(filename)
        
        
    print (origf)
   
    
    gridf=cf.read_field('one_lev_one_deg.nc')
    #print  gridf
    
    print('about to regrid')
   
    # assume tripoplar grid
    if modelname == 'CESM1.0.5':
       regridf=origf.regrids(gridf,
                          method='bilinear',
                          src_axes={'X': 'ncdim%x','Y':'ncdim%y'},
                          src_cyclic=True)
    else:   
        regridf=origf.regrids(gridf,
                          method='bilinear',
                          src_axes={'X': 'ncdim%x','Y':'ncdim%y'},
                          src_cyclic=True)
    print('regridded')
    
    # see if we can remove auxillary coordinates
    if (modelname=='IPSLCM5A' or modelname=='IPSLCM5A2'):
        regridf.remove_item(description='T',role='a')
   
    
   
    # write to a temporary file so that we can read in as an iris cube and do all our iris analysis 
    # in exactly the same way as before
    
    print(filename)
    print(filenameout)
    cf.write(regridf,filenameout,fmt='NETCDF4_CLASSIC')
 
   
  
    
   
  

#############################################################################
def getnames(modelname,fieldname,exptname): 
    """
    parameters: 
        modelname - the model
        fieldname - the name of the field ie 'tos'
        exptname - the name of the experiment ie [Eoi400]
        
    returns:
        fielduse - the name of the field in the file
        filename - the input file
        filenameout = the output file
    
    this program will get the names of the files and the field for each
    of the models   
    
    """
   
    if modelname == 'ACCESS':
        filename = (filestartin + 'ACCESS/' + exptnamein + '/' + field +
                    'Omon_ACCESS_ESM1-5_')
        if exptnamein == 'LP':
            filename = filename + 'LP_09021231-10011231.nc'
        if exptnamein == 'PI':
            filename = filename + 'pre-industrial_01010116-02001216.nc'
     
        if field  == 'tos':
            fielduse = 'Conservative temperature'
        
        
    return [fielduse,filename]


##########################################################
# main program

filename=' '
linux_win='l'
modelname='ACCESS' # IPSLCM5A, IPSLCM5A2
                   #                      IPSLCM6A
                   # 'CESM.1.0.5 HadGEM3
                   
#modelname = 'IPSLCM5A'

exptname = {
        "E280" : "E280",
        "Eoi400" : "EOI400",
        "E400":"E400",
        "E560": "E560"}

exptname_l = {
        "E280" : "e280",
        "Eoi400" : "eoi400",
        "E400":"e400",
        "E560": "e560"}

# only need fields that are on a tripolar grid like orca
fieldname = {
        "tos": "SST"
        }


# this is regridding where all results are in a single file
fieldnamein=['tos'] # ocean tempeature or sst

exptnamein=['PI']
filestart='/uolstore/Research/b/b0164/Data/PlioMIP3/'
filestartin = filestart
filestartout='/nfs/hera1/earjcti/PLIOMIP3/'
    
    


for expt in range(0,len(exptnamein)):
    for field in range(0,len(fieldnamein)):

        # call program to get model dependent names
        # fielduse, and  filename 
        fielduse, filename, filenameout = (getnames
                                           (modelname,fieldnamein[field],
                                            exptnamein[expt]))
        
        fieldnameout = fieldname.get(fieldnamein[field])
        exptnameout = exptname.get(exptnamein[expt])

        print ('filename',filename)
        
        regrid_data(fieldnamein[field],fieldnameout,exptnamein[expt],exptnameout,
                    filename,modelname,linux_win,fielduse,filenameout)

#sys.exit(0)
