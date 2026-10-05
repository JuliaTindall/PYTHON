#
#  This program will plot the EP-LP SST on the same scale for all experiments
#  It will also produce a multimodel mean
#
#
import matplotlib.pyplot as plt
import numpy as np
#import cf
import iris
from iris.cube import CubeList
import iris.quickplot as qplt
import iris.plot as iplt
import sys

#modelnames = ['AWI-CM3']
modelnames = ['AWI-CM3','CCSM4_Utr','COSMOS','iCESM','HadCM3','MIROC4m','NorESM1-F']

run_length = {'AWI-CM3':'spun up (orbit_err)',
              'CCSM4_Utr':'~500???',
              'COSMOS':'500yrs',
              'iCESM':'~500????',
              'HadCM3':'4000 (final)',
              'MIROC4m':'4000 (final)',
              'NorESM1-F':'2200 (final)'}
filelocs = '/home/earjcti/PLIOMIP3_results/interim_other_models/'


def get_AWICM3_cube():
    """
    gets the EP-LP SST anomaly for AWICM3
    """
    filename = filelocs + 'AWI/AWI-CM3_EP_SST_annclim_anomaly_vs_LP_360x180.nc'
    icube = iris.load_cube(filename)
    cube = iris.util.squeeze(icube)

    return cube

def get_COSMOS_cube():
    """
    gets the EP-LP SST anomaly for COSMOS
    """
    filename = filelocs + 'AWI/COSMOS_EP_SST_annclim_anomaly_vs_LP_360x180.nc'
    icube = iris.load_cube(filename)
    cube = iris.util.squeeze(icube)

    return cube


def get_CCSM4_Utr_cube():
   
    EP_file = filelocs + 'CCSM4_Utr/Eoi490_Ocean_clim_regrid.nc'
    LP_file = filelocs + 'CCSM4_Utr/Eoi400_Ocean_clim_regrid.nc'
   
    EPallmonths_cube = iris.util.squeeze(iris.load_cube(EP_file))
    LPallmonths_cube = iris.util.squeeze(iris.load_cube(LP_file))

    EP_cube = EPallmonths_cube.collapsed('time',iris.analysis.MEAN)
    EP_cube.data = np.ma.where(EP_cube.data.mask, -99999., EP_cube.data)
    LP_cube = LPallmonths_cube.collapsed('time',iris.analysis.MEAN)
    cube = EP_cube.copy(data=EP_cube.data - LP_cube.data)
    #iris.save([EP_cube,LP_cube,cube],'temporary.nc',fill_value=-99999.)
    cube.remove_coord('time')
     
    return cube


      
def get_iCESM_cube():

    EP_file = filelocs + 'iCESM/SST.iCESM13.EP.100yrmean.nc'
    LP_file = filelocs + 'iCESM/SST.iCESM13.LP_Control.100yrmean.nc'
    EP_cube = iris.load_cube(EP_file,'TEMP')
    LP_cube = iris.load_cube(LP_file,'TEMP')
    cube = EP_cube - LP_cube

    return cube

def get_HadCM3_cube():
    EP_file = filelocs + 'HadCM3/xqbwg_Annual_Average_#pf_SST_3900_4000.nc'
    LP_file = filelocs + 'HadCM3/xqbwd_Annual_Average_#pf_SST_3900_4000.nc'
    grid=iris.load_cube('one_lev_one_deg.nc')

    EP_cube = iris.util.squeeze(iris.load_cube(EP_file))
    LP_cube = iris.util.squeeze(iris.load_cube(LP_file))

    EP_cube_r = EP_cube.regrid(grid,iris.analysis.Linear())
    LP_cube_r = LP_cube.regrid(grid,iris.analysis.Linear())

    cube = EP_cube_r - LP_cube_r

    return cube
    
def get_MIROC_cube():
    EP_file = filelocs + 'MIROC4m/tos_MIROC4m_pliomip3-EP_average_1deg1deg.nc'
    LP_file = filelocs + 'MIROC4m/tos_MIROC4m_pliomip3-LP_average_1deg1deg.nc'

    EP_cube = iris.util.squeeze(iris.load_cube(EP_file))
    LP_cube = iris.util.squeeze(iris.load_cube(LP_file))
    cube = EP_cube - LP_cube
    return cube
    
def get_NorESM1_cube():

    EP_file = filelocs + 'NorESM1-F/EP490_sst_climo.nc'
    LP_file = filelocs + 'NorESM1-F/LP400_sst_climo.nc'

    EPallmonths_cube = iris.util.squeeze(iris.load_cube(EP_file))
    LPallmonths_cube = iris.util.squeeze(iris.load_cube(LP_file))

    EP_cube = EPallmonths_cube.collapsed('time',iris.analysis.MEAN)
    LP_cube = LPallmonths_cube.collapsed('time',iris.analysis.MEAN)
    cube = EP_cube - LP_cube
    return cube


        

def plot_cube(cube):
    """
    plots the cube
    """
    vals = np.arange(-5,5.5,0.5)
    qplt.contourf(cube,levels=vals,extend='both',cmap='RdBu_r')
    plt.gca().coastlines()
   
    
def Utr_reformat_tripolar(exptname):

    """
    we are going to load the data as an iris cube and add some
    auxillary coordinates and then write out to a file called temporary.nc
    
    the data is currently in filename
    """

    #1. reformat it to a better format
    filename = filelocs + 'CCSM4_Utr/' + exptname + '_Ocean_clim.nc'
    fileout = filelocs + 'CCSM4_Utr/' + exptname + '_Ocean_clim_regrid.nc'
    
    origcube = iris.load_cube(filename, 'Sea Surface Temperature')
    latcube = iris.load_cube(filename, 'array of t-grid latitudes')
    loncube = iris.load_cube(filename, 'array of t-grid longitudes')
    

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
    iris.util.promote_aux_coord_to_dim_coord(origcube, 'y')
    iris.util.promote_aux_coord_to_dim_coord(origcube, 'x')
    
 
    
    
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
    
    
    
    iris.save(origcube, 'temporary.nc',fill_value=2.0E20)


    #2.  convert to rectilinear

    origf = cf.read('temporary.nc')[0]
    gridf = cf.read('one_lev_one_deg.nc')[0]

    regridf = origf.regrids(
        gridf,
        method='bilinear',
        src_axes={'X': 'ncdim%x', 'Y': 'ncdim%y'},
        src_cyclic=True)

    #print('regridf',regridf)
    #print('regridf.properties()',regridf.properties())
    #print('regridf.get_property(_FillValue, None)',
    #      regridf.get_property('_FillValue', None))
    cf.write(regridf,fileout,fmt='NETCDF4')
   






allSSTcubes = CubeList([])

# get data

for model in modelnames:
    if model == 'AWI-CM3':
        anomcube = get_AWICM3_cube()
    if model == 'COSMOS':
        anomcube = get_COSMOS_cube()
    if model == 'iCESM':
        anomcube = get_iCESM_cube()
    if model == 'HadCM3':
        anomcube = get_HadCM3_cube()
    if model == 'MIROC4m':
        anomcube = get_MIROC_cube()
    if model == 'NorESM1-F':
        anomcube = get_NorESM1_cube()
    if model == 'CCSM4_Utr':
        #Utr_reformat_tripolar('Eoi490') # if regridding required
        anomcube = get_CCSM4_Utr_cube()
        anomcube.attributes.pop('invalid_units',None)

    allSSTcubes.append(anomcube)

print('allSST')
print(allSSTcubes)

    

# plot data
plt.figure(figsize=[14,5])
for i,cube in enumerate(allSSTcubes):
    plt.subplot(2,4,i+1)
    plot_cube(cube)
    plt.title(modelnames[i]+ ':' + run_length.get(modelnames[i]))

plt.savefig('allplots.png')
plt.savefig('allplots.eps')
plt.close()

#plt.show()
#sys.exit(0)

