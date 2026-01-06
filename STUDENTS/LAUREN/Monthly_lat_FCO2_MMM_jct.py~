
# coding: utf-8

# # Seasonal and latitudinal line plots for $FCO_2$
# 
# This code is for creating line plots showing the seasonal and latitudinal variation in FCO2 for SAT, SST and precipitation

# #### Importing libraries

# In[1]:

import numpy as np
import numpy.ma as ma

import matplotlib as mpl
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt
import iris.plot as iplt
import sys
#iris.FUTURE.netcdf_promote=True


# #### Bringing in data and creating cubes

# In[2]:

def individual_model(model):
    FILENAME_EOI400 = ('/nfs/hera1/earjcti/regridded/' + model + '/EOI400.' + FIELD + '.mean_month.nc')
    FILENAME_E400 = ('/nfs/hera1/earjcti/regridded/' + model + '/E400.' + FIELD + '.mean_month.nc')
    FILENAME_E280 = ('/nfs/hera1/earjcti/regridded/' + model + '/E280.' + FIELD + '.mean_month.nc')

    # making and correcting cubes
    if model == 'NorESM1-F' and FIELD == 'SST':
        print('in NorESM')
        cube_EOI400 = iris.load_cube(FILENAME_EOI400, 'Ocean surface temperature')
        cube_E400 = iris.load_cube(FILENAME_E400, 'SST')
        cube_E280 = iris.load_cube(FILENAME_E280, 'Ocean surface temperature')
        cube_EOI400.coord('time').rename('month')
        cube_EOI400.coord('month').points=[1,2,3,4,5,6,7,8,9,10,11,12]
        cube_E280.coord('time').rename('month')
        cube_E280.coord('month').points=[1,2,3,4,5,6,7,8,9,10,11,12]
        cube_EOI400.coord('month').attributes=None
        cube_E280.coord('month').attributes=None
        cube_E400.coord('month').attributes=None
        cube_EOI400.coord('month').units=None
        cube_E280.coord('month').units=None
        cube_E400.coord('month').units=None
       
    else:
        cube_EOI400 = iris.load_cube(FILENAME_EOI400, FIELD)
        cube_E400 = iris.load_cube(FILENAME_E400, FIELD)
        cube_E280 = iris.load_cube(FILENAME_E280, FIELD)
    
    cube_E400 = cube_E280.copy(data=cube_E400.data)

    return cube_EOI400, cube_E400, cube_E280


# #### Setting coords and latitude bands

# In[3]:

def get_means(cube):
    try:
        cube.coord('latitude').guess_bounds() 
        cube.coord('longitude').guess_bounds()
    except:
        print('no bounds')
    
    grid_areas = iris.analysis.cartography.area_weights(cube)
    mean_cube_list = iris.cube.CubeList([])

    for i, upper in enumerate(UPPER_LIMIT):
        grid_areas_region = iris.analysis.cartography.area_weights(cube)
        print(UPPER_LIMIT[i], LOWER_LIMIT[i], i)
        for j, lat in enumerate (cube.coord('latitude').points):
            if lat > upper or lat < LOWER_LIMIT[i]:   # > then < for combined
                grid_areas_region[:, j, :] = 0.0

        zonalmean = cube.collapsed(['longitude', 'latitude'],
                            iris.analysis.MEAN,
                            weights = grid_areas_region)
        mean_cube_list.append(zonalmean)

    
    return mean_cube_list


# #### Plotting data

# In[4]:

def plot_area_means(cube_list):
    for i,cube in enumerate(cube_list):
        print(cube.coord('month').points)
        if LOWER_LIMIT[i] < 0:
            cube_label = np.str(UPPER_LIMIT[i]*(-1)) + '°S - ' + np.str(LOWER_LIMIT[i]*(-1)) + '°S'
        else:
            cube_label = np.str(LOWER_LIMIT[i]) + '°N - ' + np.str(UPPER_LIMIT[i]) + '°N'    # N for NH or S for SH
        qplt.plot(cube, label=cube_label)

    plt.title('$FCO_2$ on SST')
    plt.legend(bbox_to_anchor=(1,1.02), loc=2)
    #plt.grid(True, linestyle='-')
    plt.xlabel('Month') 
    plt.ylabel('$FCO_2$')
    plt.ylim(0,1,0.1)
    
    iplt.show()
    #plt.savefig('Monthly_lat_FCO2_MMM_excl_NorESM1-F.png', dpi=300, bbox_inches='tight')
    plt.close()


# #### Main program

# Defining models and upper/lower limits

# In[5]:

MODELNAMES = ['HadCM3', 'CESM2', 'MIROC4m', 'CCSM4-UoT', 'COSMOS', 'NorESM1-F']
#MODELNAMES = ['NorESM1-F']
#MODELNAMES = ['CCSM4-UoT']
# NorESM1-F
FIELD = 'SST' 

UPPER_LIMIT = [90, 60, 30, 0, -30, -60]
LOWER_LIMIT = [60, 30, 0, -30, -60, -90]

EOI400_cubes = iris.cube.CubeList([])
E400_cubes = iris.cube.CubeList([])
E280_cubes = iris.cube.CubeList([])


for model in MODELNAMES:
    cube_EOI400, cube_E400, cube_E280 = individual_model(model)
    EOI400_cubes.append(cube_EOI400)
    E400_cubes.append(cube_E400)
    E280_cubes.append(cube_E280)


# Calculating MMM EOI400

count=0
newdata=0

for cube in EOI400_cubes:
    newdata = newdata + cube.data
    count=count + 1
newdata = newdata / count
MMM_EOI400 = EOI400_cubes[0].copy(data=newdata)

# calculating MMM E400

count=0
newdata=0

for cube in E400_cubes:
    newdata = newdata + cube.data
    count=count + 1

newdata = newdata / count

MMM_E400 = E400_cubes[0].copy(data=newdata)

# calculating MMM E280

count=0
newdata=0

for cube in E280_cubes:
    newdata = newdata + cube.data
    count=count + 1

newdata = newdata / count

MMM_E280 = E280_cubes[0].copy(data=newdata)

#print(MMM_E280)
  
#plt.subplot(311)
#qplt.contourf(MMM_EOI400[0,:,:])
#plt.subplot(312)
#qplt.contourf(MMM_E400[0,:,:])
#plt.subplot(313)
#qplt.contourf(MMM_E280[0,:,:])
#plt.show()
#sys.exit(0)
  


# In[11]:

MMM_EOI400_list = get_means(MMM_EOI400)
print('EOI400',MMM_EOI400_list[0].data)

MMM_E280_list = get_means(MMM_E280)
print('E280',MMM_E280_list[0].data)

MMM_E400_list = get_means(MMM_E400) #E400 or #EOI280
print('E400',MMM_E400_list[0].data)

#plt.subplot(311)
#plt.plot(MMM_EOI400_list[0].data)
#plt.subplot(312)
#plt.plot(MMM_E280_list[0].data)
#plt.subplot(313)
#plt.plot(MMM_E400_list[0].data)
#plt.show()
#sys.exit(0)



    
# calculating FCO2
fco2_cube_list_MMM = iris.cube.CubeList([])
for i, cube in enumerate(MMM_EOI400_list):
    fco2_MMM = (MMM_E400_list[i].data - MMM_E280_list[i].data) / (cube.data - MMM_E280_list[i].data)
    newcube_MMM = cube.copy(data = fco2_MMM)
    fco2_cube_list_MMM.append(newcube_MMM)
    
# making anomaly (e.g. EOI400-E280)
anomaly_cube_list_MMM = iris.cube.CubeList([])
for i, cube in enumerate(MMM_E400_list):
    anomaly_MMM = (cube.data - MMM_E280_list[i].data)
    new_cube_MMM = cube.copy(data = anomaly_MMM)
    anomaly_cube_list_MMM.append(new_cube_MMM)

# plotting (mean_cube_list_X, fco2_cube_list, anomaly_cube_list)
plot_area_means(fco2_cube_list_MMM)




