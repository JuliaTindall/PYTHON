
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
#iris.FUTURE.netcdf_promote=True


# #### Bringing in data and creating cubes

# In[2]:

def individual_model(model):
    FILENAME_EOI400 = ('/nfs/hera1/earjcti/regridded/' + model + '/EOI400.' + FIELD + '.mean_month.nc')
    FILENAME_E400 = ('/nfs/hera1/earjcti/regridded/' + model + '/E400.' + FIELD + '.mean_month.nc')
    FILENAME_E280 = ('/nfs/hera1/earjcti/regridded/' + model + '/E280.' + FIELD + '.mean_month.nc')

    # making and correcting cubes
    if model == 'NorESM1-F' and FIELD == 'SST':
        cube_EOI400 = iris.load_cube(FILENAME_EOI400, 'Ocean surface temperature')
        cube_E400 = iris.load_cube(FILENAME_E400, 'SST')
        cube_E280 = iris.load_cube(FILENAME_E280, 'Ocean surface temperature')
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

    #print(mean_cube_list)
    
    return mean_cube_list


# #### Plotting data

# In[4]:

def plot_area_means(cube_list):
    for i,cube in enumerate(cube_list):
        #print(cube.data)
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
# NorESM1-F
FIELD = 'SST' 

UPPER_LIMIT = [90, 60, 30, 0, -30, -60]
LOWER_LIMIT = [60, 30, 0, -30, -60, -90]

for model in MODELNAMES:
    cube_EOI400, cube_E400, cube_E280 = individual_model(model)


# Creating cube lists before calculating MMM

# In[6]:

EOI400_cubes = iris.cube.CubeList([])
E400_cubes = iris.cube.CubeList([])
E280_cubes = iris.cube.CubeList([])


# In[7]:

EOI400_cubes.append(cube_EOI400)
E400_cubes.append(cube_E400)
E280_cubes.append(cube_E280)


# Calculating MMM

# In[8]:

count=0
newdata=0

for cube in EOI400_cubes:
    newdata = newdata + cube.data
    count=count + 1

newdata = newdata / count

MMM_EOI400 = EOI400_cubes[0].copy(data=newdata)

#print(MMM_EOI400)


# In[9]:

count=0
newdata=0

for cube in E400_cubes:
    newdata = newdata + cube.data
    count=count + 1

newdata = newdata / count

MMM_E400 = E400_cubes[0].copy(data=newdata)

#print(MMM_E400)


# In[10]:

count=0
newdata=0

for cube in E280_cubes:
    newdata = newdata + cube.data
    count=count + 1

newdata = newdata / count

MMM_E280 = E280_cubes[0].copy(data=newdata)

#print(MMM_E280)


# In[11]:

MMM_EOI400_list = get_means(MMM_EOI400)
MMM_E280_list = get_means(MMM_E280)
MMM_E400_list = get_means(MMM_E400) #E400 or #EOI280
    
# calculating FCO2
influence_cube_list_MMM = iris.cube.CubeList([])
for i, cube in enumerate(MMM_EOI400_list):
    influence_MMM = (MMM_E400_list[i].data - MMM_E280_list[i].data) / (cube.data - MMM_E280_list[i].data)
    newcube_MMM = cube.copy(data = influence_MMM)
    influence_cube_list_MMM.append(newcube_MMM)
    
# making anomaly (e.g. EOI400-E280)
anomaly_cube_list_MMM = iris.cube.CubeList([])
for i, cube in enumerate(MMM_E400_list):
    anomaly_MMM = (cube.data - MMM_E280_list[i].data)
    new_cube_MMM = cube.copy(data = anomaly_MMM)
    anomaly_cube_list_MMM.append(new_cube_MMM)

# plotting (mean_cube_list_X, influence_cube_list, anomaly_cube_list)
plot_area_means(influence_cube_list_MMM)




