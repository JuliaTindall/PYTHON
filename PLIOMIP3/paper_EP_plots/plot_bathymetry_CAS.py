import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import cartopy.crs as ccrs
import cmocean
import sys

      


filename = '/home/earjcti/um/xpsie/pg/xpsieo#pg000001798c1+.nc'
allcubes = iris.load(filename)

for cube in allcubes:
    print(cube.var_name)
    if cube.var_name == 'insitu_T':
        temp_cube = iris.util.squeeze(cube)
        temp_depth = temp_cube.coord('depth_1').points
    if cube.var_name == 'W':
        W_cube = iris.util.squeeze(cube)
        W_depth = W_cube.coord('depth').points
        W_lats = W_cube.coord('latitude').points
        W_lons = W_cube.coord('longitude').points



nlats = len(W_lats)
nlons = len(W_lons)
W_data = W_cube.data
temp_data = temp_cube.data

#print(W_data.mask)
bottom_data = np.zeros((nlats,nlons))

for j in range(0,nlats):
    for i in range(0,nlons):
        # find where depth becomes nan
        #print(W_data.mask[:,j,i])
        if temp_data.mask[19,j,i]==False:
            bottom_data[j,i]=5500.
        for k in range(0,19):
            if temp_data.mask[k+1,j,i] == True and temp_data.mask[k,j,i] == False:
                #print('here',k,W_depth[k])
                bottom_data[j,i] = W_depth[k]
                #print('found',bottom_data[j,i])
                
    
bottom_depth_cube = (temp_cube[0,:,:].copy(data=bottom_data))

##########################
# now plot CAS region  AI generated
# CAS region (adjust as needed)

boundaries=temp_depth[6:]
#boundaries=temp_depth

      
  
#cmap_terrain = plt.cm.get_cmap('terrain',len(boundaries))
#colors=list(cmap_terrain(np.arange(len(boundaries))))
#cmap_terrain.set_over(colors[-1])  # set to last color
#cmap_terrain.set_under('green')  # set to first color

cmap_ocean = plt.cm.get_cmap('rainbow',len(boundaries)-1)
#cmap_ocean = cmocean.cm.deep.resampled(len(boundaries)-1)
colors=list(cmap_ocean(np.arange(cmap_ocean.N)))
cmap_ocean.set_over(colors[-1])  # set to last color
cmap_ocean.set_under('lightgrey')  # set to first color


reg_EP_cube = bottom_depth_cube.extract(iris.Constraint(latitude = lambda cell: 5< cell < 18, longitude=lambda cell: 265 < cell < 285))
print('data',reg_EP_cube.data)
#print('boundaries',boundaries)
   
ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())

cs=iplt.pcolormesh(reg_EP_cube, cmap=cmap_ocean,
                   norm = matplotlib.colors.BoundaryNorm(boundaries, 
                                            ncolors=len(boundaries)-1, 
                                                         clip=False))
cbar=plt.colorbar(cs,orientation='vertical', extend='neither',spacing='uniform')
cbar.ax.invert_yaxis()


cbar.set_ticks(W_depth[6:])
cbar.set_label('Depth (m)')
#cbar.set_ticks(temp_depth)
#qplt.contour(e280_lsmcube,[0.5],colors='red',linewidths=0.75)
#qplt.contour(eoi400_lsmcube,[0.5],colors='black',linewidths=0.75,linestyle='dashed')
gl=ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels =False
gl.xlocator = mticker.FixedLocator([-93.75, -90.0,-86.25,-82.5, -78.75, -70, -60])
gl.ylocator = mticker.FixedLocator([5, 7.5,10, 12.5, 15,])
ax.coastlines()
#plt.show()
plt.savefig('bathymetry.png')
sys.exit(0)
 
    
