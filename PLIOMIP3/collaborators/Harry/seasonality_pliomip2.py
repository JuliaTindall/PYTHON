import numpy as np
import iris
import matplotlib.pyplot as plt
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.cm as cm
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import cartopy.crs as ccrs
import sys

def make_cmap(colors, position=None, bit=False):
    '''
    I didn't write this I found it on the web.
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    import matplotlib as mpl
    import numpy as np
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap


def customise_cmap2():
    """
    as customise_cmap but 19 colors only + 2 white in middle added by Julia
    """

    colors = [(5,48,97),(22,75,124),(39,102,151),(56,130,178),
              (80,154,199),(114,173,209),(147,192,219),(181,211,228),
              (215,230,238),(255,255,255),(255,255,255),(255,255,255),
              (242,220,217),(236,192,185),
              (229,163,153),(223,135,121),(216,107,89),(195,80,69),
              (164,53,56),(133,26,43),(103,0,31)]
    my_cmap = make_cmap(colors, bit=True)
    return my_cmap


def read_harry_data():
    """
    reads harrys data returns an array of lat lon and seasonalinty
    """

    with open('mPWP_seasonality.txt') as f:
        lines = f.readlines()
        
    lats = np.zeros(len(lines)-1)
    lons = np.zeros(len(lines)-1)
    seasonality = np.zeros(len(lines)-1)
    site = np.zeros(len(lines)-1,'U15')
        
    for i,line in enumerate(lines):
        if i >0:  # to avoid title line 
            a,b,c,d = line.split()
            lats[i-1] = b
            lons[i-1] = c
            site[i-1] = a
            seasonality[i-1]=d
       
    return lats,lons,seasonality,site



def get_pliomip2_mean(latreq,lonreq):
    """
    gets the pliomip2 mean for the latreq,lonreq
    """
    filename = ('/uolstore/home/users/earjcti/hera1/regridded/' +
                'SST_multimodelmean_month.nc')
    cube = iris.load_cube(filename,'SSTmean_plio')
    #jancube = cube[0,:,:]
    febcube = cube[1,:,:]
    #marcube = cube[2,:,:]
    #aprcube = cube[3,:,:]
    #maycube = cube[4,:,:]
    #juncube = cube[5,:,:]
    #julcube = cube[6,:,:]
    augcube = cube[7,:,:]

    seascube = augcube.copy(data=np.abs(augcube.data - febcube.data))
    cubelats = seascube.coord('latitude').points
    cubelons = seascube.coord('longitude').points

    modelSST = np.zeros(np.shape(latreq))

    for i,lat in enumerate(latreq):
        if lonreq[i] < 0:
            lon=lonreq[i] + 360.
        else:
            lon=lonreq[i]

        # find nearest longitude and latitude to the value
        latix = (np.abs(cubelats-lat)).argmin()
        lonix = (np.abs(cubelons-lon)).argmin()

        # get data from this location
        data_slice  =  seascube.extract(iris.Constraint(
            latitude = cubelats[latix],longitude = cubelons[lonix]))
        #print(data_slice.data,lonreq[i],lon,lat)
        modelSST[i] = data_slice.data

    return modelSST,seascube  
        
        
    
def plot(seascube,datalon,datalat,dataseas,modelseas):
    """
    plots the cube and overplots the data
    """

    # plot model

    V= np.arange(0.0,16.0,2.0)
    print(V)
    #mycmap = customise_cmap2()
    mycmap=plt.cm.Reds
  

    qplt.contourf(seascube,levels=V,extend='max',cmap=mycmap) 
    #qplt.contourf(seascube,levels=V,cmap=mycmap,extend='max')
    plt.gca().coastlines()
    plt.title('LP: ABS( SST (Aug-Feb))')
    

    # overplot data
  
    norm = colors.BoundaryNorm(boundaries=V,ncolors=mycmap.N)
  
    plt.scatter(datalon, datalat, c='black',  marker='o',
                s=60, transform=ccrs.Geodetic())

    plt.scatter(datalon, datalat, c=dataseas,  marker='o', s=30,
                norm = norm , cmap=mycmap , transform=ccrs.Geodetic())

    plt.savefig('data_model_comparison_abs.png')


# main program
# 1 read data

datalat,datalon,dataseas,datasite = read_harry_data()

# 2 get model output (we want seasonality to return absolute value)

modelseas,seascube = get_pliomip2_mean(datalat,datalon)

# 3 print out

for i,lon in enumerate(datalon):
    print('lon=',lon,'lat=',datalat[i],dataseas[i],round(modelseas[i],2))

# 4 plot points

plot(seascube,datalon,datalat,dataseas,modelseas)

# 5 stats
diffarr = dataseas - modelseas
diffsqsum=0
diffsum=0
count=0
for i,diff in enumerate(diffarr):
    if np.isfinite(diff):
        diffsqsum = diffsqsum + (diff * diff)  # sum of squares
        diffsum = diffsum + diff
        count=count+1
        print(datasite[i],(diffsqsum / count),dataseas[i],round(modelseas[i],2))

rmse = np.sqrt(diffsqsum / count)

print('RMSE =', rmse)
print('mae =', diffsum / count)
