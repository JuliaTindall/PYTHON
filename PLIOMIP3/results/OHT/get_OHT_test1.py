
# redo to use theta rather than deltaT


# we are going to calculate the ocean heat transport over the full depth
# of the ocean. I think this is fine as temperature gradients are not very
# large at the deepest layers
import iris
import numpy as np
import matplotlib.pyplot as plt
import sys

expt='xqbwc'
startyear=3979
endyear=4000



def get_OHT(expt,year):
    """
    Calculates ocean heat transport in Peta Watts for the year
    """
    density=1025
    cp=4000

    filename = ('/home/earjcti/um/'+expt+'/pg/' +expt + 'o#pg00000'
                +str(year) + 'c1+.nc')

    vcube = iris.load_cube(filename,'TOTAL OCEAN V-VELOCITY      CM S**-1')
    vcube = iris.util.squeeze(vcube) / 100. # convert to m/s
    Tcube_origgrid = iris.util.squeeze(iris.load_cube(filename,
                    'POTENTIAL TEMPERATURE (OCEAN)  DEG.C'))
    Tcube_origgrid = Tcube_origgrid + 273.15 # convert to kelvin
    Tcube = Tcube_origgrid.regrid(vcube,iris.analysis.Linear())

    # get stuff we need for area
    depths = vcube.coord('depth_1').points
    lons = vcube.coord('longitude').points
    lats = vcube.coord('latitude').points
    lonres=lons[1]-lons[0]

    depdiff=np.zeros(20)
    top_lev = np.zeros(20)
    low_lev = np.zeros(20)
    depths2=np.zeros(21)

    for k in range(0,19):
        low_lev[k]=(depths[k] + depths[k+1]) / 2.0
    low_lev[19]=depths[19]+((depths[19]-depths[18])/2.0)
    top_lev[1:]=low_lev[:19]
    depdiff = low_lev - top_lev
    
    lats = vcube.coord('latitude').points
    coslats = np.cos(lats * 2.0 * np.pi / 360.) # cos latitude of vacross
    onedeg=111100.   # 111.1 km

    #print(depdiff)
    #print(coslats)
    #sys.exit(0)
  
    area_arr = np.zeros_like(vcube.data,dtype=float)
    area_arr[:]=depdiff[:,np.newaxis,np.newaxis] # add depth
    area_arr[:]=(area_arr[:] *
                 onedeg * coslats[np.newaxis,:,np.newaxis] * lonres)
    area_cube = vcube.copy(data=area_arr)
    area_cube.long_name = 'area'
    
    Qcube = density * cp * vcube * Tcube * area_cube
    Qcube = Qcube / 1.0E15
    Qcube.data = np.ma.where(Qcube.data.mask,-99999.,Qcube.data)
    Qcube.data = np.ma.masked_where(Qcube.data==-99999.,Qcube.data)
    Qcube.units = 'PW'
    Qcube.long_name = 'Ocean heat transport petawatts'
    Tcube.data = np.ma.where(Tcube.data.mask,-99999.,Tcube.data)
    Tcube.long_name = 'potT'

    # find depth and longitude integrated Q/V.  V should be zero
    #Qcube_zm = Qcube.collapsed('longitude',iris.analysis.SUM)
    #Qcube_zm.long_name = 'Ocean heat transport (PW) summed over longitude'
    #Vcube_zm = vcube.collapsed('longitude',iris.analysis.MEAN)
    #Vcube_zm.long_name = 'V zonal mean'
    #iris.save([Qcube_zm,Qcube,vcube,Tcube,area_cube,Vcube_zm],'temporary.nc',
    #           fill_value=-99999.)
  
    #Qcube_int = Qcube_zm.collapsed('depth_1',iris.analysis.SUM)
    #Qcube_int.long_name = 'Depth integrated Ocean heat transport (PW)'
    Q_int = Qcube.collapsed(['longitude','depth_1'],iris.analysis.SUM)
    Q_int.long_name = 'Q_int_mean_dep_lon'
    Q_int_zm = Qcube.collapsed(['longitude'],iris.analysis.SUM)
    Q_int_zm.long_name = 'Q_int_zonal_mean'
    Q_int_surf = Qcube[0:10,:,:].collapsed(['longitude','depth_1'],iris.analysis.SUM)
    Q_int_surf.long_name = 'Q_int_surf'
    V_int_surf = vcube[0:10,:,:].collapsed(['longitude','depth_1'],iris.analysis.SUM)
    V_area = vcube.copy(data=vcube.data * area_arr)
    #V_int = vcube.collapsed(['longitude','depth_1'],iris.analysis.MEAN,weights=area#_arr)
    V_int = vcube.collapsed(['longitude'],iris.analysis.MEAN)
    V_int.long_name = 'mean V by longitude no weights'
    Vdz =V_area.collapsed('depth_1',iris.analysis.SUM)
    Vdz.long_name = 'Vdz'
    Vdz_mean =Vdz.collapsed('longitude',iris.analysis.SUM)
    Vdz_mean.long_name = 'Vdz lon mean'
    
    iris.save([Qcube,Tcube,Vdz,Vdz_mean,Q_int,Q_int_surf,Q_int_zm,V_int,vcube,area_cube,V_area],'temporary.nc',fill_value=-99999.)
    #sys.exit(0)
    
    lats=Q_int.coord('latitude').points

    plt.subplot(121)
    plt.plot(Q_int.data,lats)
    plt.vlines(0,-90,90)
    plt.title('Qint')
    plt.subplot(122)
    plt.plot(Q_int_surf.data,lats)
    plt.vlines(0,-90,90)
    plt.title('qintsurf')
    
            
    
    #for i, lat in enumerate(lats):
    #    print(lat,Qcube_int.data[i])
    plt.show()
    sys.exit(0)
    return lats,Qcube_int.data


years = []
lats_allyears = []
OHT_allyears = []
for year in range(startyear,endyear):
    lats,OHT = get_OHT(expt,year)
    years.append(year)
    lats_allyears.append(lats)
    OHT_allyears.append(OHT)

# save to file
filestart = '/home/earjcti/um/' + expt  
fileout = filestart + '/timeseries/OHT_' + expt + '_' + str(startyear) + '_' + str(endyear) + '.tex'
    
with open(fileout, 'w') as f:
    line = f"year," + ",".join(f"{lat}" for lat in lats)
    f.write(line + "\n")
    # f.write(f"year,{lats}\n")
    for x, ys in zip(years, OHT_allyears):
        line = f"{x}," + ",".join(f"{y}" for y in ys)
        f.write(line + "\n")
    f.close()
