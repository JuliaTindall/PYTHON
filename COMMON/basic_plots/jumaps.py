#
#  contour globe will plot data on a global map
#
def contourglobe(plotdata,subplotx,subploty,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname,res):

    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid



    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(subplotx,subploty,fileno+1)

   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution=res)
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary
    x, y = map(lons, lats)
    map.drawcoastlines() 
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.),extend="both")
        cbar = plt.colorbar(cs,orientation="horizontal",extend='max')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu_r',extend="both")
            cbar = plt.colorbar(cs,orientation="horizontal",extend='max')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend="both")
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                cs = map.contourf(x,y,plotdata,V,extend="both")
                cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
#end def contourglobe
