#!/usr/bin/env python
# Mark Richardson, CEMAC
# adapted from work by Robin Steven at SEE Leeds

from pylab import *
from matplotlib.colors import ListedColormap
import netCDF4 as nc4

# --- PARAMETERS ---
# Should check if the basemap exists
try:
    from mpl_toolkits.basemap import Basemap
    dobasemap = True                     
except ImportError:                         
    dobasemap = False                    

# set a destination for the plots
Home="/nfs/see-fs-02_users/earmgr/tmpOmniGlobeContent"
# Need to choose a file
file_a='base2000_cesm111_v2_hourly_avg.nc'
# this is the reference topology JPG
etopo_jpg = '/nfs/see-fs-02_users/earmgr/Images/etopo3000.jpg'

# parse cmd line for choice of variable to inspect
# O3_SRF, ISOP_SRF (but T is 3D so need some filtering to get surface T)
VarChoice = 1

# Target visualisation system
OmniGlobe=True
FullHD =  False
Mon1680 = False
# The number of pixels in the output is related to the "inches" and dots per inch
# this will be set by user in production code
# 100 dpi seems to be default save image setting
# So values are inches to get a set number of pixels 
# when doing FullHD need to think about 80dpi (how to set, why to set?)
if OmniGlobe:
  h_inches = 30.0  # 3000 by 1500 i.e. OmniGlobe at 100dpi
  v_inches = 15.0

elif FullHD:
  h_inches = 19.2   # 1920 by 1080 i.e. FullHD@100dpi
  v_inches = 10.8
  #h_inches = 24.0 v_inches = 13.5 # 1920 by 1080 i.e. FullHD@80dpi
elif Mon1680:
  h_inches = 16.8   # 1680x1050@100dpi
  v_inches = 10.5

# --- INITIALIZE ---
lightlevel = zeros([96,145])  # CAUTION size has been hardwired
#
# sunlight level 0.0 = no sunlight, 1.0 is maximum sunlight (tbd)
night = 0.0
light = 1.0

# At midnight 0/01/2001, adjust these per hour (shift by 6 cells to left - subtracting)
# AGAIN CAUTION DUE TO HARDWIRED SIZE
refdusk = 108  
refdawn = 36
CellsPerDay = 144
hr_shift = CellsPerDay/24   # for different resolution will auto adjust to one hour shift

# this is a mechanism for generating light and dark data
lux = np.linspace(0,144,145,endpoint=True)   # make sure lux is a list
# turn on the light (there is a more efficient method for this I am sure
for i in np.arange(CellsPerDay+1):
  lux[i] = light         # cells 0 to CellsPerDay (+1 for faces )

# Should we allow for choosing different variables within the code or preset it?
if VarChoice == 1:
  ChoiceOfVar = 'O3_SRF'
  VarName = "/o3_surf_"
elif VarChoice == 2 :
  ChoiceOfVar = 'ISOP_SRF'
  VarName = "/isop_s_"
elif VarChoice == 3 :
  ChoiceOfVar = 'T'
  VarName = "/degK_"

if dobasemap:
   map = Basemap(resolution='c', llcrnrlon=0.0,llcrnrlat=-90, urcrnrlon=360.0,urcrnrlat=90)

##### OUR DATA for displaying ###############
ncfile_a = nc4.Dataset(file_a, 'r')

# Extract latitude, longitude
lons = ncfile_a.variables['lon'][:]
lats = ncfile_a.variables['lat'][:]

# Cater for lack of 360.0 information (periodic)
nlats = len(lats)
nlons = len(lons)
lons_p = zeros(nlons+1)
lons_p[:-1] = lons
lons_p[-1] = 360.

# Make a copy of colormap and make it scale from transparent to solid
c_map = cm.Reds
m_map = c_map(np.arange(c_map.N) )
m_map[:,-1] = np.linspace(0,1,c_map.N)
m_map = ListedColormap(m_map)

if dobasemap:
   x,y =  map(*meshgrid(lons_p,lats))

# Assume lats lons persist - avoid reading them repeatedly
# For each timeslice of interest (established in the loop limits and stride)
#   Here I expect to see 30 images at 24 hour interval at midday GMT
start = 702
finish = 727
stride = 1
for simhour in arange(start,finish,stride):
  dlight = mod(simhour,24)  # number of hours from midnight, requires first dataset to midnight
  
  # Extract chosen var surface dat (2D geom)
  chosen_a = ncfile_a.variables[ChoiceOfVar][simhour,:,:]
  
  # Min and Max of this field
  dat_lb = amin(chosen_a)
  dat_ub = amax(chosen_a)
  print "Min, Max values this timeslice are",simhour,dat_lb,dat_ub

  # Cater for lack of 360.0 information (periodic)
  chosen_p = zeros( [nlats, nlons+1] ) # P the plottable version has extra column
  chosen_p[:,:-1] = chosen_a          # copy the content of A into P
  chosen_p[:,-1] = chosen_a[:,0]      # insert first column of A into final column of P
  
  # Need an artificial sunlight representation to sync with "simhour"
  # reference where the terminators exist (naive north-south line)
  # Dusk
  dusk=fmod((refdusk-dlight*hr_shift),CellsPerDay)
  if dusk < 0 : 
    dusk = CellsPerDay + dusk
  # Dawn
  dawn=fmod((refdawn-dlight*hr_shift),CellsPerDay)
  if dawn < 0 : 
    dawn = CellsPerDay + dawn
  
  ###print 'Dlight=',dlight,'Dusk=',dusk,' Dawn=',dawn
  # between dusk and dawn it is dark
  if (dusk < dawn):
    lux[0:dusk] = light
    lux[dusk:dawn] = night
    lux[dawn:CellsPerDay] = light
  else:
    lux[0:dawn] = night
    lux[dawn:dusk] = light
    lux[dusk:CellsPerDay] = night
  
  # now propogate over all latitudes (96 in this specific sim)
  for j in np.arange(nlats):
    for i in np.arange(CellsPerDay+1):
      lightlevel[j][i] = lux[i]
  # END OF setting up data now do the figure planning
  
  # the inches were determined near start of script
  figure(figsize=(h_inches,v_inches))
  
  # make canvas as big as figure area (for OmniGlobe PNGs)
  subplots_adjust(0,0,1,1)
  
  if dobasemap:
    # try put the geography under data
    map.warpimage( image=etopo_jpg )

    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25)
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary()
    # draw lat/lon grid lines every 30 degrees. But omit labels
    ###map.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,0],fontsize=10,linewidth=0.25)
    ###map.drawparallels(np.arange(-90.,120.,30.),labels=[0,0,0,0],fontsize=10,linewidth=0.25)

    # now plot the quantity on tthe map
###    map.pcolormesh(lons_p, lats, chosen_a, cmap=m_map,alpha=0.6)
    map.pcolormesh(lons_p, lats, chosen_a, cmap=m_map, edgecolors=face)
    # indicate the night and day regions
    map.pcolormesh(lons_p, lats, lightlevel , cmap=cm.gray,alpha=0.2)

  else:
    print "Basemap not found plotting anyway"
    # Now to plot the actual data
    pcolormesh(lons_p, lats, chosen_p , cmap=cm.Purples)
    # indicate the night and day regions
    pcolormesh(lons_p, lats, lightlevel , cmap=cm.gray,alpha=0.2)
  
  plot()    # Render the image (internally without drawing it
  
  # should do better to change this and introduce checking
  TimeStamp = "2001-01-01"   # this should be in the nc file
  Frame = str("%04d"%simhour)
  MetaData = VarName+TimeStamp+"_"+Frame
  fname_out = Home+MetaData+".png"
  
  savefig(fname_out,bbox_inches="tight",pad_inches=0.0)

  # Iteration end one figure per timestep
  print "Saved figure ",fname_out
  if simhour > finish-stride-2:
    show()                        # THIS IS NOT SENSIBLE in production
    #draw()    # not sure if this is right invocation
  clf()    # clear the figure but leave the canvas
  # still trying to decide how to use next as I would like to see at least one figure
  close()  # could be explicit with figure name or number

###show() # Show final figure (I hope)
# END OF SCRIPT
