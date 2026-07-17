#python program
#
# all plots for figure 4 of the paper
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import iris
import iris.quickplot as qplt
import iris.plot as iplt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import sys


##################################################################
# set up figure and subplot
fig = plt.figure(figsize=[8,8])

spec = gridspec.GridSpec(ncols=2, nrows=4,hspace=0.1,
                         wspace=0.1, width_ratios=[1, 1],
                         height_ratios=[1,0.2,1,0.5])


period = {'xqbwd':'LP','xqwbj':'LP$_{490}','xqbwe':'EP$_{400}$',
          'xqbwg':'EP','xqbwc':'PI'}
filestart = '/home/earjcti/um/'

# read in data
file_LP = (filestart + 'xqbwd/MOC/xqbwdpkmean_3900_3999.nc')
file_LP490 = (filestart + 'xqbwj/MOC/xqbwjpkmean_3900_3999.nc')
file_EP = (filestart+ 'xqbwg/MOC/xqbwgpkmean_3900_3999.nc')
file_EP400 = (filestart + 'xqbwe/MOC/xqbwepkmean_3900_3999.nc')

field = 'Meridional Overturning Stream Function (Atlantic)'

LPmoc_cube = iris.load_cube(file_LP,field)
EPmoc_cube = iris.load_cube(file_EP,field)
LP490moc_cube = iris.load_cube(file_LP490,field)
EP400moc_cube = iris.load_cube(file_EP400,field)

#temperature xqbwg-xqbwj Atlantic
ax0=fig.add_subplot(spec[0])

constraint = iris.Constraint(latitude = lambda cell: cell > -30)
moc_cube = LPmoc_cube.extract(constraint)
max_strength = str(np.round(np.max(moc_cube.data)))
print('LP max strength=',max_strength)
cs=iplt.contourf(moc_cube,cmap='RdBu_r',
              levels=np.arange(-25.,30.,5.),extend='both')

ax0.set_title('a) LP',fontsize=14)
plt.ylabel('depth (m)',fontsize=14)
#plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax0.set_xticklabels([])
ax0.set_xticks([-30,0,30,60,90])


# do LP490

ax0=fig.add_subplot(spec[1])

constraint = iris.Constraint(latitude = lambda cell: cell > -30)
moc_cube = LP490moc_cube.extract(constraint)
max_strength = str(np.round(np.max(moc_cube.data)))
print('LP490 max strength=',max_strength)

cs=iplt.contourf(moc_cube,cmap='RdBu_r',
              levels=np.arange(-25.,30.,5.),extend='both')

ax0.set_title('b) LP$_{490}$',fontsize=14)
#plt.ylabel('depth (m)',fontsize=14)
#plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax0.set_yticklabels([])
ax0.set_xticks([0,30,60,90])

# do EP

ax0=fig.add_subplot(spec[4])

constraint = iris.Constraint(latitude = lambda cell: cell > -30)
moc_cube = EPmoc_cube.extract(constraint)
for j,lat in enumerate(moc_cube.coord('latitude').points):
    if lat < 8.:
        moc_cube.data[:,j] = np.nan
max_strength = str(np.round(np.max(moc_cube.data)))
print('eP max strength=',max_strength)

cs=iplt.contourf(moc_cube,cmap='RdBu_r',
              levels=np.arange(-25.,30.,5.),extend='both')

ax0.set_title('c) EP',fontsize=14)
plt.ylabel('depth (m)',fontsize=14)
plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax0.set_xticklabels([])
ax0.set_xticks([0,30,60,90])


# do EP400

ax0=fig.add_subplot(spec[5])

constraint = iris.Constraint(latitude = lambda cell: cell > -30)
moc_cube = EP400moc_cube.extract(constraint)
for j,lat in enumerate(moc_cube.coord('latitude').points):
    if lat < 8.:
        moc_cube.data[:,j] = np.nan
max_strength = str(np.round(np.max(moc_cube.data)))
print('eP max strength=',max_strength)

cs=iplt.contourf(moc_cube,cmap='RdBu_r',
              levels=np.arange(-25.,30.,5.),extend='both')

ax0.set_title('d) EP$_{400}$',fontsize=14)
#plt.ylabel('depth (m)',fontsize=14)
plt.xlabel('latitude ($^\circ$)',fontsize=14)
plt.yticks(fontsize=14)
plt.xticks(fontsize=14)
ax0.set_xticks([0,30,60,90])
ax0.set_yticklabels([])



#colorbar
cax_holder=fig.add_subplot(spec[3,0:3])
cax_holder.axis("off")
cax = inset_axes(
    cax_holder,
    width="100%",   # fraction of the subplot width
    height="30%",  # fraction of the subplot height
    loc="lower center"
)
cbar=plt.colorbar(cs, cax=cax, orientation="horizontal")
cbar.set_label("Sv", fontsize=14)
cbar.ax.tick_params(labelsize=14)

cax_holder.spines[:].set_visible(False)
cax_holder.tick_params(axis='both', which='both', bottom=False, top=False,
                  right=False, left=False,labelbottom=False,
                  labeltop=False,labelleft=False,labelright=False)

plt.tight_layout()
#plt.show()
#sys.exit(0)
plt.savefig('AMOC_plots.png')


# Save the plot to a file
#plt.savefig('seaice/sea-ice_30yr_2_2999.png')
