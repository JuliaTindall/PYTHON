#python program
#
# this program will plot the diagnostics associated with ocean heat transport
#
# will do
# 1. mean OHT over a given number of years and anomaly between them for all
# basins.  will put in one figure
# 
import matplotlib.pyplot as plt
import iris
from iris.cube import CubeList
import numpy as np
from matplotlib.gridspec import GridSpec
import sys

#exptnames = ['xpsic','xpsid','xpsig','xpsie']
exptnames = ['xpsie','xpsig']
timeseries_req='n'
yearstart=1400
yearend=1600



# Function to read data from a netcdf file
def get_data(filestart,yearstart,yearend):
    """
    this will read the data from a netcdf files according to yearstart,yearend
    note that the netcdf files will be in chunks of 1000 years
    I am assume that the date in the file corresponds to the calendar year
    """

    Atl_allcubes = CubeList([])
    Pac_allcubes = CubeList([])
    Glob_allcubes = CubeList([])

    if yearstart < 2000 and yearend > 1000:
        try:
            Atlcube = iris.load_cube(filestart + '_1000_2000Tref_is_zero.nc',
                                  'Ocean heat transport (atlantic)')
            Paccube = iris.load_cube(filestart + '_1000_2000Tref_is_zero.nc',
                                  'Ocean heat transport (pacific)')
            Globcube = iris.load_cube(filestart + '_1000_2000Tref_is_zero.nc',
                                  'Ocean heat transport (global)')
        except:
            Atlcube = iris.load_cube(filestart + '_1000_2000.nc',
                                  'Ocean heat transport (atlantic)')
            Paccube = iris.load_cube(filestart + '_1000_2000.nc',
                                  'Ocean heat transport (pacific)')
            Globcube = iris.load_cube(filestart + '_1000_2000.nc',
                                  'Ocean heat transport (global)')

        Atl_allcubes.append(Atlcube)
        Pac_allcubes.append(Paccube)
        Glob_allcubes.append(Globcube)
    else:
        print('add more years - already set up to cross files')
        sys.exit(0)
  

    time_constraint = iris.Constraint(t=lambda
                      cell: yearstart <= cell.point.year <= yearend)


    fullcube = Atl_allcubes.concatenate_cube()
    Atl_timecube = fullcube.extract(time_constraint)
    Atl_meancube = Atl_timecube.collapsed('t',iris.analysis.MEAN)

    fullcube = Pac_allcubes.concatenate_cube()
    Pac_timecube = fullcube.extract(time_constraint)
    Pac_meancube = Pac_timecube.collapsed('t',iris.analysis.MEAN)

    fullcube = Glob_allcubes.concatenate_cube()
    Glob_timecube = fullcube.extract(time_constraint)
    Glob_meancube = Glob_timecube.collapsed('t',iris.analysis.MEAN)

    return Atl_meancube, Pac_meancube, Glob_meancube



##################################################################
period = {'xpsid':'LP','xpsij':'LP490','xpsie':'EP$_{400}$','xpsig':'EP',
          'xpsic':'PI','xqbwg':'EP','xqbwj':'LP490'}


filestart = '/home/earjcti/um/xpsig/timeseries/OHT_xpsig'
Atl_cube_g, Pac_cube_g, Glob_cube_g =get_data(filestart,yearstart,yearend)

filestart = '/home/earjcti/um/xpsie/timeseries/OHT_xpsie'
Atl_cube_e, Pac_cube_e, Glob_cube_e =get_data(filestart,yearstart,yearend)




fig = plt.figure(figsize=[10.5,9])
spec = GridSpec(ncols=4, nrows=2,
                         width_ratios=[1, 1,1,1], wspace=0.1,
                         hspace=0.5, height_ratios=[1, 1])

# plot Atlantic

ax=fig.add_subplot(spec[0,0:2])
latitudes = Atl_cube_e.coord('latitude').points
values_g = Atl_cube_g.data
values_e = Atl_cube_e.data
values_anom = values_g - values_e

ax.plot(latitudes, values_g, label=period.get('xpsig'),linewidth=2)
ax.plot(latitudes, values_e, label=period.get('xpsie'),linewidth=2,linestyle='--')

ax2=ax.twinx()
ax2.plot(latitudes,values_anom,label=period.get('xpsig') +  ' - ' +
         period.get('xpsie'),linewidth=2,color='grey',
         linestyle='dotted')

#anomvals = (valuse_anom / values) * 100.
#for ix,anom in enumerate(anomvals):
#if -85 < latitudes[ix] < 0.0:
#print(latitudes[ix],round(anom,0),round(values[ix],2),round(cntl_values[ix],2),round((values[ix]-cntl_values[ix]),2))
#                sys.exit(0)
                
                
ax.set_xlabel('Latitude',fontsize=16)
ax.set_ylabel('OHT (PW)',fontsize=16)
ax.tick_params(axis='both',which='major',labelsize=14)

   
ax.axhline(y=0, color='black', linewidth=2)
   
# sort out legend
    
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

lines = lines1 + lines2
labels = labels1 + labels2


ax.set_ylim(-1.5,1.5)
ax2.set_ylim(-0.15,0.15)
ax2.set_yticklabels([])
#ax2.tick_params(axis='y', colors='grey',labelsize=14)
#ax2.spines['right'].set_color('grey')
#ax2.set_ylabel('anomaly (PW)',fontsize=16,color='grey')

ax.grid(True)

ax.set_title('a) Atlantic Ocean Heat Transport',fontsize=16)
ax.set_xlim(-30,90)
#ax.legend(lines, labels, loc='lower right',fontsize=16)



# plot Pacific

ax=fig.add_subplot(spec[0,2:4])
latitudes = Pac_cube_e.coord('latitude').points
values_g = Pac_cube_g.data
values_e = Pac_cube_e.data
values_anom = values_g - values_e

ax.plot(latitudes, values_g, label=period.get('xpsig'),linewidth=2)
ax.plot(latitudes, values_e, label=period.get('xpsie'),linewidth=2,linestyle='--')

ax2=ax.twinx()
ax2.plot(latitudes,values_anom,label=period.get('xpsig') +  ' - ' +
         period.get('xpsie'),linewidth=2,color='grey',
         linestyle='dotted')

                
ax.set_xlabel('Latitude',fontsize=16)
#ax.set_ylabel('OHT (PW)',fontsize=16)
ax.tick_params(axis='both',which='major',labelsize=14)

   
ax.axhline(y=0, color='black', linewidth=2)

ax.set_yticklabels([])
ax.set_ylim(-1.5,1.5)
ax2.set_ylim(-0.15,0.15)
ax2.tick_params(axis='y', colors='grey',labelsize=14)
ax2.spines['right'].set_color('grey')
ax2.set_ylabel('anomaly (PW)',fontsize=16,color='grey')

ax.grid(True)

ax.set_title('b) Pacific Ocean Heat Transport',fontsize=16)
ax.set_xlim(-30,90)
#ax.legend(lines, labels, loc='lower right',fontsize=16)




# plot Global

ax=fig.add_subplot(spec[1,0:3])
latitudes = Glob_cube_e.coord('latitude').points
values_g = Glob_cube_g.data
values_e = Glob_cube_e.data
values_anom = values_g - values_e

ax.plot(latitudes, values_g, label=period.get('xpsig'),linewidth=2)
ax.plot(latitudes, values_e, label=period.get('xpsie'),linewidth=2,linestyle='--')

ax2=ax.twinx()
ax2.plot(latitudes,values_anom,label=period.get('xpsig') +  ' - ' +
         period.get('xpsie'),linewidth=2,color='grey',
         linestyle='dotted')

                
ax.set_xlabel('Latitude',fontsize=16)
ax.set_ylabel('OHT (PW)',fontsize=16)
ax.tick_params(axis='both',which='major',labelsize=14)

   
ax.axhline(y=0, color='black', linewidth=2)

ax.set_ylim(-1.5,1.5)
ax2.set_ylim(-0.15,0.15)
ax2.tick_params(axis='y', colors='grey',labelsize=14)
ax2.spines['right'].set_color('grey')
ax2.set_ylabel('anomaly (PW)',fontsize=16,color='grey')

ax.grid(True)

ax.set_title('c) Global Ocean Heat Transport',fontsize=16)
ax.set_xlim(-90,90)
ax.legend(lines, labels, loc='lower right',fontsize=16,bbox_to_anchor=[1.50,0.6])



plt.tight_layout() 
#plt.show()
#sys.exit()
plt.savefig('supp_figure3.png')
plt.close()
sys.exit(0)

