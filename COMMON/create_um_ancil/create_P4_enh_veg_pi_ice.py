#NAME
#    create_P4_enh_veg_pi_ice.py
#PURPOSE 
#
#  This program will create a netcdf file that can be inputted into xancil to create
#  a vegetation ancil file.
#  The file will have P4 vegetation types but E280 ice sheets.

# Import necessary libraries

import numpy as np
import sys
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")
 
    
#=================================================================
# MAIN PROGRAM STARTS HERE

# read in the data
EOI400_vegfile = '/nfs/hera1/earjcti/ancil/P4_enh/P4_enh_mb_qrfrac.type.nc'
E280_vegfile = '/nfs/hera1/earjcti/ancil/preind2/qrfrac.type.nc'

cubes = iris.load(E280_vegfile)

EOI400_cube_orig = iris.load_cube(EOI400_vegfile)
E280_cube_orig = iris.load_cube(E280_vegfile)

EOI400_cube = iris.util.squeeze(EOI400_cube_orig)
E280_cube = iris.util.squeeze(E280_cube_orig)

pi_icecube = E280_cube[0,:,:]  # ice is on level 1 because of how it was setup

EOI400_data = EOI400_cube.data
nlev,nlat,nlon = np.shape(EOI400_data)

# find a location where the PI_ice  is different to the Plio_ice
for j in range(0,nlat):
    for i in range(0,nlon):
        if EOI400_data[8,j,i] !=pi_icecube.data[j,i]:
            if pi_icecube.data[j,i] == 1.0:
                EOI400_data[8,j,i] = pi_icecube.data[j,i]
                EOI400_data[0:7,j,i] = 0.0
            else:
                print('to change',i,j,EOI400_data[8,j,i],pi_icecube.data[j,i])
                sys.exit(0)
        if EOI400_data[8,j,i] > 1.0E36:
            EOI400_data[:,j,i] = -999.999


# save new vegetation file
new_cube = EOI400_cube.copy(data = EOI400_data)
iris.save(new_cube,'qrfrac_P4_lsm_pi_ice.nc',fill_value=-999.999)

#



