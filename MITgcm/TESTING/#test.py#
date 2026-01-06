#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 09:57:27 2020
testing mitgcm

@author: earjcti
"""

import os
from MITgcmutils import mds
import matplotlib.pyplot as plt
import numpy as np



def plot_output():
    XC = mds.rdmds('XC')
    YC = mds.rdmds('YC')
    Tend = mds.rdmds('T', 3600)
    Tmid = mds.rdmds('T', 3420)
    print(XC)
    print(np.shape(XC))
    print(np.shape(Tend))
    print(Tend[0,:,:])
    plt.contourf(XC, YC, Tend[0,:,:], cmap='RdBu_r')
    plt.colorbar()
    plt.show()
    
def plot_input():
     #MITgcm likes its binary in big endian, float ('>f') or double ('>d')
    
    raw = np.fromfile('input/bathymetry.bin', dtype='>f')
    print(np.shape(raw))
    bathy = np.reshape(raw, (40,90))
    
    print(np.shape(bathy))
    print(bathy)
    lsm = np.where(bathy < 0.0, 1.0, 0.0)

    os.chdir('run')
    XC = mds.rdmds('XC')
    YC = mds.rdmds('YC')
        

    cs = plt.contourf(XC, YC, lsm)
    cbar = plt.colorbar(cs)
    plt.show()

    # Temp has 12 time records, 15 vertical levels, 40 in lat, 90 in lon
    #temp = np.reshape(raw, (12,15,40,90))
    
    
os.chdir('/nfs/hera1/earjcti/MITgcm/verification/tutorial_global_oce_latlon/')
plot_input()
