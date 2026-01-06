#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 29.03.2022 by Julia

We will use the joint pdf's creates by CFMIP to try and attempt to 
"""
import numpy as np
import iris
#from iris.experimental.equalise_cubes import equalise_attributes
import iris.quickplot as qplt
import matplotlib.pyplot as plt

# stuff for kmeans clustering
#from kneed import KneeLocator
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
#from sklearn.metrics import silhouette_score
#from sklearn.preprocessing import StandardScalar

import sys

  
def get_pdfcube():
    """
    gets the pdfcube
    """
    cubes = iris.load(FILEIN)
    for cube in cubes:
        if cube.var_name=='clisccp':
            jointpdf_cube = cube
     
    cloudtop_press = jointpdf_cube.coord('air_pressure').points
    thickname = 'atmosphere_optical_thickness_due_to_cloud'
    optical_depth = jointpdf_cube.coord(thickname).points
          
    return jointpdf_cube.data,cloudtop_press, optical_depth

def get_alt_albedo(joint_pdf_cube,alt_albedo_cube):
    """
    this will estimate the albedo from the jointpdf cube instead of the file
    """
    # albedos corresponding to the optical depths in the pdf boxes
    # from williams and webb 2008
    albedos = [0.028, 0.107, 0.232, 0.407, 0.626, 0.828, 0.950]
    ctp = joint_pdf_cube.coord('air_pressure').points
    tau = joint_pdf_cube.coord('atmosphere_optical_thickness_due_to_cloud').points
    # when calculating albedo it is unclear whether we should average over all points in the histogram or only over those where there is cloud.  
    # I have decided to only average over those where there is cloud, because it is cloud albedo not total albedo

    nz,ny,nx = joint_pdf_cube.data.shape # nz =ctp, ny =optical depth, nx = loc
    for i in range(0,nx): # over locs'
        joint_pdf = joint_pdf_cube.data[:,:,i]
        meanalbedo=0.0
        for j in range(0,nz): # over ctp
            for k in range(0,ny): # over optical depth
                if joint_pdf[j,k] > 0.0:
                    meanalbedo = (meanalbedo +
                                  (joint_pdf[j,k] * albedos[k]))
        alt_albedo_cube.data[i] = meanalbedo / np.sum(joint_pdf)
    
    return alt_albedo_cube

def get_pdfcube_tropics():
    """
    gets the pdfcube and extracts the tropical points only
    """
    cubes = iris.load(FILEIN)
    for cube in cubes:
        print(cube.var_name)
        if cube.var_name=='clisccp':
            jointpdf_cube = cube
        if cube.var_name == 'latitude':
            latitude_cube = cube
        if cube.var_name == 'pctisccp':
            meanctp_cube = cube
        if cube.var_name == 'albisccp':
            meanalbedo_cube = cube
        if cube.var_name == 'cltisccp':
            meantcc_cube = cube
        if cube.var_name == 'boxptopisccp':
            subcol_ctp_cube = cube
    
    subcol_ctp_cube.data.mask = np.where(subcol_ctp_cube.data < 0, 1.0, 0.0)
    
    subcol_ctp_mean_cube = subcol_ctp_cube.collapsed('subcolumn indices',iris.analysis.MEAN)
    
    alt_albedo_cube = meanalbedo_cube.copy()
    alt_albedo_cube = get_alt_albedo(jointpdf_cube,alt_albedo_cube)
    tropical_cubes = iris.cube.CubeList([])
    tropical_ctp  = iris.cube.CubeList([])
    tropical_albedo= iris.cube.CubeList([])
    tropical_tcc= iris.cube.CubeList([])
    tropical_tcc_alt = iris.cube.CubeList({})
    tropical_ctp_alt = iris.cube.CubeList({})
    tropical_albedo_alt = iris.cube.CubeList({})

    
    count=0
    for i, lat in enumerate(latitude_cube.data):
        if -20.0 < lat < 20.0:
            count=count+1
            tropical_cubes.append(jointpdf_cube[:,:,i])
            tropical_ctp.append(meanctp_cube[i])
            tropical_albedo.append(meanalbedo_cube[i])
            tropical_tcc.append(meantcc_cube[i])
            alternative_tcc = jointpdf_cube[:,:,i].collapsed(['air_pressure','atmosphere_optical_thickness_due_to_cloud'],iris.analysis.SUM)
            tropical_tcc_alt.append(alternative_tcc)
            tropical_ctp_alt.append(subcol_ctp_mean_cube[i])
            tropical_albedo_alt.append(alt_albedo_cube[i])

    jointpdf_cube = tropical_cubes.merge_cube()
    mean_ctp_cube = tropical_ctp.merge_cube()
    mean_ctp_cube.data = mean_ctp_cube.data / 100000.
    mean_tcc_cube = tropical_tcc.merge_cube()
    mean_tcc_cube.data = mean_tcc_cube.data / 100.
    mean_albedo_cube = tropical_albedo.merge_cube()
    # here we are getting an alternative tcc. Where we add up all the fractions
    # in the histogram.  The one that we took directly from the file
    # did not look correct.  The sum from the histogram looks better.  
    mean_tcc_alt_cube = tropical_tcc_alt.merge_cube()
    mean_tcc_alt_cube.data = mean_tcc_alt_cube.data / 100.
    mean_ctp_alt_cube = tropical_ctp_alt.merge_cube()
    mean_ctp_alt_cube.data = mean_ctp_alt_cube.data / 100000.
    mean_albedo_alt_cube = tropical_albedo_alt.merge_cube()
 
    

   
    mean_ctp_cube.data = np.where(mean_ctp_cube.data < 0, 
                                  np.nan, mean_ctp_cube.data)
    mean_tcc_cube.data = np.where(mean_tcc_cube.data < 0, 
                                  np.nan, mean_tcc_cube.data)
    mean_albedo_cube.data = np.where(mean_albedo_cube.data < 0, 
                                  np.nan, mean_albedo_cube.data)
    cloudtop_press = jointpdf_cube.coord('air_pressure').points
    thickname = 'atmosphere_optical_thickness_due_to_cloud'
    optical_depth = jointpdf_cube.coord(thickname).points
    dummy = [' ',' ',' ',' ',' ',' ',' ']
    
    #plt.plot(mean_tcc_cube.data)
    #plt.plot(mean_albedo_cube.data)
    #plt.plot(mean_ctp_cube.data)
    #plt.xlim(-5.0,100.0)
    #plt.vlines(0, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.vlines(10, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.vlines(20, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.vlines(30, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.vlines(40, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.vlines(50, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.vlines(60, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.vlines(70, ymin=0,ymax=1, linestyle='dotted',color='red')
    #plt.show()
    #sys.exit(0)
    fig = plt.figure(figsize=[11.4,11.4])
  
    for k in range(0, 90,10):
        pdf = jointpdf_cube.data[k,:,:]
        ax=plt.subplot(3,3,np.int(np.floor(k/10)+1))
        cs = ax.pcolormesh(pdf)
        if k==0 or k==30 or k==60 or k==150:
            plt.yticks([1,2,3,4,5,6,7], cloudtop_press/100., size=7)
            plt.ylabel('cloud top pressure')
        else:
            ax.set_yticklabels(dummy)
        if k>=60:
            plt.xticks([1,2,3,4,5,6,7], optical_depth, size=7,rotation=90)
            plt.xlabel('cloud optical depth')
        else:
            ax.set_xticklabels(dummy)
        plt.colorbar(cs)
        #plt.title(np.str(k) + ': ctp='+ "{:.2f}".format(mean_ctp_cube.data[k]),fontsize=7) 
        plt.title('alb='+"{:.2f}".format(mean_albedo_alt_cube.data[k]) + ': ctp='+ "{:.2f}".format(mean_ctp_alt_cube.data[k])+"\n" + 'tcc=' + "{:.2f}".format(mean_tcc_alt_cube.data[k]),fontsize=8)
    plt.subplots_adjust()
    plt.savefig('tropics_sample20.eps')
    plt.savefig('tropics_sample20.png')
    plt.close()

    return (jointpdf_cube.data,cloudtop_press, optical_depth,
            mean_ctp_alt_cube, mean_tcc_alt_cube, mean_albedo_alt_cube)
################################################################
def get_cloud_regiemes(joint_pdf_data, NK, cloud_regiemes, ctp, tau):
    """
    We are subsetting the joint_pdf_data into NK cloud regiemes
    We have some cloud_regiemes to start off with in cloud_regiemes
     
    For each pdf in joint_pdf_data
    1. decide which cloud regieme is nearest to this pdf
    2. add the pdf to this cloud regieme
    3. recalculate the cloud regieme with the addition of this pdf
    """
    fig = plt.figure(figsize=[11.4,11.4])
    dummy = [' ',' ',' ',' ',' ',' ',' ']
    nz, ny, nx = np.shape(joint_pdf_data)
    for k in range(0, 200,40):
        pdf = joint_pdf_data[:,:,k]
        ax=plt.subplot(4,5,k+1)
        cs = ax.pcolormesh(pdf)
        if k >=15:
            plt.xticks([1,2,3,4,5,6,7], tau, size=7, rotation=90)
        else:
            ax.set_xticklabels(dummy)
        if k==0 or k==5 or k==10 or k==15:
            plt.yticks([1,2,3,4,5,6,7], ctp, size=7)
        else:
            ax.set_yticklabels(dummy)
        plt.colorbar(cs)
        plt.title(np.str(k))
    plt.savefig('first20_histograms.eps')
    plt.close()


def test_sklearn(joint_pdf_data, tau, ctp):
    """
    testing stuff from the tutorial
    """
    nclusters=8
    nx, ny, nz = joint_pdf_data.shape
    joint_pdf_reshape = joint_pdf_data[:, :, :].reshape((nx * ny, nz))
    joint_pdf = np.transpose(joint_pdf_reshape)
    print(joint_pdf_reshape.shape)
    plt.subplot(2,1,1)
    for i in range(0,20):
        plt.plot(joint_pdf_reshape[:,i])

    kmeans = KMeans(init="random", n_clusters=nclusters, n_init=10)
    which_cluster=kmeans.fit_predict(joint_pdf)
    print(which_cluster[0:20])

    print(kmeans.cluster_centers_.shape)
    cluster_centers = kmeans.cluster_centers_.reshape((nclusters,nx, ny))
    plt.subplot(2,1,2)
    for i in range(0,nclusters):
        plt.plot(kmeans.cluster_centers_[i,:])
    plt.savefig('2d_clusters.eps')
    plt.close()
  
    fig = plt.figure(figsize=[11.4,11.4])
    print('cluster centers shape', cluster_centers.shape)

    for i in range(0,nclusters):
        ax = plt.subplot(3,3,i+1)
        cs = ax.pcolormesh(cluster_centers[i, :, :])
        plt.xticks([1,2,3,4,5,6,7], tau, size=7, rotation=90)
        plt.yticks([1,2,3,4,5,6,7], ctp, size=7)
        plt.colorbar(cs)
    plt.savefig('test_clusters.eps')
    plt.savefig('test_clusters.png')
    print('kmeans niter',kmeans.n_iter_)
    plt.close()
    # elbow method to see how many clusters
    #see=[]
    #for k in range(1,31):
    #    kmeans = KMeans(n_clusters=k)
    #    kmeans.fit(joint_pdf)
    #    see.append(kmeans.inertia_)
    #plt.plot(range(1,31),see)
    #plt.show()

       

def test_sklearn_tropics(joint_pdf_data, tau, ctp):
    """
    testing stuff from the tutorial
    """


    nclusters=9
    nx, ny, nz = joint_pdf_data.shape
    joint_pdf = joint_pdf_data[:, :, :].reshape((nx,nz * ny))
    print(joint_pdf.shape)
    
    scalar = StandardScaler()
    scale = scalar.fit_transform(joint_pdf)
    joint_pdf = scale
    print(joint_pdf_data)
    plt.subplot(2,1,1)
    for i in range(0,20):
        plt.plot(joint_pdf[i,:])

    kmeans = KMeans(init="k-means++", n_clusters=nclusters, n_init=10)
    which_cluster=kmeans.fit_predict(joint_pdf)
    print(which_cluster)
    total_in_cluster = np.zeros(nclusters)
    largest_cluster=0
    for i in range(0,nclusters):
        total_in_cluster[i] = np.sum((np.where(which_cluster == i, 1.0, 0.0)))
        if total_in_cluster[i] == np.max(total_in_cluster):
            largest_cluster=i

    print(kmeans.cluster_centers_.shape)
    cluster_centers = kmeans.cluster_centers_.reshape((nclusters,nz, ny))
    #plt.subplot(2,1,2)
    #for i in range(0,nclusters):
    #    plt.plot(kmeans.cluster_centers_[i,:])
    #plt.savefig('2d_clusters.eps')
    #plt.close()
  
    fig = plt.figure(figsize=[11.4,11.4])
    print('cluster centers shape', cluster_centers.shape)

    for i in range(0,nclusters):
        ax = plt.subplot(3,3,i+1)
        cs = ax.pcolormesh(cluster_centers[i, :, :])
        plt.xticks([1,2,3,4,5,6,7], tau, size=7, rotation=90)
        plt.yticks([1,2,3,4,5,6,7], ctp, size=7)
        plt.title("{:.1f}".format(total_in_cluster[i] *100./ np.sum(total_in_cluster)) + '%')
        plt.colorbar(cs)
    plt.savefig('test_clusters.eps')
    plt.savefig('test_clusters.png')
    print('kmeans niter',kmeans.n_iter_)
    print('largest cluster',largest_cluster)
    plt.close()
    # elbow method to see how many clusters
    #see=[]
    #for k in range(1,31):
    #    kmeans = KMeans(n_clusters=k)
    #    kmeans.fit(joint_pdf)
    #    see.append(kmeans.inertia_)
    #plt.plot(range(1,31),see)
    #plt.show()

 ##########################################################
def tropical_cloud_regiemes(CTP_data, TCC_data, albedo_data,
                            joint_pdf_data, ctpall, tauall):
    """
    estimates the cloud regieme for each tropical gridpoint
    """
    REGIEME_NAMES_TR = {0:"Shallow cumulus", 1:"Congestus",
                 2:"Thin cirrus", 3:"Stratocu./Cu. Transition",
                 4:"Anvil cirrus", 5:"Deep Convection",
                        6:"Stratocumulus", 7:"ClearSky"}
    REG_CHARS = {0:[0.261, 0.652, 0.314], 
                 1:[0.339, 0.483, 0.813],
                 2:[0.211, 0.356, 0.740],
                 3:[0.338, 0.784, 0.640],
                 4:[0.313, 0.327,0.944],
                 5:[0.532, 0.285,0.979],
                 6:[0.446, 0.722, 0.824]}

    regieme=np.zeros(len(CTP_data))
    min_dists=np.zeros(len(CTP_data))
    nz, ny, nx = np.shape(joint_pdf_data)
    regieme_pdf = np.zeros((8, ny, nx))
    regieme_count = np.zeros((8))
    dists = np.zeros(7)
    for i,ctp in enumerate(CTP_data):
        features = np.array([albedo_data[i], ctp, TCC_data[i]])
        if TCC_data[i] < 0.05:
            regieme[i]=7
            regieme_pdf[7,:,:] = joint_pdf_data[i,:,:]
            regieme_count[7] = regieme_count[7] + 1.0
        else:
            for reg in range(0,7):
                dists[reg] = np.linalg.norm(features - REG_CHARS.get(reg))
            regieme[i] = np.argmin(dists)
            regieme_pdf[int(regieme[i]), :, :] = joint_pdf_data[i,:,:]
            regieme_count[int(regieme[i])]=regieme_count[int(regieme[i])] + 1.0
            min_dists[i] = dists[np.argmin(dists)]
            if i==0 or i==10 or i==20 or i==30 or i==40 or i==50 or i==170:
                print(i,REGIEME_NAMES_TR.get(regieme[i]),dists,features)

    for reg in range(0,8):
        print(reg,np.sum(np.where(regieme==reg, 1.0, 0.0))*100./275.,regieme_count[reg]*100./275.)

    # create average regieme for each type and plot
    for reg in range(0,8):
        regieme_pdf[reg,:,:] = regieme_pdf[reg,:,:] / regieme_count[reg]
    
    fig = plt.figure(figsize=[11.4,11.4])
   
    for i in range(0,8):
        ax = plt.subplot(3,3,i+1)
        cs = ax.pcolormesh(regieme_pdf[i, :, :])
        plt.xticks([1,2,3,4,5,6,7], tauall, size=7, rotation=90)
        plt.yticks([1,2,3,4,5,6,7], ctpall, size=7)
        plt.title(REGIEME_NAMES_TR.get(i)+np.str(regieme_count[i]))
        plt.colorbar(cs)
    plt.savefig('williams_webb_tropical_clusters.eps')
    
    
########################################################
# main program

FILEIN = '/nfs/see-fs-02_users/earjcti/PYTHON/PROGRAMS/COLLABORATORS/TAMARA/CFMIP/COSPv2.0_juliatest/driver/data/outputs/UKMO/cosp2_output_um.nc'

# note joint pdf data(7, 7, 1236) dim3=loc/time, dim2=optical depth, dim1=ctp
#(joint_pdf_data, cloudtop_press, optical_depth) = get_pdfcube()

# subset the data into NK cloud_regiemes the start cloud regiemes will 
# be the first NK locations of joint_pdf_data
#NK=2
#start_cloud_regiemes = joint_pdf_data[0:NK, :, :]
#cloud_regiemes = get_cloud_regiemes(joint_pdf_data, NK, start_cloud_regiemes,
#                                    cloudtop_press,optical_depth)

print('julia0')
(joint_pdf_data, cloudtop_press, optical_depth, mean_ctp_cube, 
mean_tcc_cube, mean_albedo_cube) = get_pdfcube_tropics()

# attempt to put mean_ctp,mean_tcc and mean_albedo onto the cloud regieme
# averages
print('juli1')
cloud_regiemes = tropical_cloud_regiemes(mean_ctp_cube.data, mean_tcc_cube.data, mean_albedo_cube.data, joint_pdf_data, cloudtop_press, optical_depth)

print('julia2')
#plt.subplot(311)
#plt.plot(mean_ctp_cube.data)
#plt.subplot(312)
#plt.plot(mean_tcc_cube.data)
#plt.subplot(313)
#plt.plot(mean_albedo_cube.data)
#plt.show()

test_sklearn_tropics(joint_pdf_data, optical_depth, cloudtop_press)

