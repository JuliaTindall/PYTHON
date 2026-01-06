
#!/usr/bin/env python2.7
#NAME
#    PLOT_MOC
#PURPOSE
#    This program will do some plots on the MOC (which were calculated in the 
#    ~earjcti/MOC directory
#
# search for 'main program' to find end of functions
# Julia 14/1/2017



import os
import numpy as np
import scipy as sp
import scipy.signal as sig
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset, MFDataset
import sys
from mpl_toolkits.basemap import Basemap, shiftgrid
import sklearn.decomposition as sk
from sklearn.preprocessing import normalize



#functions are:
#  def plotdata
#  def annmean
#  def seasmean

# functions start here
def plotdata(plotdata,fileno,lat,dep,titlename,minval,maxval,valinc,V,uselog,cbarname,pcplot):
    lats, depths = np.meshgrid(lat,dep)

    if pcplot == 'avg':
        plt.subplot(1,1,1)
    else:
        if pcplot == 'pc':
            plt.subplot(2,2,fileno+1)
        else:
            plt.subplot(5,2,fileno+1)

    V=np.arange(minval,maxval,valinc)
    if pcplot =='pc':
        cs = plt.contourf(lats,depths,plotdata,V,cmap='RdBu_r',extend='both')
    else:
        cs = plt.contourf(lats,depths,plotdata,V,extend='both')
    plt.gca().invert_yaxis()
    plt.title(titlename)

    if pcplot != 'n' or fileno >9:
        cbar = plt.colorbar(cs,orientation="horizontal")
        cbar.set_label(cbarname,labelpad=-40)



#end def plotdata

def plotmap(plotdata,fileno,lon,lat,titlename,minval,maxval,valinc,V,uselog,cbarname):
    lons, lats = np.meshgrid(lon,lat)
    plt.subplot(2,2,fileno+1)

   # this is good for a NAO region
   # map=Basemap(width=12000000,height=8000000,projection='stere',\
   #                 resolution='c',lat_ts=50,lat_0=50,lon_0=0)
   # this is good for the globe
    map=Basemap(llcrnrlon=-180.0,urcrnrlon=180.0,llcrnrlat=-90.0,urcrnrlat=90.0,projection='cyl',resolution='c')
    #map.drawmapboundary(fill_color='aqua')
    map.drawmapboundary
    x, y = map(lons, lats)
    map.drawcoastlines()
    if V == 0:
        V=np.arange(minval,maxval,valinc)
    if uselog =='y':
        cs = map.contourf(x,y,plotdata,V,norm=mp.colors.PowerNorm(gamma=1./3.))
        cbar = plt.colorbar(cs,orientation="horizontal",extend='both')
    else:
        if uselog =='la':
            cs = map.contourf(x,y,plotdata,V,norm=mp.colors.SymLogNorm(linthresh=2.0,linscale=2.0,vmin=-32,vmax=32),cmap='RdBu')
            cbar = plt.colorbar(cs,orientation="horizontal",extend='both')

        else:
            if uselog =='a':
                cs = map.contourf(x,y,plotdata,V,cmap='RdBu',extend='both')
                cbar = plt.colorbar(cs,orientation="horizontal")
            else:
                if uselog =='at':
                    cs = map.contourf(x,y,plotdata,V,cmap='RdBu_r',extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")
                else:
                    print(np.shape(x),np.shape(y),np.shape(plotdata))
                    cs = map.contourf(x,y,plotdata,V,extend='both')
                    cbar = plt.colorbar(cs,orientation="horizontal")

    plt.title(titlename)
    cbar.set_label(cbarname,labelpad=-40)
#end def plotmap


def indexplot(toplot,xdata,fileno,data_sm,xmin,xmax,expt,xlabel,ymin,ymax,titlename):
    plt.subplot(4,2,fileno*2)
    print(fileno)

    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    datasize=len(toplot)

    
    # plot data
    plt.plot(xdata,toplot)
    plt.title(titlename)
    plt.xlabel(xlabel)
    # overplot smoothed data
    plt.plot(xdata,data_sm,'-')
    # overplot zero line and +-0.5deg line
    plt.plot(xdata,np.zeros(datasize))
    #bar_width=1.0/12.0
    #plt.bar(xdata,elninoarr,bar_width,color='red',edgecolor="none")
    #plt.bar(xdata,laninaarr,bar_width,color='blue',edgecolor="none")
   
# 

# end def indexplot

##################################################
def plot_avg_moc(expt_name,extra,yearstart,yearend):

    nyears=yearend-yearstart+1
    figcount=0
    for year in range(yearstart,yearend):
        if year >= 10:
            datasetname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pk2/'+expt_name+'o@pg'+extra+str(year)+'c1.nc'
        else:
            datasetname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pk2/'+expt_name+'o@pg'+extra+'0'+str(year)+'c1.nc'
        f=Dataset(datasetname)
        lat = f.variables['latitude'][:]
        depth = f.variables['depth'][:]
        ndepth=len(depth)
        nlat=len(lat)
        AMOC=f.variables['Merid_Atlantic'][:]
        AMOC=np.squeeze(AMOC)
        if year == yearstart:
            allAMOC=np.zeros((nyears,ndepth,nlat))

        allAMOC[year-yearstart,:,:]=AMOC
        f.close()

    avgAMOC=np.mean(allAMOC,axis=0)
    print(np.shape(avgAMOC))
    titlename='AMOC avg: '+expt_name
    plotdata(avgAMOC,-99,lat,depth,titlename,-18,20,2.0,0.0,'n','Sv','avg')
        
    plt.tight_layout()
    plt.show()
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/'+expt_name+'_avg.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()
    
    return(lat,depth,avgAMOC)
    
    
#end def plot_avg_MOC
#======================================================
def amoc_diff(expt_name,cntl_name,AMOC_e,AMOC_c,lat,depth):
# this function will difference two AMOCs
    AMOC_diff=AMOC_e-AMOC_c
    titlename='diff '+expt_name+'  and '+cntl_name
    plotdata(AMOC_diff,-99,lat,depth,titlename,-4,5,1.0,0.0,'n','Sv','avg')
    plt.tight_layout()
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/'+expt_name+'-'+cntl_name+'_avg.eps' 
    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()
    
#end def amoc_diff
#==========================================
def plot_all_moc(expt_name,extra,yearstart,yearend,plotallyears):

    nyears=yearend-yearstart+1
    figcount=0
    for year in range(yearstart,yearend):
        if year >= 10:
            datasetname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pk2/'+expt_name+'o@pg'+extra+str(year)+'c1.nc'
        else:
            datasetname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pk2/'+expt_name+'o@pg'+extra+'0'+str(year)+'c1.nc'
        f=Dataset(datasetname)
        lat = f.variables['latitude'][:]
        depth = f.variables['depth'][:]
        ndepth=len(depth)
        nlat=len(lat)
        AMOC=f.variables['Merid_Atlantic'][:]
        AMOC=np.squeeze(AMOC)

        if plotallyears == 'y':
            titlename='AMOC yr:'+str(year)
            plotdata(AMOC,figcount,lat,depth,titlename,-30,30,5.0,0.0,'n','Sv')
            if figcount == 9:
                plt.tight_layout()
                plt.show()
                sys.exit()
            figcount=(figcount+1)%10
          
        if year == yearstart:
            allAMOC=np.zeros((nyears,ndepth,nlat))

        allAMOC[year-yearstart,:,:]=AMOC
        f.close()
    

    sys.exit()
    
    
#end def plot_all_MOC
#######################################
# do principal component analysis to find the periods over which the MOC varies
def MOC_PC(expt_name,extra,yearstart,yearend):

# read in data

    nyears=yearend-yearstart+1
    figcount=0
    for year in range(yearstart,yearend+1):
        if year >= 10:
            datasetname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pk2/'+expt_name+'o@pg'+extra+str(year)+'c1.nc'
        else:
            datasetname='/nfs/hera1/earjcti/um/HadGEM_data/'+expt_name+'/pk2/'+expt_name+'o@pg'+extra+'0'+str(year)+'c1.nc'
        f=Dataset(datasetname)
        lat = f.variables['latitude'][:]
        depth = f.variables['depth'][:]
        ndepth=len(depth)
        nlat=len(lat)
        AMOC=f.variables['Merid_Atlantic'][:]
        AMOC=np.squeeze(AMOC)

        # mask out region south of equator
        ix1=(lat >=0)
        lats_reg=lat[ix1]
        nlat=len(lats_reg)
        lat=lats_reg
        AMOC=AMOC[:,ix1]
    

        if year == yearstart:
            allAMOC=np.zeros((nyears,ndepth,nlat))
            
        print(yearend-yearstart,year-yearstart,nyears)
        allAMOC[year-yearstart,:,:]=AMOC
        f.close()


    

  
    # remove time average from all AMOC 
    for j in range(0,nlat):
        for i in range(0,ndepth):
            allAMOC[:,i,j]=allAMOC[:,i,j]-np.mean(allAMOC[:,i,j])
    
    # ideally we will need to multiply by a weighting factor to account for the fact that all latitudes have a different number of gridboxes at different sizes


    # reshape and transpose the data to the correct dimension
    print('amoc shape',np.shape(allAMOC),nlat,ndepth,nyears)
    rs_allAMOC_nt=np.reshape(allAMOC,(nyears,nlat*ndepth))
    rs_allAMOC=np.transpose(rs_allAMOC_nt)

    # do a PC analysis using sklearn
    neofs=2
    AMOCpca=sk.PCA(n_components=neofs)
    AMOCpca.fit(rs_allAMOC)
    expl_var=AMOCpca.explained_variance_ratio_
    EOFs=AMOCpca.transform(rs_allAMOC)
    
    # scale so each EOF has unit length
    EOFs=normalize(EOFs,axis=0)

    # get principal components
    PCs=np.mat(rs_allAMOC_nt) * np.mat(EOFs)

    for i in range(0,neofs):
        EOF_temp=EOFs[:,i]
        EOF_plot=np.reshape(EOF_temp,(ndepth,nlat))
        stdevpc=np.std(PCs[:,i])
        PCs[:,i]=PCs[:,i]/stdevpc
        EOF_plot=EOF_plot * stdevpc
        
        titlename='PC'+str(i+1)+':'+expt_name+' '+str(np.ceil(expl_var[i]*100.))+'%'
        plotdata(EOF_plot,(i*2),lat,depth,titlename,-2.0,2.2,0.2,0,'n','Sv','pc')
        
        toplot=PCs[:,i]
        datasize=len(toplot)
        xdata=np.arange(datasize)
        indexplot(toplot,xdata,(2*i)+1,toplot,0,nyears,expt_name,'year',-2.0,2.0,'index')

        # do a spectral analysis on the index
        # to see over which periods it is varying

        print('toplot is',toplot)
        toplot=np.squeeze(toplot)

        Pxx_f, Pxx_den=sig.periodogram(toplot,1.0)
        Pxx_den=np.squeeze(Pxx_den)
        print('PXx',Pxx_f,'size',np.shape(Pxx_f),np.shape(toplot))
        print('pxx2',Pxx_den,'size',np.shape(Pxx_den))
        indexplot(Pxx_den,Pxx_f,(2*i)+2,Pxx_den,0,0.5,expt_name,'cycles per year',0,10,'index power spectrun')
        

    plt.tight_layout()
  
    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/'+expt_name+'_PCs.eps' 
    plt.savefig(fileout, bbox_inches='tight')  
    plt.close()


#  write both indexes out to a file

    fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/'+expt_name+'_PCs.txt' 
    f1=open(fileout,'w+')
    f1.write('strength of principal components')
    f1.write('nyears='+str(nyears)+' \n')
    f1.write('extra year  PC1    PC2 \n')
    for year in range (yearstart,yearend+1):
        yearuse=year
        extrause=extra
        if yearuse >= 100:
            yearuse=yearuse-100
            extrause=chr(ord(extra)+1)

        if yearuse < 0:
            yearuse=yearuse+100
            extrause=chr(ord(extra)-1)

        yearfname=str(yearuse)

        if yearuse < 10:
            yearfname='0'+str(yearuse)

        print(year,yearstart,yearend,np.shape(PCs))
        f1.write(extrause+';'+str(yearfname)+';'+str(PCs[year-yearstart,0])+';'+str(PCs[year-yearstart,1])+'\n')
    f1.close()
        


# 

#end def annmean

# relate principal components to climate
def PCs_to_climate_telecon(exptname,fieldlocation,fileext,fieldname,seasname,monthnames,lonname,latname):
 
# this will plot the teleconnections associated with the MOC principal 
# components by taking the
# most exteme 5% of the MOC years and plotting the climate anomaly


    # read in the MOC principal components

    filein='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/'+exptname+'_PCs.txt' 
    f=open(filein,'r')
    # get titleline which contains number of years
    textline=f.readline()
    b=textline.split()  # split text by removing newline
    print(b)
    c=b[3]
    b=c.split('=')  # split text by removing equals sign
    nyears=int(b[1])
    # discard second titleline
    textline=f.readline()

   
    # read over the rest of the data
    extraindex=np.empty(nyears,dtype=np.dtype('S1'))
    yearindex=np.zeros(nyears)
    PC1index=np.zeros(nyears)
    PC2index=np.zeros(nyears)
    extremePC1=np.zeros(nyears)  # here we mark the 5% most extreme values
    extremePC2=np.zeros(nyears)  # here we mark the 5% most extreme values
    

    count=0
    for line in f:
        # extract extra year and NAOindex
        linesplit=line.split(';')   # the data in the file is split by ;
        extraindex[count]=linesplit[0]
        yearindex[count]=linesplit[1]
        PC1index[count]=linesplit[2]
        PC2index[count]=linesplit[3]
        count=count+1


    # next we want to find the years that have the largest positive and negative
    # nao index

    num_extr=int(np.ceil(nyears*0.05))

    # get highest 5 values.
    # note that 'zip' zips the arrays together forming a multi dim list
    # sort sorts the list on the first element
    # reverse will reverse the sort
    lowPC1=sorted(zip(PC1index,extraindex,yearindex))[:num_extr]
    uppPC1=sorted(zip(PC1index,extraindex,yearindex),reverse=True)[:num_extr]
    lowPC2=sorted(zip(PC2index,extraindex,yearindex))[:num_extr]
    uppPC2=sorted(zip(PC2index,extraindex,yearindex),reverse=True)[:num_extr]
    

    # we will now put our indices to one side and get the data we are 
    # interested in 
    # the field we are using will be passed in the calling program

    nmonths=len(monthnames)

    dirname='/nfs/hera1/earjcti/um/HadGEM_data/'+exptname+'/'+fieldlocation+'/'
    os.chdir(dirname)

    # firstly get average field over the season
    allfiles=[]
    for month in monthnames:
        allfiles.append(dirname+exptname+'a@pd*'+month+fileext)

    for monthno in range(0,nmonths):
        print(allfiles)
        f=MFDataset(allfiles[monthno])
        lat = f.variables[latname][:]
        lon = f.variables[lonname][:]
    
        if len(fieldname) ==1 :
            atemp=f.variables[fieldname[0]][:]
            atemp=np.squeeze(atemp)

        if len(fieldname) == 2:
            atemp=f.variables[fieldname[0]][:]
            atemp=np.squeeze(atemp)
            atemp2=f.variables[fieldname[1]][:]
            atemp2=np.squeeze(atemp2)
        
        if len(fieldname) > 2:
            print('length of fieldname is', len(fieldname))
            print('you are requesting too many variables')
            sys.exit()
        f.close()
        ntimes,ny,nx=np.shape(atemp)

        #average across the time dimension
        temp_m1=np.mean(atemp,axis=0)

        # set array for storing average
        if monthno == 0:
            temp_avg=temp_m1
        else:
            temp_avg=temp_avg+temp_m1

        
        if len(fieldname) ==2:
            temp_m1=np.mean(atemp2,axis=0)
            if monthno == 0:
                temp2_avg=temp_m1
            else:
                temp2_avg=temp2_avg+temp_m1

                

    temp_avg=temp_avg/nmonths
    if len(fieldname) == 2:
        temp2_avg=temp2_avg/nmonths

    

    # now get data for the highest and lowest years


    for ex in range(0,2):     # loop for lowest or highest year
        for npcs in range (0,2):   # we will have to loop over PCS
            if ex == 0:
                if npcs == 0: 
                    extremedata=lowPC1
                if npcs == 1:
                    extremedata=lowPC2
            if ex == 1:
                if npcs == 0:
                    extremedata=uppPC1
                if npcs == 1:
                    extremedata=uppPC2


            for time in range (0,num_extr):
                for monthno in range (0,nmonths):
                    singleline=extremedata[time]
                    yearuse=int(singleline[2])
                    extrause=singleline[1]
                    month=monthnames[monthno]

                    if seasname == 'djf-1' and month == 'dc':
                        yearuse=yearuse-1

                    if seasname == 'djf+1' and month == 'ja':
                        yearuse=yearuse+1

                    if seasname == 'djf+1' and month == 'fb':
                        yearuse=yearuse+1

                    if yearuse >= 100:
                        yearuse=yearuse-100
                        extrause=chr(ord(extrause)+1)
                                    
                    if yearuse < 0:
                        yearuse=yearuse+100
                        extrause=chr(ord(extrause)-1)
                                        
                    yearfname=str(yearuse)

                    if yearuse < 10:
                        yearfname='0'+str(yearuse)
                        

                    fname=dirname+exptname+'a@pd'+extrause+yearfname+month+fileext
                
                    f=Dataset(fname,mode='r')
                    lat = f.variables[latname][:]
                    latsize=len(lat)
                    lon = f.variables[lonname][:]
                    lonsize=len(lon)
                    lontemp=lon

                    if len(fieldname) ==1 :
                        atemp=f.variables[fieldname[0]][:]
                        atemp=np.squeeze(atemp)

                    if len(fieldname) == 2:
                        atemp=f.variables[fieldname[0]][:]
                        atemp=np.squeeze(atemp)
                        atemp2=f.variables[fieldname[1]][:]
                        atemp2=np.squeeze(atemp2)
                        
                    if len(fieldname) > 2:
                        print('you are requesting too many variables')
                        sys.exit()
    

                    # set array for storing average
                    if monthno == 0 and time == 0:
                        temp_extreme=atemp
                        count=1
                        count2=1
                        temp2_extreme=0.
                        if len(fieldname)==2:
                            temp2_extreme=atemp2
                    else:
                        temp_extreme=temp_extreme+atemp
                        count=count+1
                        if len(fieldname)==2:
                            temp2_extreme=temp2_extreme+atemp2
                            count2=count2+1

            # put temperature data in lower or higher catogry
            if ex == 0:
                if npcs == 0: 
                    temp_low_PC1=temp_extreme/count
                    temp2_low_PC1=temp2_extreme/count2
                if npcs == 1: 
                    temp_low_PC2=temp_extreme/count
                    temp2_low_PC2=temp2_extreme/count2
            if ex == 1:
                if npcs == 0: 
                    temp_high_PC1=temp_extreme/count
                    temp2_high_PC1=temp2_extreme/count2
                if npcs == 1: 
                    temp_high_PC2=temp_extreme/count
                    temp2_high_PC2=temp2_extreme/count2


    # we have finished with the loop        
    # shiftdata for plot
    lontemp=lon
    temp_low_PC1,lon = shiftgrid(180.,temp_low_PC1,lon,start=False)    
    lon=lontemp
    temp_high_PC1,lon = shiftgrid(180.,temp_high_PC1,lon,start=False)    
    lon=lontemp
    temp_low_PC2,lon = shiftgrid(180.,temp_low_PC2,lon,start=False)    
    lon=lontemp
    temp_high_PC2,lon = shiftgrid(180.,temp_high_PC2,lon,start=False)    
    lon=lontemp
    temp_avg,lon = shiftgrid(180.,temp_avg,lon,start=False)    



    if fieldname[0] == 'temp_1':  # temperature data
        plotmap(temp_low_PC1-temp_avg,0,lon,lat,'low PC1 Tanom',-2.0,2.2,0.2,0,'at','degC')
        plotmap(temp_high_PC1-temp_avg,1,lon,lat,'high PC1 Tanom',-2.0,2.2,0.2,0,'at','degC')
        plotmap(temp_low_PC2-temp_avg,2,lon,lat,'low PC2 Tanom',-2.0,2.2,0.2,0,'at','degC')
        plotmap(temp_high_PC2-temp_avg,3,lon,lat,'high PC2 Tanom',-2.2,2.0,0.2,0,'at','degC')
       

        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_PCs_tempanom'+exptname+'_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        plt.show()
        plt.close()


        titlename=exptname+' low seas temp'
        plotmap(temp_low_PC1-273.15,0,lon,lat,'low PC1 temp',-40.0,40.0,10.0,0,'n','degC')
        plotmap(temp_high_PC1-273.15,1,lon,lat,'high PC1 temp',-40.0,40.0,10.0,0,'n','degC')
        plotmap(temp_low_PC2-273.15,2,lon,lat,'low PC2 temp',-40.0,40.0,10.0,0,'n','degC')
        plotmap(temp_high_PC2-273.15,3,lon,lat,'high PC2 temp',-40.0,40.0,10.0,0,'n','degC')

        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_PCs_temp'+exptname+'_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        plt.close()

    


      


    if fieldname[0] == 'p':  # mean se level pressure data
        plotmap((temp_low_PC1-temp_avg)/100.,0,lon,lat,'low PC1 mslp anom',-5.,6.,1.,0,'at','hPa')
        plotmap((temp_high_PC1-temp_avg)/100.,1,lon,lat,'high PC1 mslp anom',-5.,6.,1.,0,'at','hPa')
        plotmap((temp_low_PC2-temp_avg)/100.,2,lon,lat,'low PC2 mslp anom',-5.,6.,1.,0,'at','hPa')
        plotmap((temp_high_PC2-temp_avg)/100.,3,lon,lat,'high PC2 mslp anom',-5.,6.,1.,0,'at','hPa')
       

        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_PCs_mslpanom'+exptname+'_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight')  
        plt.close()

        plotmap(temp_low_PC1/100.,0,lon,lat,'low PC1 mslp',980.,1040.,5,0,'n','hPa')
        plotmap(temp_high_PC1/100.,1,lon,lat,'high PC1 mslp',980.,1040.,5,0,'n','hPa')
        plotmap(temp_low_PC2/100.,2,lon,lat,'low PC2 mslp',980.,1040.,5,0,'n','hPa')
        plotmap(temp_high_PC2/100.,3,lon,lat,'high PC2 mslp',980.,1040.,5,0,'n','hPa')


        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_PCs_mslp_'+exptname+'_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight') 
        plt.show()
        plt.close()

    


       

    if fieldname[0] == 'precip_1':  # precipitation data
        titlename=exptname+' low seas precip'
        plotmap(temp_low_PC1*60.*60.*24.*30.,0,lon,lat,'low PC1 precip',-0,275,25,0,'n','mm/month')
        plotmap(temp_high_PC1*60.*60.*24.*30.,1,lon,lat,'high PC1 precip',0,275,25,0,'n','mm/month')
        plotmap(temp_low_PC2*60.*60.*24.*30.,2,lon,lat,'low PC2 precip',-0,275,25,0,'n','mm/month')
        plotmap(temp_high_PC2*60.*60.*24.*30.,3,lon,lat,'high PC2 precip',0,275,25,0,'n','mm/month')

        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_precip'+exptname+'_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight')  

        plt.close()


        plotmap((temp_low_PC1-temp_avg)*60.*60.*24.*30.,0,lon,lat,'low PC1 panom',-30,35,5,0,'a','mm/month')
        plotmap((temp_high_PC1-temp_avg)*60.*60.*24.*30.,1,lon,lat,'high PC1 panom',-30,35,5,0,'a','mm/month')

        plotmap((temp_low_PC2-temp_avg)*60.*60.*24.*30.,2,lon,lat,'low PC2 panom',-30,35,5,0,'a','mm/month')
        plotmap((temp_high_PC2-temp_avg)*60.*60.*24.*30.,3,lon,lat,'high PC2 panom',-30,35,5,0,'a','mm/month')

     
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_precipanom'+exptname+'_'+seasname+'.eps' 
        plt.savefig(fileout, bbox_inches='tight')  

        plt.close()

        # if precip also do percentage change
        pcent_low_PC1=((temp_low_PC1-temp_avg)/temp_avg)*100.
        pcent_high_PC1=((temp_high_PC1-temp_avg)/temp_avg)*100.
        pcent_low_PC2=((temp_low_PC2-temp_avg)/temp_avg)*100.
        pcent_high_PC2=((temp_high_PC2-temp_avg)/temp_avg)*100.
        titlename=exptname+' low PC1 precip anom'
        plotmap(pcent_low_PC1,0,lon,lat,titlename,-50,60,5,0,'a','%')
        plotmap(pcent_high_PC1,1,lon,lat,'high PC1 Precip anomaly',-50,60,5,0,'a','%')
        plotmap(pcent_low_PC2,2,lon,lat,'low PC2 precip anomaly',-50,60,5,0,'a','%')
        plotmap(pcent_high_PC2,3,lon,lat,'high PC2 Precip anomaly',-50,60,5,0,'a','%')
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_precip_pcent'+exptname+'_'+seasname+'.eps' 




    if len(fieldname)==2 and fieldname[0] == 'u':  # winds dat
        ms='m/s'
        lon=lontemp
        temp2_low,lon = shiftgrid(180.,temp2_low,lon,start=False)    
        lon=lontemp
        temp2_high,lon = shiftgrid(180.,temp2_high,lon,start=False)    
        lon=lontemp
        temp2_avg,lon = shiftgrid(180.,temp2_avg,lon,start=False)    

        titlename=exptname+' low seas winds'
        plotquiver(temp_low,temp2_low,lon,lat,0,titlename,0,400,40.0,0.0,'n',ms)
        plotquiver(temp_high,temp2_high,lon,lat,1,'high seas winds',1,400,40.0,0.0,'n',ms)       
        #plotquiver(temp_avg,temp2_avg,lon,lat,1,'avg winds',1,400,40.0,0.0,'n',ms)       
        plotquiver(temp_low-temp_avg,temp2_low-temp2_avg,lon,lat,2,'low seas uvanom',0,400,40.0,0.0,'n',ms)
        plotquiver(temp_high-temp_avg,temp2_high-temp2_avg,lon,lat,3,'high seas uvanom',1,400,40.0,0.0,'n',ms)       
     
        fileout='/nfs/see-fs-02_users/earjcti/PYTHON/PLOTS/HadGEM2/plot_MOC/teleconnections_PC/tele_winds'+exptname+'_'+seasname+'.eps' 

    


    plt.savefig(fileout, bbox_inches='tight')  

    plt.close()
 
#end def PCs_to_climate

################################
# main program


##################################
# plot average of the MOC

#retdata=plot_avg_moc('xkvje','n',1,99)
#lat=retdata[0]
#depth=retdata[1]
#xkvje_AMOC=retdata[2]

#retdata=plot_avg_moc('xkvjf','n',1,98)
#xkvjf_AMOC=retdata[2]

#retdata=xkvjg_AMOC=plot_avg_moc('xkvjg','n',1,99)
#xkvjg_AMOC=retdata[2]

################################################
# difference of two AMOCs
#amoc_diff('xkvjg','xkvje',xkvjg_AMOC,xkvje_AMOC,lat,depth)
#amoc_diff('xkvjf','xkvje',xkvjf_AMOC,xkvje_AMOC,lat,depth)


##################################################
# plot all years of the MOC

#plt.figure(figureno)
#plotallyears='y'
#plot_all_moc('xkvjg','n',1,99,plotallyears)
#figureno=figureno+1

###########################################################################
# Principal component analysis to determine periods of variability in the MOC

#MOC_PC('xkvje','n',1,98)
#MOC_PC('xkvjg','n',1,98)

#############################################################
# Relate extreme values of PC's to climate


# usage Pcs_to_climate(expt_name',dirname,fileextension,fieldname,monthnames)
PCs_to_climate_telecon('xkvjg','temp_data','_temp.nc',['temp_1'],'djf+1',\
            ['dc','ja','fb'],\
            'longitude','latitude')

PCs_to_climate_telecon('xkvjg','precip_data','_precip.nc',['precip_1'],'djf+1',\
            ['dc','ja','fb'],\
            'longitude','latitude')

PCs_to_climate_telecon('xkvjg','mslp_data','_mslp.nc',['p'],'djf+1',\
            ['dc','ja','fb'],\
            'longitude','latitude')





sys.exit(0)

####

