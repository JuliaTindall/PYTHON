import numpy as np



def readproxy():
    """
    reads in the excel spreadsheet of the proxy dataset
    returns arrays of the latitude and longitude and sstanom of the proxy dataset
    """
    
    lats = np.zeros(37)
    lons = np.zeros(37)
    sstanom  = np.zeros(37)
    f1 = open(FILEPROXY,'r') # to read
    count=0
    for line in f1.readlines():
        if count >= 37:
            pass
        else:
            if line[0:4] == 'Loca': # titleline ignore
                print('titleline is',line)
                pass
            else:
                vals = line.split(',')
             
                if float(vals[2]) > 0:
                    lons[count] = float(vals[2])
                else:
                    lons[count] = float(vals[2]) + 360.
                lats[count] = float(vals[1])
                sstanom[count] = float(vals[15]) # this is plio - noaa
                count = count + 1
          
            
       
    f1.close()
   
    return lats, lons, sstanom




def get_model_output(lat,lon):
    """
    gets the model output for this lat and this long)


############################################
#
# get data


FILEPROXY = 'LP_data/cs_mp_sst_data_30k_plusNOAA.csv'

# some file stuff.
#self.filename = DATASTART + 'pliovar_uk37_ori_vs_bayspline.xlsx'
#            self.bsloc = 8
#        if datatype == 'MGCA':
#            self.filename = DATASTART +  'pliovar_mgca_OrivsBaymag.xlsx'
#            self.bsloc = 7
#        self.metafile = DATASTART + 'pliovar_metadata_global_02102019.csv'
#        sel


proxylat, proxylon, proxy_sst_anom = readproxy()

print(proxylat)
print(proxylon)
print(proxy_sst_anom)
print(len(proxylat))

#####################################################
#
# get model output

modelnames=['ACCESS', 'AWI-CM2', 
            'AWI-ESM2', 'COSMOS', 'HadCM3', 'MIROC4m','iCESM3']

all_model_output = []

for model in modelnames:
    model_points = get_model_output(proxylat,proxylon)
    all_model_output.append(model_points)

nmodels=len(modelnames_eoi400) # we also have HadISST in e280
