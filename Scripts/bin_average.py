import xarray as xr
import numpy as np
import os

def bin_av(datadir,filename):
    """
    Compute bin averge every 1m depth and overwrite the datafile
    """
    data = xr.open_dataset(os.path.join(datadir,filename))
    keys_tot = [i for i in data.keys() if i!='DEPTH']
    datanew_tot = {i:[] for i in keys_tot}
    #print(datanew_tot)
    depth_max = int(data.DEPTH.max()) # as integer
    #print(depth_max)
    for i in range(1,depth_max):
        for key in keys_tot:
            datanew_tot[key].append(float(data[key][(data.DEPTH < i+0.5) & \
                                              (data.DEPTH >= i-0.5)].mean().values))
    #print(datanew_tot)
    datanew = data.copy()
    datanew=datanew.rename({'scan':'scan_old'})
    #print(datanew)
    datanew['DEPTH'] = ('scan', np.arange(1,depth_max))
    for key in keys_tot:
        datanew[key] = ('scan', datanew_tot[key])
    #print(datanew)
    print('Bin averaging finished')
    print('overwriting data')
    os.rename(os.path.join(datadir,filename),os.path.join(datadir,filename+'.bak'))
    datanew.to_netcdf(os.path.join(datadir,filename))
     

if __name__ == '__main__':
    "Has been done now"
    bin_av('Data/ctd_files/processed2nc/Skagerak','SK_20181210_13.nc')
    
