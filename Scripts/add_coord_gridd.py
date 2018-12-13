import numpy as np
import xarray as xr
import pandas as pd
import os

def add_coordinates(datadir, metafilename):
    """
    Add coordinates and save into *.nc, moves old files to *.nc.bak
    """
    meta=pd.read_csv(metafilename,header=[0])
    ls = os.listdir(datadir)
    ls.sort()
    for filename in ls:
        if filename[-2:] == 'nc':
            data =  xr.open_dataset(os.path.join(datadir,filename))
            meta_data = meta[meta['filename'] == filename]
            data['lat'] = np.ones(data.DEPTH.shape) * meta_data.lat.values[0]
            data['lon'] = np.ones(data.DEPTH.shape) * meta_data.lon.values[0]
            data.swap_dims({'scan':'DEPTH'},inplace=True)            
            print('Moving {0} to {0}.bak'.format(filename))
            os.rename(os.path.join(datadir,filename), \
                      os.path.join(datadir,filename+'.bak'))
            print('Saving file {}'.format(filename))
            data.to_netcdf(os.path.join(datadir,filename))

def sort_coords(datadir):
    """
    Resort coordinates and save into *.nc, moves old files to *.nc.bak
    """
    ls = os.listdir(datadir)
    ls.sort()
    for filename in ls:
        if filename[-2:] == 'nc':
            data =  xr.open_dataset(os.path.join(datadir,filename))
            data.rename({'PRES':'DEPTH'},inplace=True)
            data.swap_dims({'scan':'DEPTH'},inplace=True)            
            print('Moving {0} to {0}.bak'.format(filename))
            os.rename(os.path.join(datadir,filename), \
                      os.path.join(datadir,filename+'.bak'))
            print('Saving file {}'.format(filename))
            data.to_netcdf(os.path.join(datadir,filename))

def read_ascii(datadir,filename):
    df=pd.read_csv(os.path.join(datadir,filename),sep='\s+',header=None)
    df.columns=['timeS','lon','lat','TEMP','PRES','Flr','density',
               'depth','o2','psal']
    return df
    
def save_ascii2nc(datadir):
    ls = os.listdir(datadir)
    ls.sort()
    for filename in ls:
        if filename[-3:] == 'asc':
            data=read_ascii(datadir,filename)
            ds=data.to_xarray()
            filename=filename[:-3]+'nc'
            print('Saving file {}'.format(filename))
            # raise(NotImplementedError)
            ds.to_netcdf(os.path.join(datadir,filename))




def gridd(datadir, MAXDEPTH=120):
    """
    TODO (Romain)
    """
    raise(NotImplementedError)
            
if __name__ == '__main__':
    #Has been done
    #add_coordinates('Data/ctd_files/', 'Data/Trygve/TB20181210_meta_edit.csv')
    save_ascii2nc('/home/wizard/Documents/observing_the_ocean/analysis/Data/Skagerak/SK20181210/SK20181210_CTD/SK20181210_Processed_data/')
