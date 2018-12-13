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
            data['lat'] = ('scan',np.ones(data.DEPTH.shape) * meta_data.lat.values[0])
            data['lon'] = ('scan',np.ones(data.DEPTH.shape) * meta_data.lon.values[0])
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

def extract_meta_from_nc(datadir):
    ls = os.listdir(datadir)
    ls.sort()
    lat=[]
    lon=[]
    station=[]
    for filename in ls:
        if filename[-2:] == 'nc':
            data=xr.open_dataset(os.path.join(datadir,filename))
            if 'LATITUDE' in str(data.data_vars.keys):
                data=data.rename({'LATITUDE':'lat'})
                data=data.rename({'LONGITUDE':'lon'})
            lat.append(data.lat[0])
            lon.append(data.lon[0])
            station.append(filename)
    lon=np.asarray(lon)
    lat=np.asarray(lat)
    station=np.asarray(station)
    meta={'station':station,'lon':lon,'lat':lat}
    df=pd.DataFrame(data=meta)
    ds=df.to_xarray()
    #save to nc
    ds.to_netcdf('meta_SK.nc')

def gridd_all(datadir, MAXDEPTH=120):
    """
    Apply gridd to all files in datadir
    """
    ls = os.listdir(datadir)
    ls = [i for i in ls if i[-3:] == '.nc']
    ls.sort()
    for filename in ls:
        gridd(datadir, filename)

        
def gridd(datadir, filename, MAXDEPTH=120):
    """
    Regridding of all data along the same depth grid on *datadir*.
    The grid will be range(1,MAXDEPTH+1)

    The data should not have higher depth than MAXDEPTH

    Also change 'scan' attribute to 'DEPTH'
    Also change 'index' to 'DEPTH'
    Also drop the DEPTH of variable list
    Also change 'LATITUDE' to 'lat' and 'LONGITUDE' to 'lon'
    """
    #raise(NotImplementedError)
    data =  xr.open_dataset(os.path.join(datadir,filename))
    if MAXDEPTH < data.DEPTH.max():
        raise(ValueError('Data have a depth that is higher than MAXDEPTH'))
    datanew = data.copy()
    if 'LATITUDE' in str(data.data_vars.keys):
        datanew = datanew.rename({'LATITUDE':'lat'})
    if 'LONGITUDE' in str(data.data_vars.keys):
        datanew = datanew.rename({'LONGITUDE':'lon'})
    for old_str in ['scan', 'index']:
        if old_str in str(data.data_vars.keys):
            datanew=datanew.rename({old_str:'DEPTH'})
    datanew = datanew.drop('DEPTH')
    datanew=datanew.rename({'DEPTH':'DEPTH_old'})
    keys_tot = [i for i in datanew.keys() if i!='DEPTH_old']
    for key in keys_tot:
        line = data[key]
        line_new = np.zeros(MAXDEPTH) * np.NaN
        line_new[:line.shape[0]] = line
        datanew[key] = ('DEPTH', line_new)
    datanew.to_netcdf(os.path.join(datadir, filename[:-3]+'_gridd.nc'))
    
if __name__ == '__main__':
    #Has been done
    #add_coordinates('Data/ctd_files/', 'Data/Trygve/TB20181210_meta_edit.csv')
<<<<<<< HEAD
    save_ascii2nc()
=======
    gridd_all('Data/ctd_files/')
>>>>>>> branchB
