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
            if not meta_data.empty:
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
    df=df.loc[:,:9]
    df.columns=['timeS','lon','lat','TEMP','PRES','Flr','density',
                'DEPTH','o2','PSAL']
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
        print('Gridding ' + filename)
        gridd(datadir, filename)

        
def gridd(datadir, filename, MAXDEPTH=120):
    """
    Regridding of all data along the same depth grid on *datadir*.
    The grid will be range(0,MAXDEPTH+1)

    The data should not have higher depth than MAXDEPTH

    Also change 'scan' attribute to 'DEPTH'
    Also change 'index' to 'DEPTH'
    Also drop the DEPTH of variable list
    Also change 'LATITUDE' to 'lat' and 'LONGITUDE' to 'lon'
    """
    #raise(NotImplementedError)
    data =  xr.open_dataset(os.path.join(datadir,filename))
    datanew = data.copy()
    vars_and_coords = [i for i in datanew.data_vars.keys()] + \
                      [i for i in datanew.coords.keys()]
    if 'DEPTH' in vars_and_coords:
        datanew = datanew.drop('DEPTH')
        print('Removing DEPTH')
    if 'LATITUDE' in vars_and_coords:
        datanew = datanew.rename({'LATITUDE':'lat'})
    if 'LONGITUDE' in vars_and_coords:
        datanew = datanew.rename({'LONGITUDE':'lon'})
    #print(datanew)
    for i in datanew.coords.keys():
        # changing coordinates to variables
        line = datanew[i]
        datanew = datanew.drop(i)
        datanew[i] = ('DEPTH', line)
    for old_str in ['scan', 'index']:
        if old_str in datanew.dims:
            datanew=datanew.rename({old_str:'DEPTH'})
            print('renaming scan/index')
    data_maxd = data.DEPTH.max()
    data_mind = data.DEPTH.min()
    if MAXDEPTH < data_maxd:
        print('Skiping data, depth is higher than MAXDEPTH')
        return 0
    datanew=datanew.rename({'DEPTH':'DEPTH_old'})
    #print([i for i in datanew.keys()])
    #print(datanew)
    vars_and_coords = [i for i in datanew.data_vars.keys()] + \
                      [i for i in datanew.coords.keys()]
    keys_tot = [i for i in vars_and_coords if i!='DEPTH_old']
    for key in keys_tot:
        line = data[key]
        line_new = np.zeros(MAXDEPTH) * np.NaN
        line_new[int(data_mind):line.shape[0]+int(data_mind)] = line
        datanew[key] = ('DEPTH', line_new)
        #print(key)
    #raise NotImplementedError
    vars_and_coords = [i for i in datanew.data_vars.keys()] + \
                      [i for i in datanew.coords.keys()]
    #print(vars_and_coords)
    #print(datanew)
    if 'DEPTH_old' in vars_and_coords:
        datanew = datanew.drop('DEPTH_old')
        print('Removing DEPTH_old')
    #print(datanew.DEPTH)
    datanew.to_netcdf(os.path.join(datadir, filename[:-3]+'_grid.nc'))
    
if __name__ == '__main__':
    #Has been done
    #add_coordinates('Data/ctd_files/processed2nc','Data/ctd_files/meta/TB20181210_meta.csv')
    #save_ascii2nc('Data/Skagerak/SK20181210/SK20181210_CTD/SK20181210_Processed_data')
    gridd_all('Data/ctd_files/processed2nc/Trygve')
