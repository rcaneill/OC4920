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

def read_ascii_reprocessed_ctd(datadir,filename):
    print('Loading ' + filename)
    df=pd.read_csv(os.path.join(datadir,filename),sep='\s+',header=0, encoding="ISO-8859-1")
    # renaming columns
    for i in range(df.columns.values.shape[0]):
        if df.columns.values[i] == 'Longitude':
            df.columns.values[i] = 'lon'
        elif df.columns.values[i] == 'Latitude':
            df.columns.values[i] = 'lat'
        elif df.columns.values[i] == 'DepSM':
            df.columns.values[i] = 'DEPTH'
        elif df.columns.values[i] == 'Sal00':
            df.columns.values[i] = 'PSAL'
        elif df.columns.values[i] == 'T090C':
            df.columns.values[i] = 'TEMP'
        df.columns.values[i] = df.columns.values[i].replace('/','per')
    df = df.copy()
    # binning every meter compared to every half meter
    half=df[(df.DEPTH % 1) == 0.5] # every half meter
    half_up = half.copy()
    half_down = half.copy()
    # We will take the mean of the half meter -> half_new
    half_up.index = half_up.index-1
    half_down.index = half_down.index+1
    half_new=(half_up+half_down)/2
    # removing NaN
    half_new = half_new.loc[~np.isnan(half_new.TEMP)]
    # every int meter
    meter=df[(df.DEPTH % 1) == 0]
    meter.loc[half_new.index] = (meter.loc[half_new.index]+half_new)/2
    return meter
    
def save_ascii2nc(datadir, read_func, new_datadir=None):
    if new_datadir==None:
        new_datadir = datadir
    ls = os.listdir(datadir)
    ls.sort()
    for filename in ls:
        if filename[-3:] == 'asc':
            filename_new = filename[:2].upper() + filename[2:-3] + 'nc'
            filename_new = filename_new.replace(' ','_')
            #be sure to have SK... and not sk...
            data=read_func(datadir,filename)
            ds=data.to_xarray()
            print('Saving file {}'.format(filename_new))
            # raise(NotImplementedError)
            ds.to_netcdf(os.path.join(new_datadir,filename_new))

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

def gridd_all(datadir, MAXDEPTH=120, new_datadir=None):
    """
    Apply gridd to all files in datadir
    """
    if new_datadir == None:
        new_datadir = datadir
    ls = os.listdir(datadir)
    ls = [i for i in ls if i[-3:] == '.nc']
    ls.sort()
    for filename in ls:
        if filename[-2:] == 'nc':
            print('Gridding ' + filename)
            gridd(datadir, filename, new_datadir=new_datadir)

        
def gridd(datadir, filename, MAXDEPTH=120, new_datadir=None):
    """
    Regridding of all data along the same depth grid on *datadir*.
    The grid will be range(0,MAXDEPTH+1)

    The data should not have higher depth than MAXDEPTH

    Also change 'scan' attribute to 'DEPTH'
    Also change 'index' to 'DEPTH'
    Also drop the DEPTH of variable list
    Also change 'LATITUDE' to 'lat' and 'LONGITUDE' to 'lon'
    """
    if new_datadir == None:
        new_datadir = datadir
    #raise(NotImplementedError)
    data =  xr.open_dataset(os.path.join(datadir,filename))
    datanew = data.copy()
    #print(datanew)
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
    data['DEPTH'] = np.round(data['DEPTH']) # being sure that depth are integer meters
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
    datanew.to_netcdf(os.path.join(new_datadir, filename[:-3]+'_grid.nc'))
    #print(datanew.DEPTH)


def fix_left_out_files(datadir,filename,lon,lat):
    data =  xr.open_dataset(os.path.join(datadir,filename))
    data['lat'] = ('scan',np.ones(data.DEPTH.shape) * lat)
    data['lon'] = ('scan',np.ones(data.DEPTH.shape) * lon)
    print('Moving {0} to {0}.bak'.format(filename))
    os.rename(os.path.join(datadir,filename), \
    os.path.join(datadir,filename+'.bak'))
    print('Saving file {}'.format(filename))
    data.to_netcdf(os.path.join(datadir,filename))

if __name__ == '__main__':
    #Has been done
    #add_coordinates('Data/ctd_files/processed2nc','Data/ctd_files/meta/TB20181210_meta.csv')
    #save_ascii2nc('Data/Skagerak/SK20181210/SK20181210_CTD/SK20181210_Processed_data', read_func=read_ascii)
    #gridd_all('Data/ctd_files/processed2nc/Trygve')

    # fix_left_out_files('Data/ctd_files/processed2nc/Trygve','TB_20181211_cal_down.nc',11.543796,58.319344)
    # fix_left_out_files('Data/ctd_files/processed2nc/Trygve','TB_2018121cal_down.nc',11.332095, 58.250534)
    # fix_left_out_files('Data/ctd_files/processed2nc/Trygve','SK_20181210_01.nc',11.543796,58.319344)
    # fix_left_out_files('Data/ctd_files/processed2nc/Trygve','TB_20181210b_down.nc',11.261770,  58.207741)

    # Reprocessed data (Marcus)
    #save_ascii2nc('Data/ctd_files/ascii', \
    #              read_func=read_ascii_reprocessed_ctd, \
    #              new_datadir='Data/ctd_files/processed2nc')
    #add_coordinates('Data/ctd_files/processed2nc','Data/ctd_files/meta/TB20181211_meta.csv')
    #add_coordinates('Data/ctd_files/processed2nc','Data/ctd_files/meta/TB20181210_meta.csv')
    #gridd_all('Data/ctd_files/processed2nc', new_datadir='Data/ctd_files/gridded')
