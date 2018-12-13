import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
import xarray as xr
import gsw
from matplotlib import rcParams
from load_data import load_nc
import os, sys
import pandas as pd
import cartopy.crs as ccrs

plt.style.use('seaborn')
rcParams['text.usetex'] = True

def prof(filename):
    """
    plot profile of the data.

    Arguments
    ---------
        filename: name of the file, string
            Assumes data are in 'Data/ctd_files/
    """
    name = filename[:-3]
    datadir = "Data/ctd_files/"
    data = load_nc(filename)
    fig,ax = plt.subplots(1,2,sharey=True)
    temp = data['TEMP']
    psal = data['PSAL']
    depth = data['DEPTH']
    sigma_t = data['sigma_t']

    axe0 = ax[0]
    axe1 = ax[1]
    axe2 = axe0.twiny()
    axe0.plot(psal, depth, 'r', label='Practical salinity')
    axe2.plot(temp, depth, 'b', label='Temperature')
    axe1.plot(sigma_t, depth, 'g', label=u'$\sigma_t$')
    fig.suptitle(name.replace('_',' '))
    for axe in (axe0, axe1, axe2):
        axe.invert_yaxis()
        axe.legend()
        axe.set_ylabel('Depth (m)')
    axe0.set_xlabel('Salinity (psu)')
    axe2.set_xlabel(u'Temperature ($^{\circ}$C)')
    axe1.set_xlabel(u'Density anomaly $kg/m^3$')
    fig.savefig('Figures/Profiles/profile_' + name + '.pdf')
    fig.savefig('Figures/Profiles/profile_' + name + '.png')
    plt.close()

def ts(datadir):
    """
    Plot T S for all files
    """
    all=xr.open_mfdataset(datadir+'*.nc')
    scatt = plt.scatter(all.PSAL,all.TEMP,c=all.DEPTH,s=3,cmap='viridis_r')
    cb = plt.colorbar(scatt)
    cb.set_label('Depth (m)')
    plt.xlabel('Practical salinity (psu)')
    plt.ylabel(u'Temperature ($^{\circ}$C)')
    plt.savefig('Figures/Raw/ts.png')
    plt.savefig('Figures/Raw/ts.pdf')

def stations(datadir):
    """
    plot station locations
    """
    #load bathymetry data
    etopo=xr.open_dataset(datadir+'etopo1.nc')
    df=pd.read_csv(datadir+'Trygve/TB20181210_meta_edit.csv',header=[0])
    ax=plt.axes(projection=ccrs.PlateCarree())
    ax.scatter(df.lon,df.lat,transform=ccrs.PlateCarree())
    ax.etopo.Band1.contour()
    plt.show()
    plt.savefig('Figures/Raw/stations.png')
    plt.savefig('Figures/Raw/stations.pdf')
    
def section(data_tot,meta):
    """
    plot temp, salt, density contours overlain, sections across fjord

    Arguments
    ---------
        data_all: array containing the data
            [data1, data2, data3, ...]
    """
    
    
    sal = []
    temp = []
    density = []
    depth = []
    depth_m = []
    
    for data in data_tot:
        sal.append(data.PSAL)
        temp.append(data.TEMP)
        density.append(data.sigma_t)
        depth.append(data.DEPTH)
        depth_m.append(np.max(data.DEPTH))
    print(len(data_tot))
    # make all data same depth
    depth_max = int(np.max(depth_m))
    sal_gridded = np.ones((len(data_tot), depth_max)) * np.NaN
    for i in range(len(data_tot)):
        print(sal_gridded[i,0:sal[i].shape[0]].shape, sal[i].shape)
        sal_gridded[i,0:sal[i].shape[0]] = sal[i]
    
    fig,ax=plt.subplots(2,1,sharex=True)
    ax[0].imshow(sal_gridded.T, cmap='viridis', aspect='auto')
    ax[1].contourf(sal_gridded.T, cmap='viridis')
    ax[1].invert_yaxis()
    plt.show()

    


    
if __name__ == '__main__':
    if sys.argv[1] == 'profile':
        all_names = os.listdir('Data/ctd_files')
        for filename in all_names:
            #prof('TB_20181210_03_down.nc')
            prof(filename)
    elif sys.argv[1] == 'ts':
        ts('Data/ctd_files/')
    elif sys.argv[1] == 'stations':
        stations(datadir='Data/')
    elif sys.argv[1] == 'section':
        data_all = []
        for filename in os.listdir('Data/ctd_files'):
            data_all.append(load_nc(filename))
        print(len(data_all))
        meta=pd.read_csv('Data/Trygve/TB20181210_meta_edit.csv',header=[0])
        section(data_all[:4], meta)        

