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
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.feature as cfeature
from cartopy.io import shapereader
import cmocean as cm

plt.style.use('seaborn')
rcParams['text.usetex'] = True

def load_sst_data(datadir):
    ds=xr.open_mfdataset(os.path.join(datadir,'sst/ostia_2018121*_subset.nc'))
    sst_lon=np.asarray(ds.lon)
    sst_lat=np.asarray(ds.lat)
    sst=np.asarray((ds.analysed_sst.mean(dim='time')))

    return sst,sst_lon,sst_lat

def load_wind_data(datadir):
    dsE=xr.open_mfdataset(os.path.join(datadir,'winds/cdas1.2018121*.pgrbh.east.nc'))
    dsN=xr.open_mfdataset(os.path.join(datadir,'winds/cdas1.2018121*.pgrbh.north.nc'))

    uvel=dsE.u.mean(dim='time').squeeze()
    vvel=dsN.v.sel(lon=slice(-2,15)).mean(dim='time').squeeze()

    return uvel,vvel

def plot_sst_winds(datadir):
    [sst,lon,lat]=load_sst_data(datadir)
    [uvel,vvel]=load_wind_data(datadir)

    fig = plt.figure(figsize=(8, 8))
    fig.subplots_adjust(bottom=0.15)

    ax = plt.axes(projection=ccrs.PlateCarree(),zorder=0)
    ax.set_extent([2, 15, 50, 60],ccrs.PlateCarree())

    x,y=np.meshgrid(lon,lat)
    im = plt.pcolormesh(x, y, sst-273.15,cmap=cm.cm.thermal, zorder=1)
    cb = fig.add_axes([0.13, 0.1, 0.55, 0.02])
    cbar = plt.colorbar(im, cax=cb, orientation='horizontal', extend='both')
    cbar.set_label(label='SST [$^{\circ}$C]',size=10, labelpad=13)

    q = ax.quiver(uvel.lon, uvel.lat, uvel, vvel,color='k', scale=4e2, zorder=50, alpha=0.8)  # didn't use transform, but looks ok...
    ax.quiverkey(q, 0.89, -0.16, 10, label='Wind speed [10 m/s]')

    gl=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,alpha=0.2,linestyle='--',zorder=4)
    gl.xlabels_top=False
    gl.ylabels_right=False
    ax.add_feature(cfeature.LAND, facecolor='0.9',zorder=3)
    plt.savefig('Figures/Raw/synoptic_conditions.png')
    plt.savefig('Figures/Raw/synoptic_conditions.pdf')
    # plt.show()



if __name__ == '__main__':
    if sys.argv[1] == 'load_sst_data':
        datadir='Data'
        load_sst_data(datadir)
    elif sys.argv[1] == 'load_wind_data':
        datadir='Data'
        load_wind_data(datadir)  
    elif sys.argv[1] == 'plot_sst_winds':
        datadir='Data'
        plot_sst_winds(datadir)          
