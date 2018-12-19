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
from cartopy.io import shapereader
import cmocean as cm


plt.style.use('seaborn')
rcParams['text.usetex'] = True


def prof(filename):
    """
    plot profile of the data.

    Arguments
    ---------
        filename: name of the file, string
            Assumes data are in 'Data/ctd_files/gridded
    """
    name = filename[:-3]
    datadir = "Data/ctd_files/gridded"
    data = xr.open_dataset(os.path.join(datadir,filename))
    fig,ax = plt.subplots(1,2,sharey=True)
    temp = data['TEMP']
    psal = data['PSAL']
    depth = data['DEPTH']
    #sigma_t = data['sigma_t']

    axe0 = ax[0]
    axe1 = ax[1]
    axe2 = axe0.twiny()
    axe0.plot(psal, depth, 'r', label='Practical salinity')
    axe2.plot(temp, depth, 'b', label='Temperature')
    #axe1.plot(sigma_t, depth, 'g', label=u'$\sigma_t$')
    fig.suptitle(name.replace('_',' '))
    for axe in (axe0, axe1, axe2):
        #axe.invert_yaxis()
        axe.legend()
        axe.set_ylabel('Depth (m)')
        axe.set_ylim(120,0)
    axe0.set_xlabel('Salinity (psu)')
    axe2.set_xlabel(u'Temperature ($^{\circ}$C)')
    axe1.set_xlabel(u'Density anomaly $kg/m^3$')
    fig.savefig('Figures/Profiles/profile_' + name + '.pdf')
    fig.savefig('Figures/Profiles/profile_' + name + '.png')
    plt.close()

def all_prof(datadir):
    """
    Plot profile of all data on the sqme plot
    """
    fig,ax = plt.subplots(1,2,sharey=True)
    for filename in os.listdir(datadir):
        data = xr.open_dataset(os.path.join(datadir,filename))
        ax[0].plot(data.PSAL, data.DEPTH, 'k', linewidth=0.5)
        ax[1].plot(data.TEMP, data.DEPTH, 'k', linewidth=0.5)
    ax[0].set_title('Salinity')
    ax[1].set_title('Temperature')
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_xlabel('Salinity (psu)')
    ax[1].set_xlabel(u'Temperature ($^{\circ}$C)')
    for axe in ax:
        #axe.invert_yaxis()
        axe.set_ylim(120,0)
    fig.savefig('Figures/Raw/profileTot.pdf')
    fig.savefig('Figures/Raw/profileTot.png')
    plt.show()
    
def ts(datadir):
    """
    Plot T S for all files
    """
    temp = []
    sal = []
    depth = []
    for filename in os.listdir(datadir):
        data = xr.open_dataset(os.path.join(datadir,filename))
        for (t,s,d) in zip(data.TEMP.values, data.PSAL.values, data.DEPTH.values):
            temp.append(t)
            sal.append(s)
            depth.append(d)
    scatt = plt.scatter(sal,temp,c=depth,s=3,cmap='viridis_r')
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
    etopo=xr.open_dataset(datadir+'etopo1_bedrock.nc')
    
    test = shapereader.Reader(datadir+'coastline.shp')

    #load station lat and lons
    dfT1=pd.read_csv(datadir+'ctd_files/meta/TB20181210_meta.csv',header=[0])
    dfT2=pd.read_csv(datadir+'ctd_files/meta/TB20181211_meta.csv',header=[0])
    dfS=xr.open_dataset(datadir+'ctd_files/meta/meta_SK.nc')

    #define axes
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=11))
    ax.set_extent([11.2, 11.8, 58.1, 58.5],ccrs.PlateCarree())
    geometries = [i for i in test.geometries()]
    for geometry in geometries:
        ax.add_geometries([geometry], ccrs.PlateCarree(), facecolor='lightgray',
                          edgecolor='black',zorder=2)
    levels=np.linspace(-140,0,8)
#    etopo.Band1.plot.contourf(ax=ax,levels=levels,transform=ccrs.PlateCarree(),zorder=0,cbar_kwargs={'label':'Bathymetry (m)'})
    etopo.Band1.plot.pcolormesh(ax=ax,vmin=-60,vmax=20,transform=ccrs.PlateCarree(),cmap=cm.cm.solar,zorder=0,cbar_kwargs={'label':'Bathymetry (m)','fraction':0.035})
    
    gl=ax.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,alpha=0.2,zorder=1)
    ax.scatter(dfT1.lon,dfT1.lat,s=5,c='w',transform=ccrs.PlateCarree())
    ax.scatter(dfT2.lon,dfT2.lat,s=5,c='w',transform=ccrs.PlateCarree())
    ax.scatter(dfS.lon,dfS.lat,s=5,c='w',transform=ccrs.PlateCarree())
    gl.xlabels_top=False
    gl.ylabels_right=False
    plt.savefig('Figures/Raw/stations.png')
    plt.savefig('Figures/Raw/stations.pdf')
    # add labels with the cast names
    #print(dfT1)
    for label, x, y in zip(dfT1.filename, dfT1.lon, dfT1.lat):
        #print(label)
        ax.text(x,y,str(label).replace('_',' '), transform=ccrs.PlateCarree())
    for label, x, y in zip(dfT2.filename, dfT2.lon, dfT2.lat):
        #print(label)
        ax.text(x,y,str(label).replace('_',' '), transform=ccrs.PlateCarree())
    for label, x, y in zip(dfS.station.values, dfS.lon, dfS.lat):
        #print(label)
        ax.text(x,y,str(label).replace('_',' '), transform=ccrs.PlateCarree())
    plt.show()
    
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

def ship_calib(datadir,filename_arr):
    """
    Create subplots of ship calibration casts
    filename_arr: Array of calibration files
    """
    ### Load calibration files
    SK1=xr.open_dataset(os.path.join(datadir,filename_arr[0]))
    SK2=xr.open_dataset(os.path.join(datadir,filename_arr[1]))
    TB1=xr.open_dataset(os.path.join(datadir,filename_arr[2]))
    TB2=xr.open_dataset(os.path.join(datadir,filename_arr[3]))

    # temp
    fig,axes=plt.subplots(1,2,figsize=([10,6]))
    axes[0].plot(SK2.TEMP,SK2.DEPTH,label='SK day2 temp')
    axes[0].plot(TB2.TEMP,TB2.DEPTH,label='TB day2 temp')
    axes[1].plot(SK1.TEMP,SK1.DEPTH,label='SK day1 temp')
    axes[1].plot(TB1.TEMP,TB1.DEPTH,label='TB day1 temp')
    axes[1].invert_yaxis()
    axes[0].invert_yaxis()
    axes[1].legend()
    axes[0].legend()
    plt.savefig('Figures/Calib/profile_temp.png')
    plt.savefig('Figures/Calib/profile_temp.pdf')
    plt.show()
    plt.close()

    # sal
    
    fig,axes=plt.subplots(1,2,figsize=([10,6]))
    # #plot density
    axes[0].plot(SK2.PSAL,SK2.DEPTH,label='SK day2 sal')
    axes[0].plot(TB2.PSAL,TB2.DEPTH,label='TB day2 sal')
    #axes[0].scatter(TB1.PSAL,SK1.PSAL)
    axes[1].plot(SK1.PSAL,SK1.DEPTH,label='SK day1 sal')
    axes[1].plot(TB1.PSAL,TB1.DEPTH,label='TB day1 sal')
    axes[1].invert_yaxis()
    axes[0].invert_yaxis()
    axes[1].legend()
    axes[0].legend()
    plt.savefig('Figures/Calib/profile_sal.png')
    plt.savefig('Figures/Calib/profile_sal.pdf')
    plt.show()
    plt.close()

def calibrated_profs(datadir):
    """
    plot of calibrated profiles
    """
    TB1=xr.open_dataset(os.path.join(datadir,'TB_2018121cal_down_grid.nc'))
    SK1=xr.open_dataset(os.path.join(datadir,'SK_20181210_Calibration_grid.nc'))
    SK2=xr.open_dataset(os.path.join(datadir,'SK_20181211_01_grid.nc'))
    TB2=xr.open_dataset(os.path.join(datadir,'TB_20181211_cal_down_grid.nc'))

    fig,ax=plt.subplots(1,2)
    ax[0].plot(TB1.TEMP)
    ax[0].plot(SK1.TEMP)
    ax[1].plot(TB1.t_corrected)
    ax[1].plot(SK1.TEMP)

    plt.savefig('Figures/Calib/profile_temp_calib.png')
    plt.savefig('Figures/Calib/profile_temp_calib.pdf')
    plt.show()
    plt.close()


    fig,ax=plt.subplots(1,2)
    ax[0].plot(TB1.PSAL)
    ax[0].plot(SK1.PSAL)
    ax[1].plot(TB1.s_corrected)
    ax[1].plot(SK1.PSAL)

    plt.savefig('Figures/Calib/profile_salt_calib.png')
    plt.savefig('Figures/Calib/profile_salt_calib.pdf')
    plt.show()
    plt.close()


    fig,ax=plt.subplots(1,2)
    ax[0].plot(TB2.TEMP)
    ax[0].plot(SK2.TEMP)
    ax[1].plot(TB2.t_corrected)
    ax[1].plot(SK2.TEMP)

    plt.savefig('Figures/Calib/profile_temp_calib_day2.png')
    plt.savefig('Figures/Calib/profile_temp_calib_day2.pdf')
    plt.show()
    plt.close()

    fig,ax=plt.subplots(1,2)
    ax[0].plot(TB2.PSAL)
    ax[0].plot(SK2.PSAL)
    ax[1].plot(TB2.s_corrected)
    ax[1].plot(SK2.PSAL)

    plt.savefig('Figures/Calib/profile_salt_calib_day2.png')
    plt.savefig('Figures/Calib/profile_salt_calib_day2.pdf')
    plt.show()
    plt.close()



    
if __name__ == '__main__':
    if sys.argv[1] == 'profile':
        all_names = os.listdir('Data/ctd_files/gridded')
        for filename in all_names:
            #prof('TB_20181210_03_down.nc')
            prof(filename)

    elif sys.argv[1] == 'all_prof':
        all_prof('Data/ctd_files/gridded')
            
    elif sys.argv[1] == 'ts':
        ts('Data/ctd_files/gridded')

    elif sys.argv[1] == 'stations':
        stations('Data/')

    elif sys.argv[1] == 'section':
        data_all = []
        for filename in os.listdir('Data/ctd_files'):
            data_all.append(load_nc(filename))
        print(len(data_all))
        meta=pd.read_csv('Data/Trygve/TB20181210_meta_edit.csv',header=[0])
        section(data_all[:4], meta)    
        
    elif sys.argv[1] == 'ship_calib':
        filename_arr=['SK_20181210_Calibration_grid.nc', \
                      'SK_20181211_01_grid.nc', \
                      'TB_2018121cal_down_grid.nc', \
                      'TB_20181211_cal_down_grid.nc']
        ship_calib('Data/ctd_files/gridded/',filename_arr)
    
    elif sys.argv[1] == 'calibrated_profs':
        datadir='Data/ctd_files/gridded_correlated'
        calibrated_profs(datadir)
