import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
import xarray as xr
from matplotlib import rcParams, cm
import os, sys
import pandas as pd
from matplotlib.image import NonUniformImage
from scipy.interpolate import griddata
from cartopy.io import shapereader
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

plt.style.use('seaborn')
rcParams['text.usetex'] = True

def plot_map(axe):
    """
    plot the coastline on the axe
    """
    coast = shapereader.Reader('Data/topo/coastline.shp')
    geometries = [i for i in coast.geometries()]
    for geometry in geometries:
        axe.add_geometries([geometry], ccrs.PlateCarree(), facecolor='lightgray',\
                           edgecolor='black')

def grid(func, points, values, grid_lon, grid_lat):
    """
    grid data along regular grid, using function func to
    give weight to the data points as a function of distance
    """
    return 0

def bump(dist, l):
    return -1/l**2 * (dist-l)*(dist+l)
        
def plot_var(values, lon, lat, nx=200,ny=200,coastBool=True, cb_label=None, \
             cmap='viridis', **kwargs):
    """
    plot the var with a linear interpolation
    return the axe where it is plotted

    values : 1D array with values to plot
    lon, lat : arrays with coordinates of the points
    nx, ny : number of subdivision for the interpolation/gridding
    coastBool : bool to know if we plot the coastlines
    cb_label : label for the colorbar
    cmap : colormap for the colorbar

    **kwargs : args to pass to plt.contourf
    """
    # remove point that are NaN
    NAN = np.isnan(values)
    values = np.array(values)[~NAN]
    lon = np.array(lon)[~NAN]
    lat = np.array(lat)[~NAN]
    print('Removing {} values / {} that are NaN'.format(np.sum(NAN), len(NAN)))
    # compute gridding of data
    grid_lon, grid_lat = np.mgrid[11.2:11.8:nx*1j, 58.1:58.5:ny*1j]
    points = [[i,j] for (i,j) in zip(lon,lat)]
    grid_temp = griddata(points, values, (grid_lon, grid_lat), method='linear')
    # #define axes
    axe = plt.axes(projection=ccrs.PlateCarree(central_longitude=11))
    axe.set_extent([11.2, 11.8, 58.1, 58.5],ccrs.PlateCarree())
    cf = axe.contourf(grid_lon, grid_lat, \
                      grid_temp, cmap=cmap, transform=ccrs.PlateCarree(), **kwargs)
    cb=plt.colorbar(cf)
    if cb_label != None:
        cb.set_label(cb_label)
    gl=axe.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,alpha=0.1,zorder=1)
    gl.xlabels_top=False
    gl.ylabels_right=False
    if coastBool:
        plot_map(axe)
    return axe

def layer(datadir):
    """
    Plot map of layers
    """
    depth = []
    temp = []
    lon = []
    lat = []
    flag = []
    for filename in os.listdir(datadir):
        #print(i)
        d = xr.open_dataset(os.path.join(datadir,filename))
        lon.append(np.nanmean(d.lon.values))
        lat.append(np.nanmean(d.lat.values))
        depth.append(d.layer_d.values)
        temp.append(d.layer_v.values)
        flag.append(d.layer_flag.values)
    flag = np.array(flag)
    depth = np.array(depth)
    temp = np.array(temp)
    lon = np.array(lon)
    lat = np.array(lat)
    goodFlag = flag==1
    #plt.plot(depth[:,2], 'o')
    #plt.show()
    axe = plot_var(depth[:,2][goodFlag[:,2]], lon[goodFlag[:,2]], lat[goodFlag[:,2]], nx=500, ny=500, cmap='jet', levels=20)
    axe.scatter(lon[goodFlag[:,2]], lat[goodFlag[:,2]], c=depth[:,2][goodFlag[:,2]], cmap='jet',transform=ccrs.PlateCarree())
    plt.show()
    plt.plot(depth[:,2][goodFlag[:,2]]-depth[:,1][goodFlag[:,2]],'o')
    print(np.nanmean(depth[:,2][goodFlag[:,2]]-depth[:,1][goodFlag[:,2]]))
    print(np.nanstd(depth[:,2][goodFlag[:,2]]-depth[:,1][goodFlag[:,2]]))
    plt.show()
    

if __name__ == '__main__':
    layer('Data/ctd_files/fitted')
