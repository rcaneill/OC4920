import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
import gsw
from matplotlib import rcParams
import os, sys
import pandas as pd
from matplotlib.image import NonUniformImage
from scipy.interpolate import griddata
from cartopy.io import shapereader
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

plt.style.use('seaborn')
rcParams['text.usetex'] = True

def sec_show(coord, depth, data, axe, cmap='viridis', coord_type='lon'):
    """
    plot data on the axe, mix between imshow and contour
    return the colorbar, so user can add a label
    """
    depth = - np.abs(depth) # to be sure that depth is negative
    im = NonUniformImage(axe, interpolation='nearest', \
                         extent=(np.nanmin(coord), np.nanmax(coord), \
                                 np.nanmin(depth), np.nanmax(depth)), \
                         cmap=cmap)
    # we need to reverse order of array so everything goes increasing
    im.set_data(coord, depth[::-1], data.T[::-1,:])
    axe.images.append(im)
    axe.set_xlim(np.nanmin(coord), np.nanmax(coord))
    axe.set_ylim(np.nanmin(depth), np.nanmax(depth))
    axe.set_ylabel('Depth (m)')
    axe.set_xlabel(coord_type+u' ($^{\circ}$ N/E)')

    # adding colorbar
    cb = plt.colorbar(im,ax=axe)
    # printing scale with positive depth
    #axe.set_yticklabels(np.abs([int(i) for i in axe.get_yticks()]))
    return cb

def compute_surface_layer(var,depth,ref_depth=5,threshold=0.6):
    """
    Function to compute winter water surface layer
    Select temperature or density as your variable,
    Adjust threshold and ref depth accordingly
    """
    return depth[(abs((var-var[ref_depth]))>=threshold)].min()



def plot_sec(datadir, filenames, desc='', coord_type='lon', N=None):
    """
    Plot a section.

    Coord must be increasing between the data
    
    filenames is an array containing names of all data that will be plotted
    desc : description of the section
    coord : 'lat' or 'lon' coordinate along which we do the transect
    N : number of the section
    """
    lon = []
    lat = []
    temp = []
    sal = []
    depth = []
    pdens = []
    for i in filenames:
        #print(i)
        d = xr.open_dataset(os.path.join(datadir, i))
        lat.append(np.nanmean(d.lat.values))
        lon.append(np.nanmean(d.lon.values))
        temp.append(d.TEMP.values)
        sal.append(d.PSAL.values)
        pdens.append(gsw.rho(d.PSAL.values,d.TEMP.values,d.DEPTH.values))
    depth = -np.abs(d.DEPTH.values)
    if coord_type == 'lon':
        coord = lon
    elif coord_type == 'lat':
        coord = lat
    else:
        raise ValueError("coord must be 'lon' or 'lat'")
    temp = np.array(temp)
    sal = np.array(sal)
    pdens= np.array(pdens)
    
    #sl = np.array([compute_surface_layer(x,depth) for x in temp])
    fig,ax = plt.subplots(2,2, sharey=True)
    cb_temp = sec_show(coord, depth, temp, ax[0,0], coord_type=coord_type, cmap='jet')
    cb_temp.set_label(u'Temperature ($^{\circ}C$)')
    
    ax[0,1].plot(temp.T,depth)
    ax[0,1].set_xlabel(u'Temperature ($^{\circ}C$)')
    cb_sal = sec_show(coord, depth, sal, ax[1,0], coord_type=coord_type, cmap='viridis_r')
    cb_sal.set_label('Salinity (psu)')
    ax[1,1].plot(sal.T,depth)
    ax[1,1].set_xlabel('Salinity (psu)')
    test=fig.suptitle(desc)
    test.set_fontsize(test.get_fontsize()+3)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig('Figures/Transect/transect_{}.pdf'.format(N))
    fig.savefig('Figures/Transect/transect_{}.png'.format(N))
    plt.show()
    #plt.plot(sl)
    #plt.show()

def first_non_nan(array):
    """
    return first non nan element of an array
    """
    i=0
    value=np.nan
    while np.isnan(value) and i<len(array):
        value=array[i]
        i+=1
    return value
    
def plot_surface(datadir, surf_depth=3):
    """
    plot surface properties

    Experimental, not finished...
    """
    #surf_depth ref depth for plotting
    temp=[]
    sal = []
    lon=[]
    lat=[]
    depth = []
    filenames = [i for i in os.listdir(datadir) if i not in ['TB_20181210_15_down_grid.nc',\
                                                             'TB_20181210b_down_grid.nc',\
                                                             'TB_2018121cal_down_grid.nc',\
                                                             'TB_20181211_cal_down_grid.nc',
                                                             'SK_20181210_01_grid.nc']]
    filename = os.listdir(datadir)
    for i in filenames:
        #print(i)
        d = xr.open_dataset(os.path.join(datadir,i))
        #temp.append(first_non_nan(d.TEMP.values))
        temp.append(d.ptemp.values[surf_depth])
        depth.append(d.DEPTH.values[~np.isnan(d.TEMP.values)][0])
        #print(i,d.DEPTH.values[~np.isnan(d.TEMP.values)][0])
        #sal.append(first_non_nan(d.PSAL.values))
        sal.append(d.ab_sal.values[surf_depth])
        lon.append(first_non_nan(d.lon.values))
        lat.append(first_non_nan(d.lat.values))
    # compute gridding of data
    grid_lon, grid_lat = np.mgrid[11.2:11.8:200j, 58.1:58.5:200j]
    points = [[i,j] for (i,j) in zip(lon,lat)]
    values = temp
    #values = depth
    grid_temp = griddata(points, values, (grid_lon, grid_lat), method='linear')
    # #define axes
    axe = plt.axes(projection=ccrs.PlateCarree(central_longitude=11))
    axe.set_extent([11.2, 11.8, 58.1, 58.5],ccrs.PlateCarree())
    cf = axe.contourf(grid_lon, grid_lat, \
                      grid_temp, cmap='jet', transform=ccrs.PlateCarree(), levels=20)
    cb=plt.colorbar(cf)
    cb.set_label(u'Potential Temperature ($^{\circ}C$)')
    #axe.scatter(grid_lon.flatten(), grid_lat.flatten(), \
    #            c=grid_temp.flatten(), cmap='jet', transform=ccrs.PlateCarree())
    #axe.scatter(lon,lat,c=temp,s=50,cmap='jet', transform=ccrs.PlateCarree())
    coast = shapereader.Reader('Data/topo/coastline.shp')
    geometries = [i for i in coast.geometries()]
    for geometry in geometries:
        axe.add_geometries([geometry], ccrs.PlateCarree(), facecolor='lightgray',\
                           edgecolor='black')
    gl=axe.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,alpha=0.1,zorder=1)
    gl.xlabels_top=False
    gl.ylabels_right=False
    plt.title('Temperature at {} m depth'.format(surf_depth))
    plt.savefig('Figures/Surface/{}_temp.png'.format(surf_depth))
    plt.savefig('Figures/Surface/{}_temp.pdf'.format(surf_depth))
    plt.show()

    values = sal
    grid_sal = griddata(points, values, (grid_lon, grid_lat), method='linear')
    # #define axes
    axe = plt.axes(projection=ccrs.PlateCarree(central_longitude=11))
    axe.set_extent([11.2, 11.8, 58.1, 58.5],ccrs.PlateCarree())
    cf = axe.contourf(grid_lon, grid_lat, \
                      grid_sal, cmap='viridis_r', transform=ccrs.PlateCarree(), levels=20)
    cb=plt.colorbar(cf)
    cb.set_label(u'Absolute Salinity')
    
    coast = shapereader.Reader('Data/topo/coastline.shp')
    geometries = [i for i in coast.geometries()]
    for geometry in geometries:
        axe.add_geometries([geometry], ccrs.PlateCarree(), facecolor='lightgray',\
                           edgecolor='black')
    gl=axe.gridlines(crs=ccrs.PlateCarree(),draw_labels=True,alpha=0.1,zorder=1)
    gl.xlabels_top=False
    gl.ylabels_right=False
    plt.title('Salinity at {} m depth'.format(surf_depth))
    plt.savefig('Figures/Surface/{}_salt.png'.format(surf_depth))
    plt.savefig('Figures/Surface/{}_salt.pdf'.format(surf_depth))
    plt.show()

def load_filenames_section(N, sec_datadir='Data/sections'):
    """
    load and return the filenames + decription of the section

    data must be on a txt file like: secN.txt
    N must be replaced with the section number (e.g. 'sec1.txt', 'sec5.txt')
    """
    return np.loadtxt(os.path.join(sec_datadir, 'sec'+str(N)+'.txt'), \
                      dtype=str, delimiter='\n')
    
if __name__ == '__main__':
    datadir = 'Data/ctd_files/gridded_calibrated_updated'
    N=1 #section number
    section_meta = load_filenames_section(N, 'Data/sections')
    plot_sec(datadir, section_meta[1:], desc=section_meta[0], coord_type='lat', N=N)
    #plot_surface(datadir, surf_depth=40)
