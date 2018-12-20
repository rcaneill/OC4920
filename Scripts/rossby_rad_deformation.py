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

# plt.style.use('seaborn')
rcParams['text.usetex'] = True

#Note: Typical width of Gullmar fjord is 1-3km

#########TODO: Compute rossby radius of deformation

#Lr=N*H/pi*i where i is the order. Here we compute first order


def brunt_vaisala(datadir,filenames,desc=''):

    #Load required variables
    lat = []
    lon = []
    temp = []
    sal = []
    depth = []
    for i in filenames:
        d = xr.open_dataset(os.path.join(datadir, i))
        temp.append(d.TEMP.values)
        sal.append(d.PSAL.values)
        lon.append(np.nanmean(d.lon.values))
        lat.append(np.nanmean(d.lat.values))
       
    depth = -np.abs(d.DEPTH.values)
    coords=lon
    #compute alpha and beta
    [rho,alpha,beta]=gsw.rho_alpha_beta(sal,temp,depth)

    g=9.8
    N = np.sqrt(g*((alpha[:,:-1]*(np.diff(temp)/np.diff(depth)))
    -(beta[:,:-1]*(np.diff(sal)/np.diff(depth))))) ##computed this way so that the influence of temperature 
    #vs salinity can be separated
    #plot N
    xx,yy=np.meshgrid(lon,depth[:-1])
    plt.pcolormesh(xx,yy,N.T)
    cb=plt.colorbar(extend='both')
    cb.set_label('N')
    
    #plt.show()
    
    plt.savefig('Figures/Transect/transect_N.pdf')
    plt.savefig('Figures/Transect/transect_N.png')
    plt.close()
    return N,depth,lon,lat

def rossby_radii(datadir,filenames,desc=''):
    [N,depth,lon, lat]=brunt_vaisala(datadir,filenames,desc='')
    #Length scale of rossby radius of deformation 
    def coriolis_acc(lat,omega=0.729*10**-4):
        return 2*omega*np.sin(lat/360.*2*np.pi)
    f=coriolis_acc(np.asarray(lat))
    Lr=(N*depth[:-1]*-1)/(np.nanmean(f))*0.001

    plt.figure()
    xx,yy=np.meshgrid(lon,depth[:-1])
    plt.pcolormesh(xx,yy,Lr.T)
    cb=plt.colorbar(extend='both')
    plt.show()
    # cb.set_label('N')
    print('max:',np.nanmax(Lr))
    print('min:',np.nanmin(Lr))
    print('mean:',np.nanmean(Lr))

    # np.savetxt('Data/rossby_radius_deformation.txt',np.c_[, b_temp, a_sal, b_sal], \
    #            header="a_temp b_temp a_sal b_sal",comments='')

    return Lr


def load_filenames_section(N, sec_datadir='Data/sections'):
    """
    load and return the filenames + decription of the section

    data must be on a txt file like: secN.txt
    N must be replaced with the section number (e.g. 'sec1.txt', 'sec5.txt')
    """
    return np.loadtxt(os.path.join(sec_datadir, 'sec'+str(N)+'.txt'), \
                      dtype=str, delimiter='\n')
        

    
if __name__ == '__main__':
    if sys.argv[1] == 'brunt_vaisala':
        datadir = 'Data/ctd_files/gridded_calibrated'
        section_meta = load_filenames_section(0, 'Data/sections')
        brunt_vaisala(datadir, section_meta[1:], desc=section_meta[0])
    if sys.argv[1] == 'rossby_radii':
        datadir = 'Data/ctd_files/gridded_calibrated'
        section_meta = load_filenames_section(0, 'Data/sections')
        rossby_radii(datadir, section_meta[1:], desc=section_meta[0])

