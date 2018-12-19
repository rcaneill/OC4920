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
    temp = []
    sal = []
    depth = []
    for i in filenames:
        d = xr.open_dataset(os.path.join(datadir, i))
        temp.append(d.TEMP.values)
        sal.append(d.PSAL.values)
        
    depth = -np.abs(d.DEPTH.values)
    
    #compute alpha and beta
    [rho,alpha,beta]=gsw.rho_alpha_beta(sal,temp,depth)

    g=9.8
    N = np.sqrt(g*(alpha[:,:-1]*(np.diff(temp)/np.diff(depth))
    -beta[:,:-1]*(np.diff(sal)/np.diff(depth))))
    
    return N

def rossby_radii(datadir,filenames,desc=''):
    N=brunt_vaisala(datadir,filenames,desc='')
    Lr=(N*80)/(1*np.pi)
    print('max:',np.nanmax(Lr))
    print('min:',np.nanmin(Lr))
    print('mean:',np.nanmean(Lr))
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

