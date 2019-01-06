import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
import xarray as xr
import gsw
from matplotlib import rcParams
from load_data import load_nc
import os, sys

def compute_rp(filename):
    """
    Compute Rp - the ratio of the change in salinity 
    relative to temperature. Allows for identification
    of profiles which are susceptible to double diffusion or salt
    fingering
    """
    plt.style.use('seaborn')
    name = filename[:-3]
    datadir = "Data/ctd_files/gridded_calibrated_updated"
    data = xr.open_dataset(os.path.join(datadir,filename))
    # fig,ax = plt.subplots(1,2,sharey=True)
    temp = data['ptemp']
    psal = data['ab_sal']
    depth = data['DEPTH']
    sigma_t = data['Gsw_sigma0A0']
    
    [rho,alpha,beta]=gsw.rho_alpha_beta(psal,temp,depth)
    rp=(beta[:-1]*np.diff(psal))/(alpha[:-1]*np.diff(temp))
    
    fig,ax=plt.subplots(2,2)
    fig.suptitle(name.replace('_',' '))
    ax[1,1].plot(rp,depth[:-1]*-1)
    ax[1,0].plot(rho,depth*-1)
    ax[0,1].plot(psal,depth*-1)
    ax[0,0].plot(temp,depth*-1)
    plt.show()
    
if __name__ == '__main__':
    if sys.argv[1] == 'compute_rp':
        all_names = os.listdir('Data/ctd_files/gridded_calibrated_updated')
        for filename in all_names:
            # compute_rp('SK_20181211_13_grid.nc')
            compute_rp(filename)

