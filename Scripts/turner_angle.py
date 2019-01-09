import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import gsw
from math import atan2, degrees

def turner_angle(gT, gS, alpha, beta):
    """
    Return the turner angle, in an array shape if the
    gradients are arrays

    turner angle = atan2(alpha*gT + beta*gS, alpha*gT - beta*gS)


    see johnson et al. 2012
    Relative contributions of temperature and salinity to seasonal
    mixed layer density changes and horizontal density gradients
    
    gT, gS are the temperature and salinity gradients
    """
    return np.array([degrees(atan2(a1,a2)) for (a1,a2) in \
                     zip(alpha*gT + beta*gS, alpha*gT - beta*gS)])



def turner(datadir, filename):
    """
    compute turner angle for data contained in filename
    """
    data = xr.open_dataset(os.path.join(datadir, filename))
    T = data.ptemp.values
    S = data.ab_sal.values
    lon = data.lon.values[np.isreal(data.lon.values)][0]
    lat = data.lat.values[np.isreal(data.lat.values)][0]
    p = data.DEPTH.values
    alpha = gsw.alpha(S, T, p)
    beta = gsw.beta(S, T, p)
    ta = turner_angle(np.gradient(T,-p), np.gradient(S,-p), alpha, beta)
    plt.scatter(p, ta, c='gray', s=2)

if __name__ == '__main__':
    datadir = 'Data/ctd_files/gridded_calibrated_updated'
    for filename in os.listdir(datadir):
        turner(datadir, filename)
    plt.plot([0,120],[45,45],'k')
    plt.plot([0,120],[0,0],'k')
    plt.plot([0,120],[-45,-45],'k')
    plt.plot([0,120],[90,90],'k')
    plt.plot([0,120],[-90,-90],'k')
    
    plt.plot([0,0],[-180,180],'r')
    plt.plot([20,20],[-180,180],'r')
    plt.plot([60,60],[-180,180],'r')
    plt.plot([80,80],[-180,180],'r')
    plt.savefig('Figures/ForReport/turner.pdf')
    plt.show()
    
