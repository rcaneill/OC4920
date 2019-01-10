import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import os
import gsw
from math import atan2, degrees
from matplotlib import rcParams, rc

plt.style.use('seaborn')
rcParams['text.usetex'] = True

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)




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
    colorcycle = [i['color'] for i in rcParams['axes.prop_cycle']]
    c = np.array([colorcycle[0] for i in ta])
    c[ta > 0] = colorcycle[1]
    plt.scatter(p, ta, c=c, s=2)
    return ([i for i in p], [i for i in ta])

if __name__ == '__main__':
    datadir = 'Data/ctd_files/gridded_calibrated_updated'
    fig = plt.figure(figsize=(6,4))
    p = []
    ta = []
    for filename in os.listdir(datadir):
        (a,b) = turner(datadir, filename)
        p += a
        ta += b
    p = np.array(p)
    ta = np.array(ta)
    d_layers = [11.8, 60.2, 75.3]
    print('********************')
    print('Upper cline')
    print('ta mean, std', np.nanmean(ta[p<11.8]), np.nanstd(ta[p<11.8]))
    print('Warm water')
    print('ta mean, std', np.nanmean(ta[(p>11.8) & (p<60.2)]), \
          np.nanstd(ta[(p>11.8) & (p<60.2)]))
    print('Lower cline')
    print('ta mean, std', np.nanmean(ta[(p>60.2) & (p<75.3)]), \
          np.nanstd(ta[(p>60.2) & (p<75.3)]))
    print('Bottom water')
    print('ta mean, std', np.nanmean(ta[p>73.5]), np.nanstd(ta[p>73.5]))
    print('********************')
    colorcycle = [i['color'] for i in rcParams['axes.prop_cycle']]
    r = colorcycle[2]
    plt.plot([0,0],[-180,180],c=r)
    plt.plot([11.8,11.8],[-180,180],c=r)
    plt.plot([60.2,60.2],[-180,180],c=r)
    plt.plot([75.3,75.3],[-180,180],c=r)
    plt.yticks([-90,-45,0,45,90])
    plt.ylim(-91,91)
    plt.ylabel('Turner angle ($^\circ$)')
    plt.xlabel('Depth (m)')
    fig.tight_layout()
    plt.savefig('Figures/ForReport/turner.pdf')
    plt.show()
    
