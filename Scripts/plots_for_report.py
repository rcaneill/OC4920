import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os, sys
from matplotlib import rcParams, rc
from plot_profile import stations
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.style.use('seaborn')
rcParams['text.usetex'] = True
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


def prof_fit():
    datadir = 'Data/ctd_files/fitted'
    dataG = xr.open_dataset(os.path.join(datadir, 'TB_20181210_07_grid.nc')) # Good example
    dataB = xr.open_dataset(os.path.join(datadir, 'TB_20181210_19_grid.nc')) # Bad example
    fig,ax = plt.subplots(1, 2, figsize=(5,4), sharey=True)
    ax[0].plot(dataG.TEMP_f, dataG.DEPTH, label='fit')
    ax[0].plot(dataG.ptemp, dataG.DEPTH, label='measure')
    ax[1].plot(dataB.TEMP_f, dataB.DEPTH, label='fit')
    ax[1].plot(dataB.ptemp, dataB.DEPTH, label='measure')
    ax[0].legend()
    ax[1].legend()
    xlim = ax[0].get_xlim()
    ylim = ax[0].get_ylim()
    for (d,v) in zip(dataG.layer_d, dataG.layer_v):
        ax[0].scatter(v, d, c='r', s=80, marker=u'x')
    for (d,v) in zip(dataB.layer_d, dataB.layer_v):
        ax[1].scatter(v, d, c='r', s=80, marker=u'x')
    ax[0].invert_yaxis()
    ax[0].set_xlabel(u'Temperature ($^{\circ}$C)')
    ax[1].set_xlabel(u'Temperature ($^{\circ}$C)')
    ax[0].set_ylabel('Depth (m)')
    ax[0].set_title('Good fit of a profile')
    ax[1].set_title('Only upper layer fit')
    plt.tight_layout()
    plt.savefig('Figures/ForReport/explainProfFit.pdf')
    plt.close()
    
    
def normalize(d):
    d -= np.min(d)
    d /= np.ptp(d)
    return d
    
def layer_stat(datadir):
    """
    Computes stats about the layers
    """
    d = []
    v = []
    flag = []
    lon = []
    lat = []
    for filename in os.listdir(datadir):
        data = xr.open_dataset(os.path.join(datadir,filename))
        d.append(data.layer_d)
        v.append(data.layer_v)
        flag.append(data.layer_flag)
        lon.append(data.lon.values[~np.isnan(data.lon.values)][0])
        lat.append(data.lat.values[~np.isnan(data.lat.values)][0])
    d = np.array(d)
    v = np.array(v)
    flag = np.array(flag)
    lat = np.array(lat)
    lon = np.array(lon)
    d[flag==0] = np.nan
    v[flag==0] = np.nan
    dist = normalize(np.sqrt((lat-58.24)**2+(lon-11.4)**2))
    labs = ['Top of warm water', 'Bottom of warm water', \
            'Bottom of lower thermocline']
    fig,ax = plt.subplots(1,2, figsize=(6,4))
    for (i,var) in enumerate(v.T[1:]):
        ax[0].scatter(dist, var, label=labs[i])
    for (i,var) in enumerate(d.T[1:]):
        ax[1].scatter(dist, var, label=labs[i])
    ax[0].set_xlabel(u'Distance (arbitrary unit)')
    ax[0].set_ylabel(u'Interface temperature ($^{\circ}$C)')
    ax[1].set_xlabel(u'Distance (arbitrary unit)')
    ax[1].set_ylabel(u'Interface depth (m)')
    ax[1].legend()
    fig.tight_layout()
    plt.savefig('Figures/ForReport/layers.pdf')
    plt.show()
    print('Mean depth',np.nanmean(d,axis=0))
    print('Std depth',np.nanstd(d,axis=0))
    print('Mean temp',np.nanmean(v,axis=0))
    print('Std temp',np.nanstd(v,axis=0))
    print('Number of values for each depth', np.sum(flag==1,axis=0))
    print('Mean height of the warm layer',np.nanmean(d.T[2]-d.T[1]))
    print('Std height of the warm layer',np.nanstd(d.T[2]-d.T[1]))
    plt.scatter(dist, d.T[2]-d.T[1])
    plt.show()
    
if __name__ == '__main__':
    argument = sys.argv[-1]
    if argument in ['prof_fit','all']:
        # plot the fitted profile
        prof_fit()

    if argument in ['stations', 'all']:
        stations('Data/', report=True)

    if argument in ['sk', 'all']:
        # plot skagerrak
        #plot_sk()
        pass

    if argument in ['layers', 'all']:
        layer_stat(datadir='Data/ctd_files/fitted')
        
    # create hard link of the figures to report figures
    all_fig = ['ForReport/explainProfFit.pdf', \
               'ForReport/stations.pdf', \
               'ForReport/ts.pdf', \
               'ForReport/turner.pdf', \
               'ForReport/layers.pdf', \
               'ForReport/transect.pdf', \
               'Layers/layer_bottom_of_thermocline.pdf', \
               'Layers/layer_bottom_of_warm_water.pdf']
    for i in all_fig:
        if os.path.basename(i) in os.listdir('Report/Romain/Figures'):
            os.remove(os.path.join('Report/Romain/Figures',os.path.basename(i)))
        os.link(os.path.join('Figures',i), \
                os.path.join('Report/Romain/Figures',os.path.basename(i)))
        
