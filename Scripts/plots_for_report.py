import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import os, sys
from matplotlib import rcParams
from plot_profile import stations
import cartopy.crs as ccrs
import cartopy.feature as cfeature

plt.style.use('seaborn')
rcParams['text.usetex'] = True


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
    
if __name__ == '__main__':
    argument = sys.argv[1]
    if argument == 'prof_fit':
        # plot the fitted profile
        prof_fit()

    if argument == 'stations':
        stations('Data/', report=True)

    if argument == 'sk':
        # plot skagerrak
        plot_sk()

    # create hard link of the figures to report figures
    all_fig = ['ForReport/explainProfFit.pdf', \
               'ForReport/stations.pdf', \
               'ForReport/ts.pdf', \
               'ForReport/turner.pdf']
    for i in all_fig:
        if os.path.basename(i) in os.listdir('Report/Romain/Figures'):
            os.remove(os.path.join('Report/Romain/Figures',os.path.basename(i)))
        os.link(os.path.join('Figures',i), \
                os.path.join('Report/Romain/Figures',os.path.basename(i)))
        
