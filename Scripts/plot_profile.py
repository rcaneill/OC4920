import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import gsw
from matplotlib import rcParams
from load_data import load_nc

plt.style.use('seaborn')
rcParams['text.usetex'] = True

def prof(filename):
    """
    plot profile of the data.

    Arguments
    ---------
        filename: name of the file, string
            Assumes data are in 'Data/ctd_files/
    """
    name = filename[:-3]
    datadir = "Data/ctd_files/"
    data = load_nc(filename)
    fig,ax = plt.subplots(1,2,sharey=True)
    temp = data['TEMP']
    psal = data['PSAL']
    depth = data['DEPTH']
    sigma_t = data['sigma_t']

    axe0 = ax[0]
    axe1 = ax[1]
    axe2 = axe0.twiny()
    axe0.plot(psal, depth, 'r', label='Practical salinity')
    axe2.plot(temp, depth, 'b', label='Temperature')
    axe1.plot(sigma_t, depth, 'g', label=u'$\sigma_t$')
    fig.suptitle(name.replace('_',' '))
    for axe in (axe0, axe1, axe2):
        axe.invert_yaxis()
        axe.legend()
        axe.set_ylabel('Depth (m)')
    axe0.set_xlabel('Salinity (psu)')
    axe2.set_xlabel(u'Temperature ($^{\circ}$C)')
    axe1.set_xlabel(u'Density anomaly $kg/m^3$')
    fig.savefig('Figures/Profiles/profile_' + name + '.pdf')
    fig.savefig('Figures/Profiles/profile_' + name + '.png')
    plt.close()

if __name__ == '__main__':
    prof('TB_20181210_03_down.nc')
