import numpy as np
import xarray as xr
import scipy.stats as ss
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams

plt.style.use('seaborn')
rcParams['text.usetex'] = True


def linregress(x,y):
    mask = np.isfinite([x, y]).all(axis=0)
    Yclean = y[mask]
    Xclean = x[mask]
    return ss.linregress(Xclean,Yclean)

def corr_coef(datadir,SK_filename,TB_filename):
    """
    computes the coef for the linear fit between Trygve and Skagerak

    Return *a_temp* and *b_temp* like: Trygve_new = a_temp*Trygve + b_temp
    We consider Skagerak as best values
    """
    TB = xr.open_dataset(os.path.join(datadir,TB_filename))
    SK = xr.open_dataset(os.path.join(datadir,SK_filename))
    #print(TB.TEMP)
    #print(SK.TEMP)
    # Temp
    (a_temp, b_temp, r_value, p_value, std_err) = linregress(SK.TEMP,TB.TEMP)
    plt.scatter(SK.TEMP,TB.TEMP)
    plt.plot([np.nanmin(SK.TEMP),np.nanmax(SK.TEMP)], \
             [a_temp*np.nanmin(SK.TEMP)+b_temp,a_temp*np.nanmax(SK.TEMP)+b_temp], \
             label='a={}\nb={}'.format(a_temp,b_temp))
    plt.xlabel('Skagerak')
    plt.ylabel('Trygve')
    plt.legend()
    plt.title(u'Temperature ($^{\circ}$C)') 
    plt.savefig('Figures/Calib/temp.png')
    plt.savefig('Figures/Calib/temp.pdf')
    plt.close()

    # Sal
    (a_sal, b_sal, r_value, p_value, std_err) = linregress(SK.PSAL,TB.PSAL)
    plt.scatter(SK.PSAL,TB.PSAL)
    plt.plot([np.nanmin(SK.PSAL),np.nanmax(SK.PSAL)], \
             [a_sal*np.nanmin(SK.PSAL)+b_sal,a_sal*np.nanmax(SK.PSAL)+b_sal], \
             label='a={}\nb={}'.format(a_sal,b_sal))
    plt.xlabel('Skagerak')
    plt.ylabel('Trygve')
    plt.legend()
    plt.title(u'Salinity (Psu)')  
    plt.savefig('Figures/Calib/sal.png')
    plt.savefig('Figures/Calib/sal.pdf')
    plt.close()

if __name__ == "__main__":
    corr_coef('Data/ctd_files/gridded', 'SK_20181210_Calibration_grid.nc', \
              'TB_2018121cal_down_grid.nc')
    #corr_coef('Data/ctd_files/gridded', 'SK_20181211_01_grid.nc', \
    #          'TB_20181211_cal_down_grid.nc')
