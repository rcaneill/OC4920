import numpy as np
import xarray as xr
import pandas as pd
import scipy.stats as ss
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams
import gsw

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
    (a_temp, b_temp, r_value, p_value, std_err) = linregress(TB.TEMP[15:],SK.TEMP[15:])
    print(a_temp,b_temp)
    plt.scatter(TB.TEMP[15:],SK.TEMP[15:])
    plt.plot([np.nanmin(SK.TEMP),np.nanmax(SK.TEMP)], \
             [a_temp*np.nanmin(TB.TEMP[:])+b_temp,a_temp*np.nanmax(TB.TEMP[:])+b_temp], \
             label='a={}\nb={}'.format(a_temp,b_temp))
    plt.xlabel('Trygve')
    plt.ylabel('Skagerak')
    plt.legend()
    plt.title(u'Temperature ($^{\circ}$C)') 
    plt.savefig('Figures/Calib/temp_20181211.png')
    plt.savefig('Figures/Calib/temp_20181211.pdf')
    plt.close()

    # Sal
    (a_sal, b_sal, r_value, p_value, std_err) = linregress(TB.PSAL[15:],SK.PSAL[15:])
    plt.scatter(TB.PSAL[15:],SK.PSAL[15:],)
    plt.plot([np.nanmin(SK.PSAL),np.nanmax(SK.PSAL)], \
             [a_sal*np.nanmin(TB.PSAL[:])+b_sal,a_sal*np.nanmax(TB.PSAL[:])+b_sal], \
             label='a={}\nb={}'.format(a_sal,b_sal))
    plt.xlabel('Trygve')
    plt.ylabel('Skagerak')
    plt.legend()
    plt.title(u'Salinity (Psu)')  
    plt.savefig('Figures/Calib/sal_20181211.png')
    plt.savefig('Figures/Calib/sal_20181211.pdf')
    plt.close()

    np.savetxt('Data/calib_ts_20181211.txt',np.c_[a_temp, b_temp, a_sal, b_sal], \
               header="a_temp b_temp a_sal b_sal",comments='')
    return (a_temp, b_temp, a_sal, b_sal)

def convert_ts(datadir,filename):
    """
    Function to convert Practical salinity to Absolute salinity 
    and Conservative Temperature to potential temperature
    """
    ls = os.listdir(datadir)
    ls.sort()
    for filename in ls:
        if filename[:11] == 'SK_20181211':
            # print(filename)
            data =  xr.open_dataset(os.path.join(datadir,filename))
            print(filename)
            data['ab_sal']=gsw.SA_from_SP(data.PSAL,data.DEPTH,data.lon.values,data.lat.values)
            data['ptemp']=gsw.pt_from_CT(data.ab_sal,data.TEMP)
            data['ab_sal_bal']=gsw.SA_from_SP_Baltic(data.PSAL,data.lon.values,data.lat.values)
            data['ptemp_bal']=gsw.pt_from_CT(data.ab_sal_bal,data.TEMP)

            data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated',filename))



def correct_ts(datadir,filename,corr_coeff_filename):
    """
    Corrects TB Temp and salinity according to correlation coefficients and adds corrected 
    var to ncfile
    Also converts practical salinity to absolute salinity and conservative temperature to 
    potential temperature
    """
    coeffs=pd.read_csv(corr_coeff_filename,delim_whitespace=True)
    ls = os.listdir(datadir)
    ls.sort()
    for filename in ls:
        if filename[:11] == 'TB_20181210':
            # print(filename)
            data =  xr.open_dataset(os.path.join(datadir,filename))
            print(filename)
            data['t_corrected']=(coeffs.a_temp.values*data.TEMP+(coeffs.b_temp.values))
            data['s_corrected']=(coeffs.a_sal.values*data.PSAL+(coeffs.b_sal.values))
            data['ab_sal']=gsw.SA_from_SP(data.s_corrected,data.DEPTH,data.lon.values,data.lat.values)
            data['ptemp']=gsw.pt_from_CT(data.ab_sal,data.TEMP)
            data['ab_sal_bal']=gsw.SA_from_SP_Baltic(data.s_corrected,data.lon.values,data.lat.values)
            data['ptemp_bal']=gsw.pt_from_CT(data.ab_sal_bal,data.TEMP)

            data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated',filename))
    
    # Trygve_new = a_temp*Trygve + b_temp

if __name__ == "__main__":
    # corr_coef('Data/ctd_files/gridded', 'SK_20181210_Calibration_grid.nc', \
    #           'TB_2018121cal_down_grid.nc')
    #  corr_coef('Data/ctd_files/gridded', 'SK_20181211_01_grid.nc', \
    #        'TB_20181211_cal_down_grid.nc')
    correct_ts('Data/ctd_files/gridded','TB_20181210*.nc','Data/calib_ts_20181210.txt')
    # convert_ts('Data/ctd_files/gridded','SK_20181211*.nc')
