import numpy as np
import scipy.io as scio
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
import gsw
from matplotlib import rcParams
import os, sys
import pandas as pd
from matplotlib.image import NonUniformImage
from scipy.interpolate import griddata
from scipy.stats import linregress

# plt.style.use('seaborn')
rcParams['text.usetex'] = True

#calibrate SK to TB, given SK has oxygen optode issues
import xarray as xr
datadir='Data/ctd_files/gridded_calibrated_updated'
SK_filename='SK_20181211_01_grid.nc'
TB_filename='TB_20181211_01_grid.nc'
SK=xr.open_dataset(os.path.join(datadir,SK_filename))
TB=xr.open_dataset(os.path.join(datadir,TB_filename))
mask = ~np.isnan(SK['Sbeox0MLperL']) & ~np.isnan(TB['Sbeox0MLperL'])

[a_oxy,b_oxy,r_value,p_value,stderr]=linregress(SK['Sbeox0MLperL'][mask][15:-1],TB['Sbeox0MLperL'][mask][15:-1])

mn=np.nanmin(TB.Sbeox0MLperL[mask][15:-1])
mx=np.nanmax(SK.Sbeox0MLperL[mask][15:-1])
x1=np.linspace(mn,mx,500)
y1=a_oxy*x1+b_oxy

fig=plt.figure()
plt.text(7,2.4,'r$^2$={}'.format(np.round(r_value**2,2)))
plt.scatter(SK.Sbeox0MLperL[mask][15:-1],TB.Sbeox0MLperL[mask][15:-1])

plt.plot(x1,y1,'-r',label='a={}\nb={}'.format(a_oxy,b_oxy))
plt.legend()
plt.xlim(np.nanmin(SK.Sbeox0MLperL[mask][15:-1]-0.1),np.nanmax(SK.Sbeox0MLperL[mask][15:-1]+0.1))
plt.ylim((np.nanmin(TB.Sbeox0MLperL[mask][15:-1]-0.1),np.nanmax(TB.Sbeox0MLperL[mask][15:-1]+0.1)))

plt.ylabel('O$_2$ Trygve [ml/l]')
plt.xlabel('O$_2$ Skagerak [ml/l]')
plt.savefig('Figures/Calib/oxygen_SK_TB.png')
plt.close()


#######Compute new O2 values for SK

ls = os.listdir(datadir)
ls.sort()
for filename in ls:
    if filename[:11] == 'SK_20181211':
        # print(filename)
        data =  xr.open_dataset(os.path.join(datadir,filename))
        varkey = [i for i in data.data_vars.keys()]
        if 'lon' in varkey:
            print(filename)
            # mask = ~np.isnan(data.Sbeox0MLperL)

            data['oxy_corrected']=('DEPTH', \
                    (a_oxy*data.Sbeox0MLperL+b_oxy))
        data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated/oxygen_step1',filename))

ls = os.listdir(datadir)
ls.sort()
for filename in ls:
    if filename[:11] == 'SK_20181210':
        # print(filename)
        data =  xr.open_dataset(os.path.join(datadir,filename))
        varkey = [i for i in data.data_vars.keys()]
        if 'lon' in varkey:
            print(filename)
            # mask = ~np.isnan(data.Sbeox0MLperL)

            data['oxy_corrected']=('DEPTH', \
                    (a_oxy*data.Sbeox0MLperL+b_oxy))
        data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated/oxygen_step1',filename))


#########calibrate new o2 values to titrations
datadir = "Data/oxygen"
filename='o2_calib_data.csv'

dat=pd.read_csv(os.path.join(datadir,filename))

mask = ~np.isnan(dat.cO2) & ~np.isnan(dat.o2corr)

[a_oxy,b_oxy,r_value,p_value,stderr]=linregress(dat.o2corr.loc[mask],dat.cO2.loc[mask])

print(r_value**2)

oxy_corrected=(a_oxy*dat.o2corr+(b_oxy))
mn=np.nanmin(dat.cO2[0:10])
mx=np.nanmax(dat.cO2[0:10])
x1=np.linspace(mn,mx,500)
y1=a_oxy*x1+b_oxy

fig=plt.figure()
plt.text(7,2.4,'r$^2$={}'.format(np.round(r_value**2,2)))
plt.scatter(dat.cO2[0:10],dat.o2corr[0:10],label='O$_2$ optode ')

plt.plot(x1,y1,'-r',label='a={}\nb={}'.format(a_oxy,b_oxy))
plt.legend()

plt.ylabel('c(O$_2$) [ml/l]')
plt.xlabel('O$_2$ optode [ml/l]')
plt.savefig('Figures/Calib/oxygen.png')
plt.close()

####resave with corrected oxygen values
datadir='Data/ctd_files/gridded_calibrated_updated/oxygen_step1'

ls = os.listdir(datadir)
ls.sort()
for filename in ls:
    if filename[:11] == 'SK_20181211':
        # print(filename)
        data =  xr.open_dataset(os.path.join(datadir,filename))
        varkey = [i for i in data.data_vars.keys()]
        if 'lon' in varkey:
            print(filename)
            # mask = ~np.isnan(data.oxy_corrected)

            data['oxy_calibrated']=('DEPTH', 
                    (a_oxy*data.oxy_corrected+(b_oxy)))
        data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated/oxygen_step2',filename))

ls = os.listdir(datadir)
ls.sort()
for filename in ls:
    if filename[:11] == 'SK_20181210':
        # print(filename)
        data =  xr.open_dataset(os.path.join(datadir,filename))
        varkey = [i for i in data.data_vars.keys()]
        if 'lon' in varkey:
            print(filename)
            # mask = ~np.isnan(data.oxy_corrected)

            data['oxy_calibrated']=('DEPTH', 
                    (a_oxy*data.oxy_corrected+(b_oxy)))
        data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated/oxygen_step2',filename))

###repeat for trygve values
datadir='Data/ctd_files/gridded_calibrated_updated/'

ls = os.listdir(datadir)
ls.sort()
for filename in ls:
    if filename[:11] == 'TB_20181211':
        # print(filename)
        data =  xr.open_dataset(os.path.join(datadir,filename))
        varkey = [i for i in data.data_vars.keys()]
        if 'lon' in varkey:
            print(filename)
            # mask = ~np.isnan(data.Sbeox0MLperL)

            data['oxy_calibrated']=('DEPTH', 
                    (a_oxy*data.Sbeox0MLperL+(b_oxy)))
        data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated/oxygen_step2',filename))

ls = os.listdir(datadir)
ls.sort()
for filename in ls:
    if filename[:11] == 'TB_20181210':
        # print(filename)
        data =  xr.open_dataset(os.path.join(datadir,filename))
        varkey = [i for i in data.data_vars.keys()]
        if 'lon' in varkey:
            print(filename)
            # mask = ~np.isnan(data.Sbeox0MLperL)

            data['oxy_calibrated']=('DEPTH', 
                    (a_oxy*data.Sbeox0MLperL+(b_oxy)))
        data.to_netcdf(os.path.join('Data/ctd_files/gridded_calibrated_updated/oxygen_step2',filename))

###compare titrated values with recalulated values on sk




# plt.show()
# [a_oxy,b_oxy,r_value,p_value,stderr]=linregress(dat.cO2[0:10],dat['CTDO2'][0:10])
# [a_oxy,b_oxy,r_value,p_value,stderr]=linregress(dat['CTDO2'][0:10],dat.cO2[0:10])

# print(a_oxy,b_oxy)
# print(r_value**2)
# # mat = scio.loadmat(os.path.join(datadir,filename))
# mn=np.nanmin(dat['CTDO2'][0:10])
# mx=np.nanmax(dat['CTDO2'][0:10])
# x1=np.linspace(mn,mx,500)
# y1=a_oxy*x1+b_oxy

# fig=plt.figure()
# # plt.scatter(dat.cO2[0:10],dat['CTDO2'][0:10],label='O$_2$ optode ')
# plt.scatter(dat['CTDO2'][0:10],dat.cO2[0:10],label='O$_2$ optode ')

# plt.plot(x1,y1,'-r',label='a={}\nb={}'.format(a_oxy,b_oxy))
# # plt.text(0.1,0.9,'r$^2$={}'.format(r_value**2))

# # plt.legend()

# plt.ylabel('c(O$_2$) [ml/l]')
# plt.xlabel('O$_2$ optode [ml/l]')
# plt.savefig('Figures/Calib/oxygen.png')

# plt.show()



# #correct oxygen SK to titration values
# oxy_corrected=(a_oxy*data.Sbeox0MLperL+(b_oxy))
# plt.plot(data.Sbeox0MLperL)
# plt.plot(oxy_corrected)
# plt.show()

# #correct 

# ##
# # ls = os.listdir(datadir)
# # ls.sort()
# # for filename in ls:
# #     if filename[:11] == 'TB_20181211':
# #             # print(filename)
# #         data =  xr.open_dataset(os.path.join(datadir,filename))
# #         varkey = [i for i in data.data_vars.keys()]
# #         if 'lon' in varkey:
# #             print(filename)
# #             data['oxy_corrected']=('DEPTH', \
# #                 (coeffs.a_oxy.values*data.TEMP+(coeffs.b_oxy.values)))




# # Oxy

# # (a_oxy, b_oxy, r_value, p_value, std_err) = linregress(TB.Sbeox0MLperL[15:],SK.Sbeox0MLperL[15:])
# # plt.scatter(TB.Sbeox0MLperL[15:],SK.Sbeox0MLperL[15:])
# # plt.plot([np.nanmin(SK.Sbeox0MLperL),np.nanmax(SK.Sbeox0MLperL)], \
# #                 [a_oxy*np.nanmin(TB.Sbeox0MLperL[:])+b_oxy,a_oxy*np.nanmax(TB.Sbeox0MLperL[:])+b_oxy], \
# #                 label='a={}\nb={}'.format(a_oxy,b_oxy))
# # plt.xlabel('Trygve')
# # plt.ylabel('Skagerak')
# # plt.legend()
# # plt.title(u'Oxygen (ml/l)')  
# # plt.savefig('Figures/Calib/oxy_20181211.png')
# # plt.savefig('Figures/Calib/oxy_20181211.pdf')
# # plt.close()

# # oxy_SK1=mat['oxy_SK1']
# # print(oxy_SK1.shape)
# # # oxy_SK2=mat['oxy_SK2']
# # # true_oxy_SK10=mat['true_oxy_SK10']
# # # true_oxy_SK11=mat['true_oxy_SK11']

# # # plt.scatter(oxy_SK1,true_oxy_SK10)
# # # plt.show()

# # filename='oxy_SK1.mat'
# # mat = scio.loadmat(os.path.join(datadir,filename))
# # print(mat.keys())
# # oxy_SK1=mat['oxy_SK1']
# # print(oxy_SK1.shape)