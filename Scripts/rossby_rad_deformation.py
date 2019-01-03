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
    lat = []
    lon = []
    temp = []
    sal = []
    depth = []
    for i in filenames:
        d = xr.open_dataset(os.path.join(datadir, i))
        temp.append(d.TEMP.values)
        sal.append(d.PSAL.values)
        lon.append(np.nanmean(d.lon.values))
        lat.append(np.nanmean(d.lat.values))
       
    depth = -np.abs(d.DEPTH.values)
    coords=lon
    #compute alpha and beta
    [rho,alpha,beta]=gsw.rho_alpha_beta(sal,temp,depth)

    g=9.8
    N = np.sqrt(g*((alpha[:,:-1]*(np.diff(temp)/np.diff(depth)))
    -(beta[:,:-1]*(np.diff(sal)/np.diff(depth))))) ##computed this way so that the influence of temperature 
    #vs salinity can be separated
    #plot N
    xx,yy=np.meshgrid(lon,depth[:-1])
    plt.pcolormesh(xx,yy,N.T)
    cb=plt.colorbar(extend='both')
    cb.set_label('N')
    
    #plt.show()
    
    plt.savefig('Figures/Transect/transect_N.pdf')
    plt.savefig('Figures/Transect/transect_N.png')
    plt.close()
    return N,np.asarray(rho),depth,lon,lat

def rossby_radii(datadir,filenames,desc=''):
    [N,rho,depth,lon, lat]=brunt_vaisala(datadir,filenames,desc='')
    
    #Length scale of rossby radius of deformation 
    def coriolis_acc(lat,omega=0.729*10**-4):
        return 2*omega*np.sin(lat/360.*2*np.pi)
    f=coriolis_acc(np.asarray(lat))
    Lr=(N[:]*depth[:-1]*-1)/(np.nanmean(f))*0.001

    plt.figure()
    xx,yy=np.meshgrid(lon,depth[:-1])
    plt.pcolormesh(xx,yy,Lr.T)
    cb=plt.colorbar(extend='both')
    plt.show()
    # cb.set_label('N')
    print('max:',np.nanmax(Lr))
    print('min:',np.nanmin(Lr))
    print('mean:',np.nanmean(Lr))

    #compare with width of fjord at this point ~ 3km

    #compute geostropic transport
    deldens=geos_vel(lon,lat,rho,f)

    # np.savetxt('Data/rossby_radius_deformation.txt',np.c_[, b_temp, a_sal, b_sal], \
    #            header="a_temp b_temp a_sal b_sal",comments='')

    return Lr

def haversine(lat1, lon1, lat2, lon2):
    import numpy as np
    R=6378.137
    lat1, lon1, lat2, lon2 = map(np.deg2rad, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    total_m = (R * c)*1000
    return total_m

def geos_vel(lon,lat,dens,f):
    """
    computes the geostrophic velocity across a section
    """
    dist=haversine(lon[1:],lat[1:],lon[:-1],lat[:-1])
    dist_ordered=dist.cumsum()/1000
    dist_r=(np.reshape(np.repeat(dist,119),(4,119)))

    deldens=np.diff(dens[:-1,:])/dist_r
    
    ref_dens=1000
    g=9.8

    geo_vert_shear=-((g/ref_dens)*(deldens))/np.nanmean(f)
    geo_vel=(np.nancumsum(geo_vert_shear,axis=1))
    geo_sta=np.nansum(geo_vel*dist_r*1)
    
    return (geo_sta/10**6)

# def steric_height():
#     """
#     Computes the steric height anomaly along the transect
#     This is aimed for the along fjord transect
#     """
#     alpha=gsw.specvol(sal,temp,depth)
#     alpha_ref=gsw.specvol(35,0,depth)
#     alpha_ano=alpha-alpha_ref

def section_along_density_contours():
    """
    plots section, but instead of being against depth, against density
    """
    



def load_filenames_section(N, sec_datadir='Data/sections'):
    """
    load and return the filenames + decription of the section

    data must be on a txt file like: secN.txt
    N must be replaced with the section number (e.g. 'sec1.txt', 'sec5.txt')
    """
    return np.loadtxt(os.path.join(sec_datadir, 'sec'+str(N)+'.txt'), \
                      dtype=str, delimiter='\n')
        

    
if __name__ == '__main__':
    N=2 #section number
    section_meta = load_filenames_section(N, 'Data/sections')
    if sys.argv[1] == 'brunt_vaisala':
        datadir = 'Data/ctd_files/gridded_calibrated'
        section_meta = load_filenames_section(0, 'Data/sections')
        brunt_vaisala(datadir, section_meta[1:], desc=section_meta[0])
    if sys.argv[1] == 'rossby_radii':
        datadir = 'Data/ctd_files/gridded_calibrated_updated'
        section_meta = load_filenames_section(2, 'Data/sections')
        rossby_radii(datadir, section_meta[1:], desc=section_meta[0])

