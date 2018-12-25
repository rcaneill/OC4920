

import numpy as np
from glob import glob
from osgeo import gdal
import airsea

flist = 'Data/winds/cdas*'
flist = sorted(glob(flist))

# u_wnd = np.ndarray([9, 361, 720])
# v_wnd = np.ndarray([9, 361, 720])
# net_wind = np.ndarray([9, 361, 720])
# tau = np.ndarray([9, 361, 720])

# Read the GRIB file

for i, fname in enumerate(flist):
    
    grib = gdal.Open(fname)

    if i < 1:
       
        # Read an specific band
        band = grib.GetRasterBand(266)
        # Read the band as a Python array
        u_wnd[i, :, :] = band.ReadAsArray()
   
        # Read an specific band
        band = grib.GetRasterBand(267)
        # Read the band as a Python array
        v_wnd[i, :, :] = band.ReadAsArray()
        
    if i > 0:
        
        # Read an specific band
        band = grib.GetRasterBand(288)
        # Read the band as a Python array
        u_wnd[i, :, :] = band.ReadAsArray()
        
        # Read an specific band
        band = grib.GetRasterBand(289)
        # Read the band as a Python array
        v_wnd[i, :, :] = band.ReadAsArray()
    
    x1 = np.square(u_wnd[i, :, :])
    x2 = np.square(v_wnd[i, :, :])

    x = x1+x2

    net_wind[i, :, :]=np.sqrt(x)
    
    tau[i, :, :] = airsea.windstress.stress(net_wind[i, :, :], z=10., drag='largepond', rho_air=1.22, Ta=10.)