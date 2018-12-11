import xarray as xr
import gsw

def load_nc(filename):
    """
    Load data

    TODO: compile absolute salinity and potential temperature
    
    Arguments
    ---------
        filename: name of the file, string
            Assumes data are in 'Data/ctd_files/'
    """
    datadir = 'Data/ctd_files/'
    data = xr.open_dataset(datadir + filename)
    return data
