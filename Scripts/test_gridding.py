import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt

def test_gridding(datadir_raw, datadir_grid):
    """
    test gridded data to see if some have depth changes

    datadir is datadir of nc data
    """
    for filename in os.listdir(datadir_raw):
        data = xr.open_dataset(os.path.join(datadir_raw, filename))
        #data_grid = xr.open_dataset(os.path.join(datadir_raw, filename+'.bak'))
        data_grid = xr.open_dataset(os.path.join(datadir_grid, filename[:-3]+'_grid.nc'))
        print(data.DEPTH.values)
        #raise NotImplementedError
        plt.plot(data.TEMP, data.DEPTH, label='raw')
        plt.plot(data_grid.TEMP, data_grid.DEPTH, label='grid')
        plt.legend(title=filename)
        plt.show()
        #raise NotImplementedError
        plt.close()

if __name__ == "__main__":
    test_gridding('Data/ctd_files/processed2nc/Skagerak', \
                  'Data/ctd_files/gridded')
