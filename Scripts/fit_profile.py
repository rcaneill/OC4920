import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import xarray as xr
import os

def misfit(layer, depth, temp):
    #print(len(multi_linear_fjord(depth,*layer)))
    #print(len(depth))
    #print(len(temp))
    return np.sum((multi_linear_fjord(depth,*layer)-temp)**2)

def calc_fit(depth, var, plotBool=False):
    """
    computes the fit of the variable var along depth

    uses the function *multi_linear* for the fit

    Set to NaN all parameters that are extrapolated from where there are no data
    """
    depth = np.array(depth)
    var = np.array(var)
    NAN = np.isnan(depth) | np.isnan(var)
    depth = depth[~NAN]
    var = var[~NAN]
    #p0 is the a priori fit for temperature
    p0=np.array([0,20,60,80,5,10,9,7])
    #p0=np.array([0,20,60,50,5,10,9,7.5])
    #p0=np.ones(8)
    # bounds are set by looking at all transects and estimate by eye
    # first 4, values are depth, last 4 values are temp
    bounds_low = [0, 4, 55, 60, 2, 9, 9, 6]
    bounds_up = [6, 30, 65, 80, 7, 11, 10, 8.5]
    #p0 = np.array(bounds_low)
    #bounds_low = [0, 4, 55, 50, 2, 9, 9, 0]
    #bounds_up = [6, 30, 65, 80, 7, 11, 10, 8.5]
    bounds = (bounds_low, bounds_up)
    #bounds=(-np.inf,np.inf)
    try:
        #popt, pcov = so.curve_fit(multi_linear_fjord, depth, var, p0=p0, bounds=bounds, \
        #                          verbose=2, method='trf', loss='cauchy')
        result = so.minimize(misfit, x0=p0, args=(depth,var), \
                             bounds=[(i,j) for (i,j) in zip(bounds_low, bounds_up)])
        popt = result.x
        layer_d = popt[:int(len(popt)/2)]
        layer_v = popt[int(len(popt)/2):]
        bool_d = layer_d < np.nanmax(depth) # Do we trust the data
        popt = np.concatenate((layer_d,layer_v))
        print('Resulting values for layer depth, temperature and Trust flag')
        print(layer_d)
        print(layer_v)
        print(bool_d)
    except RuntimeError:
        print("Fit cannot be computed, all layer values set to NaN")
        return (p0*np.NaN, [False, False, False, False])
    #print([0,10,50,70,5,10,11,7])
    #print(popt)
    if plotBool:
        plt.plot(var, depth,'o')
        plt.plot([multi_linear(d, popt) for d in range(120)], range(120))
        plt.ylim(119,0)
        plt.show()
    return (popt, bool_d)


def multi_linear_fjord(d_tot, d0, d1, d2, d3, v0, v1, v2, v3):
    #print([d0, d1, d2, d3, v0, v1, v2, v3])
    return [multi_linear(d, layer=[d0, d1, d2, d3, v0, v1, v2, v3]) for d in d_tot]
    
def multi_linear(d, layer):
    """
    layer contains depth and then values
    so that :
    layer_depth = layer[:int(len(layer)/2)]
    layer_value = layer[int(len(layer)/2):]
    
    returns the value of the studied variable, at the depth *d*
    
    the function evolves linearly between the different depths
    first value of layer_depth must be 0 (i.e. surface value),
    layers_depth must be an increasing array
    
    if d > max(layer_depth), return value of layer_value at this point
    
    run function multi_linear_ex for an example
    """
    layer_depth=np.array(layer[:int(len(layer)/2)])
    layer_value=np.array(layer[int(len(layer)/2):])
    # removing NaN
    NAN = np.isnan(layer_depth) | np.isnan(layer_value)
    # if only NaN, we return NaN
    if NAN.all():
        print('Returning NaN for mutli linear fit')
        return np.nan
    layer_depth = layer_depth[~NAN]
    layer_value = layer_value[~NAN]
    # sorting to get depth increasing
    layer_value = layer_value.ravel()[layer_depth.argsort(axis=None).reshape(layer_depth.shape)]
    layer_depth.sort()
    # find inbetween values for d
    # these depth values are called d0 and d1
    # idem with values, v0 and v1
    if d <= np.min(layer_depth):
        return layer_value[0]
    if d >= np.max(layer_depth):
        return layer_value[-1]
    else:
        d0 = layer_depth[d >= layer_depth][-1]
        d1 = layer_depth[d <= layer_depth][0]
        v0 = layer_value[d >= layer_depth][-1]
        v1 = layer_value[d <= layer_depth][0]
    if d0 == d1:
        return v0
    else:
        return v0 + (d-d0)*(v1-v0)/(d1-d0)

def multi_linear_ex():
    """
    plots an example of multi_linear function
    """
    d_tot = np.linspace(0,100,101)
    layer_depth = [0,10,50,70]
    layer_value = [5,10,11,7]
    layer = layer_depth + layer_value
    #print(layer)
    v_tot = [multi_linear(d, layer) for d in d_tot]

    plt.plot(v_tot, d_tot)
    plt.plot(layer_value, layer_depth, 'o')
    plt.ylim(d_tot[-1], d_tot[0])
    plt.show()

def fit_all_prof(datadir, plotBool=False):
    """
    Fit all temperature profile of data in datadir
    save the data on the same netcdf file, with new attributes:
    data.TEMP_f : for the fit of temperature
    data.layer_d : the depth of the different layers
    data.layer_t : the temperature at the different depth of layers
    """
    for filename in os.listdir(datadir):
        if filename=='SK_20181210_07_grid.nc' or True:
            print(filename)
            data = xr.open_dataset(os.path.join(datadir, filename))
            (layer,layer_flag) = calc_fit(data.DEPTH.values, data.ptemp.values, \
                                          plotBool=plotBool)
            layer_d = layer[:int(len(layer)/2)]
            layer_v = layer[int(len(layer)/2):]
            TEMP_f = [multi_linear(d, layer) for d in data.DEPTH.values]
            data['TEMP_f'] = ('DEPTH', TEMP_f)
            data['layer_d'] = layer_d
            data['layer_v'] = layer_v
            data['layer_flag'] = layer_flag
            os.remove(os.path.join(datadir, filename))
            data.to_netcdf(os.path.join(datadir, filename))

        
if __name__ == '__main__':
    data = xr.open_dataset('Data/ctd_files/gridded_calibrated_updated/TB_20181210_10_grid.nc')
    depth=data.DEPTH.values
    temp=data.TEMP.values
    NAN = np.isnan(temp) | np.isnan(depth)
    depth = depth[~NAN]
    temp = temp[~NAN]
    #calc_fit(depth, temp)
    #multi_linear_ex()
    fit_all_prof('Data/ctd_files/fitted', plotBool=False)
