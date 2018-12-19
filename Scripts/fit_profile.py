import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import xarray as xr

def calc_fit(depth, var):
    """
    computes the fit of the variable var along depth

    uses the function *multi_linear* for the fit
    """
    #p0 is the a priori fit for temperature
    p0=[0,20,60,80,5,10,9,7]
    popt, pcov = so.curve_fit(multi_linear_fjord, depth, var, p0=p0)
    #print([0,10,50,70,5,10,11,7])
    #print(popt)
    plt.plot(depth,var,'o')
    plt.plot(depth, [multi_linear(d, popt) for d in depth])
    plt.show()


def multi_linear_fjord(d_tot, d0, d1, d2, d3, v0, v1, v2, v3):
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

if __name__ == '__main__':
    data = xr.open_dataset('Data/ctd_files/gridded_calibrated/TB_20181210_10_down_grid.nc')
    depth=data.DEPTH.values
    temp=data.TEMP.values
    NAN = np.isnan(temp) | np.isnan(depth)
    depth = depth[~NAN]
    temp = temp[~NAN]
    calc_fit(depth, temp)
    #multi_linear_ex()
