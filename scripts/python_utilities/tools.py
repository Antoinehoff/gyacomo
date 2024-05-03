import fourier
import numpy as np

def zkxky_to_xy_const_z(array, iz):
    # Get shape of the phi array
    Nz, Nkx, Nky, = array.shape
    if iz < 0: #outboard midplane for negative iz
        iz = Nz // 2  # Using the middle value for z

    array = array[iz,:,:]
    array = fourier.kx_to_x(array,Nkx,-2)
    array = fourier.ky_to_y(array,Nky-1,-1)
    array = np.transpose(array)
    array = np.flip(np.fft.fftshift(array))
    return array
    
def closest_index(array, v):
    # Compute absolute differences between each element of the array and v
    absolute_diff = np.abs(array - v)
    
    # Find the index of the minimum difference
    closest_index = np.argmin(absolute_diff)
    closest_index = max(closest_index,1)
    return closest_index

def is_convertible_to_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def numpy_array_to_list(d):
    """
    Recursively convert NumPy arrays to lists within a dictionary.
    """
    for key, value in d.items():
        if isinstance(value, np.ndarray):
            d[key] = value.tolist()
        elif isinstance(value, dict):
            d[key] = numpy_array_to_list(value)
    return d