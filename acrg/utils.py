import gzip
import bz2
import os 
import shutil
import json

from acrg.config.paths import Paths

acrg_path = Paths.acrg

def decompress_files(filepaths):
    """ Decompress files at paths in filepaths

        Args:
            filepaths (list): List of paths to files to decompress
        Returns:
            list: List of paths to decompressed files
    """
    _, extension = os.path.splitext(filepaths[0])

    if extension == ".bz2":
        complib = bz2
    elif extension == ".gz":
        complib = gzip
    else:
        raise ValueError("Unable to decompress files other than bz2 and gz")
        
    decompressed_files = []
   
    for compressed_path in filepaths:
        # Here decompressed_path will be the path minus the extension
        decompressed_path, extension = os.path.splitext(compressed_path)

        with complib.open(compressed_path, "rb") as f_in, open(decompressed_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

        decompressed_files.append(decompressed_path)

    return decompressed_files


def is_number(s):
    """
    Is it a number?
    """
    try:
        float(s)
        return True
    except ValueError:
        return False

def combine_diff_resolution(data_1, data_2, method='add', verbose=True):
    '''
    Combine datasets which have different time resolutions
    Avoids having to resample data to the same resolution in order to add,
    subtract, multiply, or divide them
    
    Args
        data_1, data_2 (xarray.DataArray or xarray.Dataset)
            Data to combine
            Must have resolution of 1 hour, 1 day, 1 month, or 1 year
        method (str)
            Method by which to combine the datasets
            Can be:
                'add':      data_1 + data_2 (default)
                'subtract': data_1 - data_2
                'multiply': data_1 * data_2
                'divide':   data_1 / data_2
    
    Returns
        xarray.DataArray or xarray.Dataset
            has the same shape as resolved_data
    '''
    if 'np' not in dir(): import numpy as np
    data = {dd: dat for dd, dat in enumerate([data_1, data_2])}
    # calculate the time step for each dataset
    time_step  = {dd: (dat.time[1] - dat.time[0]).values for dd, dat in data.items()}
    
    if time_step[0]==time_step[1]:
        # if the time steps are equal then apply the method as normal
        if method.lower()=='multiply':
            return data_1 * data_2
        elif method.lower()=='add':
            return data_1 + data_2
        elif method.lower()=='subtract':
            return data_1 - data_2
        elif method.lower()=='divide':
            return data_1 / data_2
        else:
            print(f'Method not recognised: {method}')
            return(None)
    
    # work out which dataset is resolved and which is integrated
    data_val = {'resolved'  : [dd for dd, step in time_step.items() if step==np.nanmin(np.array(list(time_step.values())))][0],
                'integrated': [dd for dd, step in time_step.items() if step==np.nanmax(np.array(list(time_step.values())))][0]}
    data = {res: data[dd] for res, dd in data_val.items()}
    
    # work out the time scale for each dataset
    time_scale = {res: 'Month' if time_step[dd].astype('timedelta64[M]').astype(int)>0 else 
                       'Dayofyear' if time_step[dd].astype('timedelta64[D]').astype(int)>0 else 'hour'
                  for res, dd in data_val.items()}
    time_step  = {res: time_step[dd].astype(f'timedelta64[{time_scale[res][0]}]').astype(int)
                  for res, dd in data_val.items()}
    
    if verbose:
        print(f'Method to combine datasets: {method}',
              f'\nIntegrated data has time step of {time_step["integrated"]} {time_scale["integrated"].lower()}')
    
    # function to apply to the resolved data
    def combine_method(resolved, method):
        # select the correct time slice from the integrated data
        integrated = data['integrated'].sel(time=resolved.time[0])
        
        if method.lower()=='multiply':
            return resolved * integrated
        elif method.lower()=='add':
            return resolved + integrated
        elif method.lower()=='subtract_high':
            return integrated - resolved
        elif method.lower()=='subtract_low':
            return resolved - integrated
        elif method.lower()=='divide_low':
            # divide by the low res data
            return resolved / integrated
        elif method.lower()=='divide_high':
            return interated / resolved
        else:
            print(f'Method not recognised: {method}')
            return(None)
    
    # group the data by the same time scale as the integrated da
    grouped_data = data['resolved'].groupby(f'time.{time_scale["integrated"].lower()}')
    
    # the required method depends on which order the data was given in
    method = method if method in ['multiply', 'add'] else \
             '_'.join([method, 'high']) if (data['integrated']==data_1).all() else \
             '_'.join([method, 'low']) if (data['integrated']==data_2).all() else \
             method
    
    # apply the combine method to each slice of the grouped data
    output = grouped_data.map(combine_method, method=method)
    
    return output

def load_json(filename):
    """Returns a dictionary deserialised from JSON.

    Args:
        filename: Name of JSON file
    Returns:
        dict: Dictionary created from JSON
    """
    from json import load

    with open(filename, "r") as f:
        data: Dict[str, Any] = load(f)

    return data