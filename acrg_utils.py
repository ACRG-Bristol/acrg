


def is_number(s):
    """
    Is it a number?
    """
    try:
        float(s)
        return True
    except ValueError:
        return False
    

def combine_diff_resolution(data_low_res, data_high_res, method='multiply', verbose=True):
    '''
    Combine datasets which have different time resolutions
    
    Args
        data_low_res (xarray.DataArray or xarray.Dataset)
            data with the lower time resolution
        data_high_res (xarray.DataArray or xarray.Dataset)
            data with the higher time resolution
        method (str)
            method by which to combine the datasets
            defaults to 'multiply'
            can be:
                'multiply'      : data_low_res * data_high_res
                'add'           : data_low_res + data_high_res
                'subtract_high' : data_low_res - data_high_res
                'subtract_low'  : data_high_res - data_low_res
                'divide_by_high': data_high_res / data_low_res
                'divide_by_low' : data_low_res / data_high_res
        data_vars (dict)
            variable names for the high and low res data
            i.e. {'high_res': 'flux', 'low_res': 'flux'}
        verbose (bool)
    
    Returns
        xarray.DataArray or xarray.Dataset
            has the same shape as resolved_data
    '''
    # get the time step of the integrated data
    # used for grouping the high time res data
    time_step  = (data_low_res.time[1] - data_low_res.time[0]).values
    time_scale = 'Month' if time_step.astype('timedelta64[M]').astype(int)>0 else \
                 'Dayofyear' if time_step.astype('timedelta64[D]').astype(int)>0 else 'hour'
    time_step  = time_step.astype(f'timedelta64[{time_scale[0]}]').astype(int)
    
    if verbose:
        print(f'Method to combine datasets: {method}',
              f'\nIntegrated data has time step of {time_step} {time_scale.lower()}')
    
    # define a function to apply to the resolved data
    def combine_method(resolved, method):
        # select the correct time slice from the integrated data
        integrated = data_low_res.sel(time=resolved.time[0])
        
        if method.lower()=='multiply':
            return resolved * integrated
        elif method.lower()=='add':
            return resolved + integrated
        elif method.lower()=='subtract_high':
            return integrated - resolved
        elif method.lower()=='subtract_low':
            return resolved - integrated
        elif method.lower()=='divide_by_high':
            return interated / resolved
        elif method.lower()=='divide_by_low':
            return resolved / integrated
        else:
            print(f'Method not recognised')
    
    # group the data by the same time scale as the integrated da
    grouped_data = data_high_res.groupby(f'time.{time_scale.lower()}')
    
    # apply the combine method to each slice of the grouped data
    output = grouped_data.map(combine_method, method=method)
    
    return output