from collections import OrderedDict

import numpy as np
import xarray as xray
import acrg.obs as acrg_obs
import datetime as dt
import pandas as pd

import glob
import numpy as np
import re
import random
import logging
import os
import itertools
import math
import glob

from acrg.config.paths import Paths
data_path = Paths.data

home = os.getenv("HOME")
input_directory=os.path.join(data_path,"obs_raw/GOSAT/CH4_GOS_OCPR_v7.2/")
fp_directory = os.path.join(data_path,'LPDM/fp_NAME/')
obs_directory = os.path.join(data_path,'obs/') # Where to write output nc files
name_csv_directory = os.path.join(home,"NAME_files") # Where to write output NAME csv files
name_pressure_directory = os.path.join(data_path,"LPDM/surface_pressure/")

logger = logging.getLogger("acrg.satellite.common")

def add_coords(ds,network=None,data_vars=[]):
    '''
    The gosat_add_coords function adds co-ordinates to a GOSAT dataset.
    This is valid for (at least) v6 of GOSAT CH4 Proxy Level 2 Data Product downloaded from the CCI Open Data Portal
    (e.g. filename of the form: ESACCI-GHG-L2-CH4-GOSAT-OCPR-20101231-fv6.nc)
    
    This assumes the input Dataset will be of the form:
        <xarray.Dataset>
        Dimensions:                           (m: 20, n: 1962)
        Dimensions without coordinates: m, n
        Data variables:
            xch4_quality_flag                 (n) int8 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ...
            ch4_profile_apriori               (n, m) float32 1751.95 1751.95 1751.94 ...
            xch4                              (n) float32 1644.78 1657.72 1648.04 ...
            ...
        Attributes:
            title:                     ESA CCI GOSAT OCPR CH4
            institution:               University of Leicester (UoL), UK
            ...
    
    i.e. with no coordinates assigned for the two dimensions within the Dataset.
    This function redefines "m" as "time" and "n" as "level" and adds the associated values as coords.

    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file opened with xarray as a Dataset.
        data_vars (list) : 
            Data variables which will have the time and level (where appropriate) dataset coords added
            If data_vars if left blank ([] by default) coords will be added to all data variables within ds.
                              
    Returns:
        ds (xarray.Dataset) : 
            Dataset with time and level coordinates assigned
            (only includes specified data variables unless left blank)
    '''
    if network=='TCCON':
        dims_current = ['time', 'prior_altitude']
        dims = ['time','lev']
    elif network=='GOSAT':
        dims_current = ['n','m']
        dims = ['time','lev']
    else:
        raise ValueError('Network should be specified as one of "TCCON" or "GOSAT" when calling this function.')
    
    if list(ds.dims.keys()) != dims_current and list(ds.dims.keys()) != dims_current[::-1]: # Check dimensions match what we expect (check forward and reverse list)
        logger.warning("Do not recognise dimensions of input tccon dataset. Unable to add new dimensions.")
        return None
    
    # Define co-ordinate values. Extract from Dataset if already present (e.g. time), construct if not (e.g. lev)
    coords = {}
    for i,d in enumerate(dims):
        try:
            coords[d] = ds[d].values # Extract values for each coord from Dataset if present (e.g. time)
        except (AttributeError,KeyError):
            current_dim = dims_current[i]
            dim_length = ds.dims[current_dim]
            coords[d] = np.arange(dim_length) # If no values defined, create integer list based on dimension size
    
    if not data_vars:
        data_vars = ds.data_vars
    
    data = OrderedDict()
    for name in data_vars:
        if name not in dims:
            if len(ds[name].dims) == 1:
                if ds[name].values.size == coords[dims[0]].size:
                    data[name] = (dims[0],ds[name].values)
                elif ds[name].values.size == coords[dims[1]].size:
                    data[name] = (dims[1],ds[name].values)
                else:
                    raise ValueError("Cannot find a dimension of matching size.")
            elif len(ds[name].dims) == 2:
                data[name] = ([dims[0],dims[1]],ds[name].values)
    
    ds_new = xray.Dataset(data,coords=coords,attrs=ds.attrs)
    for dv in ds_new.data_vars:
        ds_new[dv].attrs = ds[dv].attrs
    for coord in ds_new.coords:
        try:
            ds_new[coord].attrs = ds[coord].attrs
        except (AttributeError,KeyError):
            if coord == "lev" or coord == "level":
                ds_new[coord].attrs["short_description"] = "Number for each level within the vertically resolved data."
                ds_new[coord].attrs["long_name"] = "level"
    
    # Add attribute describing modification made to original data
    mod_attr = "nc_restructure"
    # ds_new = add_history_attr(ds_new,mod_attr)
    ds_new.attrs[mod_attr] = "Dataset rearranged to re-assign {} dimensions as {} coordinates, extracting values from dataset where present.".format(dims_current,dims).replace("'","")
    
    return ds_new


def add_history_attr(ds,mod_attr):
    '''
    The add_history_attr adds or modifies the "history" attributes to include details of the modifications
    made to the original data.
    When history attribute is first created a timestamp is added of the current time.
    Note: Expect an additional attribute based on mod_attr has been added which contains details of the
    modification.
    
    Args:
        ds (xarray.Dataset) :
            Input Dataset.
        mod_attr (str) :
            Name of attribute which includes details of modification. This string will be added to the
            history attribute.
    Returns:
        xarray.Dataset:
            Input dataset with history attribute added or modified.
    '''
    if "history" in ds.attrs:
        ds.attrs["history"] += "{}, ".format(mod_attr)
    else:
        now = dt.datetime.strftime(dt.datetime.now(),"%Y-%m-%d %H:%M:%S")
        ds.attrs["history"] = "File modified on {} by University of Bristol ACRG group. Modification details listed within global attributes: {}, ".format(now,mod_attr)
    
    return ds


def ds_check_internal_unique(ds,axis="time"):
    '''
    The ds_check_internal_unique function checks whether an xarray.Dataset has repeat values on the "time" axis.
    Any time values which are repeated are modified by a small time increment to make the values unique. 
    
    Args:
        ds (xarray.Dataset):
            Dataset with time axis.
        axis (str, optional):
            Name of time axis.
            Default = "time"
    
    Returns:
        xarray.Dataset:
            If no repeats are present:
                Original dataset is returned
            If repeats are present:
                The repeats of any time values are modified by applying a small random increment so the time values
                are not identical.
    '''
    
    if len(ds[axis].values) > len(np.unique(ds[axis].values)):
        axis_copy = ds[axis].values
        unique,inverse,counts = np.unique(axis_copy,return_inverse=True,return_counts=True)
        repeat_index_0 = np.where(counts > 1)[0] # Find first indices of all elements with repeats
        repeat_all_index = [np.where(inverse == i)[0] for i in repeat_index_0] # Extract all indices for these repeats
        for indices in repeat_all_index:
            for i in indices[1:]:
                # Add very small random value to repeat of a time value to avoid two times being exactly the same
                axis_copy[i] += np.timedelta64(random.randrange(-1000,1000,1),'us')
                #axis_copy[i] += np.timedelta64(1,'s')
        #ds[axis].values = axis_copy
        ds = ds.assign_coords(**{axis:axis_copy})
        
        # Add attribute describing modification made to original data
        mod_attr = "repeat_time_modified"
        if mod_attr not in list(ds.attrs.keys()):
            ds = add_history_attr(ds,mod_attr)
            ds.attrs[mod_attr] = "Repeated times were found at original indices {} (may not be the same if data has been binned). Small random increments were added to any repeat values to allow xarry to distinguish them.".format(repeat_all_index)
        else:
            ds.attrs[mod_attr] += " Also found at indices {}.".format(repeat_all_index)
    
    return ds

def find_network(site):
    '''
    The find_network function identifies the associated networks for a given site from the
    "acrg_site_info_2018.json" file.
    Args:
        site (str) :
            Site identifier. See "acrg_site_info_2018.json" for full list of options.
    Returns:
        list :
            List of networks for each site.
    '''
    try:
        networks = list(acrg_obs.read.site_info[site].keys())
    except KeyError:
        networks = ["unknown"]
    if len(networks) > 1:
        print("Multiple networks detected for site {}: {}".format(site,networks))

    return networks

def extract_dates(ds,dim="time",dtype='M8[D]'):
    '''
    The extract_dates function converts datetime objects into date strings for all values along a given dimension.
    Note: dimension should contain an array of np.datetime64 objects.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset which contains dimension specified by dim
        dim (str, optional): 
            String of dimension containing datetime objects. Default = "time"
        dtype (str, optional):
            Specific dtype value to use to cast datetime object. Default = 'M8[D]'
        
    Returns:
        np.array (str) : dates as an array of strings
    '''
    #return ds[dim].values.astype(dtype).astype(str) # Extract dates from time column (cast as 8 byte datetime format (M8) in whole days [D] then cast as a string)
    return np.datetime_as_string(ds[dim].values.astype(dtype)) # Extract dates from time column (cast as 8 byte datetime format (M8) in whole days [D] then cast as a string)

def units(name):
    '''
    The units function defines units for output variables.
    
    Data variable units defined within this function are:
        ch4_profile_apriori: ppb
        xch4: ppb
        xch4_uncertainty: ppb
        lat: degrees
        lon: degrees
        pressure_levels: hPa
        time: seconds since 1970-01-01 00:00:00
        pressure_weights: unitless
        xch4_averaging_kernel: unitless
    
    Args:
        name (str) : 
            Name of output variable
    
    Returns:
        str : 
            Unit of name (if present), None otherwise
    '''
    
    unit_data = {}
    
    unit_data["ch4_profile_apriori"] = "1e-9"
    unit_data["xch4"] = "1e-9"
    unit_data["xch4_uncertainty"] = "1e-9"
    
    unit_data["lat"] = "degrees"
    unit_data["lon"] = "degrees"
    
    unit_data["pressure_levels"] = "hPa"
    
    unit_data["time"] = "seconds since 1970-01-01 00:00:00"
    
    unit_data["pressure_weights"] = "unitless"
    unit_data["xch4_averaging_kernel"] = "unitless"
    
    try:
        unit = unit_data[name]
    except KeyError:
        return None
    else:
        return unit
    
def split_output(ds,index,mapping=None,data_vars=[],split_dim="time",ident=None,ident_sep=','):
    '''
    The gosat_split_output function creates a new dataset based on input index values along the 
    split dimension (e.g. time).
    This is used to separate data points into subsets (usually to write out to file).
    
    If an ident column is specified the values within this column are also split out into a new "id" dimension.
    
    Note: units for each variable are set by the units() function (and so to include a unit must be
    defined there).
    
    Args:
        ds (xarray.Dataset) : 
            Dataset with consistent dimension to split along.
        index (int/list) : 
            Index value or values to use for the split.
        mapping (dict/OrderedDict, optional) : 
            Mapping between the data variable names within ds and the output data variable names. 
            If no mapping is specified this will use the same data variable names in input dataset 
            for the output dataset.
        data_vars (list, optional) : 
            Data variables within ds to include in output dataset. (list)
            If both mapping and data_vars are specified, mapping will supercede data_vars.
            If no mapping is included and data_vars if left blank ([] by default) all data variables will 
            be copied.
        split_dim (str, optional) : 
            Dimension to split along (i.e. that index is relevant to). Default = 'time' (str)
        ident (str, optional) : 
            Data variable which includes identifier values for all data points within a bin e.g. "exposure_id".
            If ident is set then this data variable will be split out to contain an extra dimension ("id")
            The ident data variable is expected to contain some identifier for each data point separated 
            by the ident_sep value (e.g. commas)
            e.g. ['2010123106380100271013,2010123106380100271014,2010123106380100271015']
        ident_sep (str, optional) : 
            Separator value to use for identifier column.
            Default is a comma (','). Only used if ident is set.
    
    Returns:
        xarray.Dataset : 
            Dataset containing values specified by indices along split_dim (e.g. time)
            Identifier column will be split into an extra dimension if ident is specified
    '''
    
    if not data_vars and not mapping:
        data_vars = ds.data_vars()
    
    if not mapping:
        mapping = OrderedDict([(name,name) for name in data_vars]) # If no mapping specified, assume data variable names in created ds match input names
   
    if data_vars:
        if ident not in data_vars:
            print('WARNING: Identifier column {0} is not within input data_vars: {1}. No ident column will be included'.format(ident,data_vars))
    
    if isinstance(index,(int,np.integer)):
        indices = [index]
    else:
        indices = index
    
    coords = OrderedDict([])
    data = OrderedDict([])
    
    # Create data and coords to be included within dataset
    for name,new_name in list(mapping.items()):
        
        dims = ds[name].dims # Extract dimensions for variable from current dataset
        data_var = ds[name][indices]
        
        # Format the identifier data variable (e.g. exposure_id) to contain an extra "id" dimension to allow for multiple values
        if name == ident: # what is the content of exposure_id
            identifiers = [value.split(ident_sep) for value in data_var.values] # Split identifier value by the ident_sep value (e.g. ',')
            identifiers = np.array(list(itertools.zip_longest(*identifiers,fillvalue=np.nan))).T # Create array with consistent dimensions for "id" and fill in any gaps with np.nan values
            #identifiers = np.array(list(itertools.izip_longest(*identifiers,fillvalue=np.nan))).T # Create array with consistent dimensions for "id" and fill in any gaps with np.nan values
           
            id_dim_name = "id" # Define new dimension name
            split_dim_dim,id_dim_dim = identifiers.shape # Define dimensionality of new dimension and dimension we're splitting on
            id_dim_coord = np.arange(1,id_dim_dim+1)
            split_dim_coord = ds[split_dim][indices] # Extract associated split dimension coordinate e.g. time coord
            
            id_coords = OrderedDict([(split_dim,split_dim_coord),
                                     (id_dim_name,id_dim_coord)])
            id_dims = OrderedDict([(split_dim,split_dim_dim),
                                   (id_dim_name,id_dim_dim)])
            
            data_var = xray.DataArray(identifiers,coords=id_coords,dims=id_dims) # Reset data_var as new variable with extra dimension
            
            # Add suitable attributes for the new ident data variable
            if ident == "exposure_id":
                data_var.attrs["short_description"] = "Exposure identification number of the sounding for each data point within the bin."
                data_var.attrs["long_name"] = "exposure_id_mult"
            else:
                data_var.attrs["short_description"] = "{} label for each data point within the bin.".format(name)
                data_var.attrs["long_name"] = ident
            
            data_var[id_dim_name].attrs["short_description"] = "Number for each data point which has been combined within the bin."
            data_var[id_dim_name].attrs["long_name"] = "id"
                  
        # Define units for each data variable, if not already specified
        unit = units(new_name)
        if unit and ("units" not in data_var.attrs):
            data_var.attrs["units"] = unit
        data[new_name] = data_var
        
        # Extract relevant coords from input dataset
        for dim in dims:
            if dim not in list(coords.keys()): 
                if dim == split_dim:
                    coords[dim] = ds[dim][indices] # If dim is split_dim dimension (e.g. time), extract only the relevant values
                else:
                    coords[dim] = ds[dim] # Extract all coords for other dimensions
        
        # Add extra "id" dimension if identifier column is included.
        if name == ident:
            coords[id_dim_name] = id_coords[id_dim_name]
    
    ds_output = xray.Dataset(data,coords=coords) # Create new dataset
    
    # Add units for coordinates
    for dim in ds_output.coords:
        unit = units(dim)
        if unit and ("units" not in ds_output[dim].attrs):
            ds_output[dim].encoding["units"] = unit # For e.g. "time" dimension, "units" is a special label set via an initial encoding and so must be written this way (rather than as an attribute) or this will produce an error.
    ds_output = ds_output.assign_attrs(ds.attrs) # Copy across global attributes from input dataset
    
    
    return ds_output

def write_netcdf(ds,filename,overwrite=False):
    '''
    The write_netcdf function writes an xarray.Dataset to a netCDF file.
    The function will also checks whether the file exists and only overwrite if this option is set to True.
    
    Args:
        ds (xarray.Dataset) :
            Dataset to write to file
        filename (str) :
            Filename (including path information) for the output file
        overwrite (bool) :
            Whether to overwrite an existing file if present.
            Default = False
    
    Returns:
        None
        
        Writes dataset to file (filename)
    '''    
    if not overwrite:
        if os.path.isfile(filename):
            print('\nERROR: {0} already exists. Data has not been written to file because overwrite=False.\n'.format(filename))
        else:
            print('Writing to filename:',filename)    
            ds.to_netcdf(filename)
    else:
        print('Writing to filename:',filename)    
        ds.to_netcdf(filename,mode="w")

def output_filename(output_directory,network,instrument,date,species,inlet=None,num=None):
    '''
    The gosat_output_filename function creates an output filename of the correct format for gosat based on the 
    inputs.
    
    Filenames are of the form: "instrument"_gosat_"date (reformated)"-"num"_"species"-"inlet".nc
    e.g. /shared_data/air/shared/obs/GOSAT/GOSAT-INDIA/gosat-fts_gosat_20120920-09_ch4-column.nc
    
    Args:
        output_directory (str) : 
            Top level for output directory e.g. "/shared_data/air/shared/obs/"
        network (str) : 
            Which network is being considered e.g. "GOSAT/GOSAT-INDIA"
            This will be used to create the output path e.g. "/shared_data/air/shared/obs/GOSAT/GOSAT-INDIA/"
        instrument (str) : 
            Instrument on satellite being used. e.g. "gosat-fts"
        date (str) : 
            Date the measurements are relevant to e.g. "2010-01-01"
        species (str) : 
            Species being considered e.g. "ch4". Should be defined within "openghg/supplementary data species_info.json" file
        inlet (str/None, optional) : 
            Additional information to add to ch4 measurement e.g. column
        num (str/None, optional) : 
            Number to append to the filename if gosat data is being written out as multiple files over one day.
            
     Returns:
         str: 
             GOSAT observation filename with path information    
    '''
    #Example /shared_data/air/shared/obs/GOSAT/GOSAT-INDIA/gosat-fts_gosat_20120920-09_ch4-column.nc
    
    output_directory = os.path.join(output_directory,network)
    
    if 'gosat' in instrument.lower():
        satellite = 'gosat'
    elif 'tccon' in instrument.lower():
        satellite = 'tccon'
    else:
        satellite = 'oco2'
    
    date = date.replace('-','') # Turn date from e.g. 2012-09-20 to 20120920
    
    if num:
        date = '-'.join([date,num]) # Create composite string of datetime and num (if present) e.g. "20120920-01"
    if inlet:
        species = '-'.join([species,inlet]) # Create composite string of species and inlet (if present) e.g. "ch4-column"
  
    filename = '_'.join([instrument,satellite,date,species]) # Create filename string joined by "_" e.g. "gosat-fts_gosat_20120920-09_ch4-column"
    filename += '.nc' # Add file extension
    
    filename = os.path.join(output_directory,filename)
       
#   From cf_data_process.output_filename   
#        return join(output_directory,
#                network + "/" + \
#                network + "-" + \
#                instrument + "_" + \
#                site + "_" + \
#                year + "0101_" + \
#                species + "-" + \
#                inlet + ".nc")
#    
    return filename
  
def output(ds,site,species,network=None,
           file_per_day=False,output_directory=obs_directory,
           overwrite=False):
    '''
    The gosat_output function creates output netCDF files each containing data related to one data point per file.
    The data variables within each netCDF file are:
        "ch4_profile_apriori","lat","lon","pressure_levels","pressure_weights","xch4","xch4_averaging_kernel",
        "xch4_uncertainty"
    Units for these values are set by the units() function.
    
    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file as an xarray Dataset.
            Should have had co-ordinates and dimensions assigned with gosat_add_coords() function.
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA (should be defined within site_info.json)
        species (str) : 
            Species of interest e.g. "ch4" (should be defined within species_info.json).
        file_per_day (bool, optional) : 
            Output all results to one file per day rather than splitting out per time point. Default = False.
        output_directory (str, optional) : 
            Top level directory to write output. Full path will be based on the "network" related to "site" .
            Default = obs_directory (defined at top of file).
        overwrite (bool, optional) : 
            Whether to overwrite any files already present in the output_directory + network folder.
            Default = False.
    
    Returns:
        None
        
        Writes output to multiple .nc files (split on time axis).
    
    '''
    if network == "oco2":
        # Define data variables to be written to output
        out_data_vars = [
            f"x{species}",  # Main species variable
            f"x{species}_uncertainty",  # Uncertainty for the species
            "lat",
            "lon",
            "pressure_levels",
            "pressure_weights",
            f"x{species}_averaging_kernel",  # Averaging kernel for the species
            f"{species.split('_')[0]}_profile_apriori",  # Apriori profile for the species
            # "mode"
        ]
    else:
        out_data_vars = [
                f"x{species}",  # Main species variable
                f"x{species}_uncertainty",  # Uncertainty for the species
                "lat",
                "lon",
                "pressure_levels",
                "pressure_weights",
                f"x{species}_averaging_kernel",  # Averaging kernel for the species
                f"{species.split('_')[0]}_profile_apriori",  # Apriori profile for the species
                "mode"
            ]
    
    # Map to input dataset (from GOSAT data)
    data_vars = ["latitude" if item=="lat" else item for item in out_data_vars]
    data_vars = ["longitude" if item=="lon" else item for item in data_vars]
    data_vars = ["pressure_weight" if item=="pressure_weights" else item for item in data_vars]

    if network.lower() != "oco2":
        data_vars = ["retr_flag" if item=="mode" else item for item in data_vars]
        
    data_var_mapping = OrderedDict([(name,new_name) for name,new_name in zip(data_vars,out_data_vars)])
    
    #if site == None:
    #    site="global"
    
    # Set dimension on which data will be split into separate files
    split_dim="time"
    
    # Set name of data variable which includes the data point identifiers
    ident = "exposure_id"    
    if 'gosat' in network.lower():
        network = find_network(site)[0] # Using first site as default.
        instrument = 'gosat-fts'
    elif 'oco2' in network.lower():
        instrument = 'oco2-spectrometer'
    else:  
        instrument = f'tccon-{site}'
        
    inlet = 'column'
    species = species.lower()
    
    #all_dates = ds["time"].values.astype('M8[D]').astype(str) # Extract dates from time column (cast as 8 byte datetime format (M8) in whole days [D] then cast as a string)
    all_dates =  extract_dates(ds,dim=split_dim)
    dates = np.unique(all_dates)
    
    full_output_directory = os.path.join(output_directory,network)
    if not os.path.isdir(full_output_directory):
        os.makedirs(full_output_directory)
    
    for date in dates:
        wh_date = np.where(all_dates == date)[0] # Find indices for each date
        if file_per_day:
            ds_output = split_output(ds,index=wh_date,mapping=data_var_mapping,split_dim=split_dim,ident=ident)
            # Create filename and write dataset to file
            filename = output_filename(output_directory,network,instrument,date,species,inlet=inlet)
            ds_output.attrs["id"] = os.path.split(filename)[1]
            write_netcdf(ds_output,filename,overwrite=overwrite)
        else:
            for ID,index in enumerate(wh_date):
                ds_output = split_output(ds,index=index,mapping=data_var_mapping,split_dim=split_dim,ident=ident)

                # Create filename and write dataset to file
                ID_str = str(ID+1).zfill(3) # Number to add to filename - three digit with leading zeros
                filename = output_filename(output_directory,network,instrument,date,species,num=ID_str,inlet=inlet)
                ds_output.attrs["id"] = os.path.split(filename)[1]
                write_netcdf(ds_output,filename,overwrite=overwrite)
 
def output_name(ds,site,network,max_level=17,use_name_pressure=False,pressure_domain=None,
                      pressure_base_dir=name_pressure_directory,
                      pressure_max_days=31,pressure_day_template=True,name_directory=name_csv_directory,
                      file_per_day=False,max_points=60,overwrite=False):
    '''
    The gosat_output_name function creates a csv file in the correct format for input into NAME-III.
    Note: it is assumed input pressure values are in hPa values and so are converted to Pa for NAME (*100.)
    
    Args:
        ds (xarray.Dataset) : 
            GOSAT CH4 Level 2 Data Product file as an xarray Dataset.
        site (str) : 
            Specified sub-set defined for gosat e.g. GOSAT-INDIA. Should be defined within site_info.json.
            Note: at the moment this is just used in the filename and not in the folder structure
        max_level (int, optional) : 
            Maximum level to include from input GOSAT data (up to 20).
            At the moment NAME footprints go up to 19km which corresponds to level 17 in GOSAT.
        use_name_pressure (bool, optional) : 
            Whether to use the NAME surface pressure rather than the GOSAT surface pressure value for each data 
            point.
        pressure_domain (str, optional) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            Must be included if use_name_pressure=True.
            Check $DATA_PATH/LPDM/surface_pressure folder to see which domains currently exist.
        pressure_base_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        pressure_max_days (int, optional) : 
            Number of days tolerance to allow when using time stamp to find relevant NAME pressure values. 
            Default = 31 (days).
        pressure_day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
        name_directory (str, optional) : 
            Top level directory to write files for NAME. Full path will be based on the "network" related to "site" 
            Default defined at the top of this module.
        file_per_day (bool, optional) : 
            Output all results to one file per day rather than splitting out per time point. 
            Default=False.
        max_points (int/None, optional) :
            Only applicable if file_per_day is True.
            Maximum number of points to write if writing file per day. If number of points within the day
            exceeds max_points multiple files will be written out with a character after the date of "-A","-B",...
        overwrite (bool, optional) : 
            Allow any files already present within the full path to name_directory to be overwritten. 
            Default = False
    
    Returns:
        None
        
        Writes output to multiple .csv files (split on time axis).
    '''
    # Note column mapping here is out:in values rather than in:out because some columns in output do not map to input data variables.
    col_mapping = OrderedDict([('ID_Level','lev'),
                               ('Time','time'),
                               ('x','longitude'),
                               ('dx','dlon'),
                               ('y','latitude'),
                               ('dy','dlat'),
                               ('z','pressure_levels'),
                               ('dz',None)])
    if network.upper()=='GOSAT':           
        try:
            ds[col_mapping['dx']]
            ds[col_mapping['dy']]
        except KeyError:
            sounding = 10.5
            dlat = np.zeros(len(ds[col_mapping['y']]))+distance_lat(sounding)
            dlon = np.array([distance_lon(sounding,lat) for lat in ds[col_mapping['y']]])
            axis1 = col_mapping['Time']
            dlat_da = xray.DataArray(dlat,coords={axis1:ds[axis1]},dims={axis1:len(ds[axis1])})
            dlon_da = xray.DataArray(dlon,coords={axis1:ds[axis1]},dims={axis1:len(ds[axis1])})
            ds = ds.assign(**{col_mapping['dx']:dlon_da,col_mapping['dy']:dlat_da})
        
        if site is None:
            site = 'global'
        network = find_network(site)[0] # Using first site as default        
    
    columns=list(col_mapping.keys())
    
    split_dim = "time"
    out_index = "ID_Level"
    
    name_pressure_convert = 100. # Input from GOSAT is in hPa and need to convert to Pa for NAME
    
    all_dates = extract_dates(ds,dim=split_dim)
    dates = np.unique(all_dates)
    
    pressure_levels,dpressure = define_pressure_levels(ds,pressure_domain=pressure_domain,
                                                       use_name_pressure=use_name_pressure,
                                                       pressure_base_dir=pressure_base_dir,
                                                       max_days=pressure_max_days,
                                                       day_template=pressure_day_template)
        
    name_directory = os.path.join(name_directory,network)
    if not os.path.isdir(name_directory):
        os.makedirs(name_directory)
    
    for date in dates:
        wh_date = np.where(all_dates == date)[0] # Find indices for each different date
        for ID,index in enumerate(wh_date):
            data = OrderedDict()
            ID_str = '{num:03d}'.format(num=ID+1) # Number for each point in a day
            
            for col in columns:
                name = col_mapping[col]
                if name:
                    ds_column = ds[name]
                    if name == "pressure_levels":
                        data[col] = pressure_levels[index][:max_level]*name_pressure_convert
                    elif name == "lev":
                        
                        if file_per_day:
                            # Write out ID as a str combination of the point number and level each with leading zeros
                            data[col] = np.array([ID_str+value.zfill(2) for value in (ds_column[:max_level].values+1).astype("str")])
                        else:
                            # Write out ID (levels) as two digit numbers with a leading zero e.g. 01, 02, ... 21, 22
                            data[col] = np.array([value.zfill(2) for value in (ds_column[:max_level].values+1).astype("str")])
                    else:
                        data[col] = np.array([ds_column[index].values]*max_level) # Populate each level with the same values e.g. lat,lon,dlat,dlon
                else:
                    if col == "dz":
                        data[col] = dpressure[index][:max_level]*name_pressure_convert

            df_output = pd.DataFrame(data,columns=columns)
            df_output.set_index(out_index,inplace=True)
            
            if file_per_day:
                filename = site+"_"+date.replace('-','')
                if max_points:
                    if len(wh_date) > max_points:
                        letter_split = chr(ord("A") + ID//max_points)
                        filename = "{}-{}".format(filename,letter_split)
                    ID_1 = ID%max_points
                else:
                    ID_1 = ID
                filename += '.csv'
                filename = os.path.join(name_directory,filename)
                
                if not overwrite:
                    if ID_1 == 0 and os.path.isfile(filename):
                        print('\nERROR: {0} already exists. Data has not been written to file because overwrite=False.\n'.format(filename))
                        break
                    elif not os.path.isfile(filename):
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f')
                    else:
                        #print('Appending to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f',mode='a',header=False)    
                else:
                    if (not os.path.isfile(filename)) or (os.path.isfile(filename) and ID_1 == 0):
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f')
                    else:
                        #print('Appending to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f',mode='a',header=False)
            else:
                filename = "{site}_{date}-{ID}.csv".format(site=site,date=date.replace('-',''),ID=ID_str)
                filename = os.path.join(name_directory,filename)
                if not overwrite:
                    if os.path.isfile(filename):
                        print('\nERROR: {0} already exists. Data has not been written to file because overwrite=False.\n'.format(filename))
                        break
                    else:
                        print('Writing to filename: {}'.format(filename))
                        df_output.to_csv(filename,float_format='%f')
                else:
                        print('Writing to filename: {}'.format(filename))    
                        df_output.to_csv(filename,float_format='%f')


radius = 6371 #km
#radius = 6367.5 #km
def distance_lat(distance,radius=radius):
    '''
    Calculate equivalent difference in latitude for a distance in km.
    
    Args:
        distance (int/float) : 
            Distance in km.
        radius (int/float, optional)   : 
            Radius of the Earth (default to 6371km)
        
    Returns:
        float : 
            Latitude difference in degrees
    '''
    # Haversine distance formula is given as distance = 2*asin(sqrt(a))*radius
    a = math.sin(distance/(2.*radius))*math.sin(distance/(2.*radius))
    
    # Full formula for a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    # When dlon=0, lon1=lon2, a reduces to a = sin(dlat/2)**2
    dlat = 2*math.asin(math.sqrt(a))
    dlat = np.abs(dlat)
    
    dlat = math.degrees(dlat)
    return dlat

def distance_lon(distance,lat,radius=radius):
    '''
    Calculate equivalent difference in longitude for a distance in km at a given latitude.
    
    Args:
        distance (int/float) : 
            Distance in km.
        lat (int/float) : 
            Latitude in degrees (between but not including -90 and +90. Note: max value = 89.995)
        radius (int/float) : 
            Radius of the Earth (default to 6371km)
        
    Returns:
        float : 
            Longitude difference in degrees
    '''
    lat = math.radians(lat)
    # Haversine distance formula is given as distance = 2*asin(sqrt(a))*radius
    a = math.sin(distance/(2.*radius))*math.sin(distance/(2.*radius))
    
    # Full formula for a = sin(dlat/2)**2 + cos(lat1)*cos(lat2)*sin(dLon/2)**2
    # When dlat=0, lat1=lat2, a reduces to a = cos^2(lat)*sin^2(dlat/2)
    b = math.sqrt(a)/math.cos(lat)
    dlon = 2.*math.asin(b)
    dlon = np.abs(dlon)
    
    dlon = math.degrees(dlon)
    return dlon

      
def define_pressure_levels(ds,pressure_domain=None,p_column="pressure_levels",set_as_edges=False,
                           include_bounds=False,use_name_pressure=False,pressure_NAME=None,
                           columns=["latitude","longitude","time"],pressure_base_dir=name_pressure_directory,
                           max_days=31,day_template=True):
    '''
    The define_pressure_levels function scales the pressure_levels and dpressure values for input into NAME.
    Note: values to be output to obs/ files *should NOT* be run through this function.
    
    The pressure_levels provided with the GOSAT data (in particular) define instantaneous pressure points and not the boundaries or 
    midpoints of each of the layers.
    The first pressure level (pressure_0) is treated as the surface pressure value.
 
    To create layers from these pressure levels the boundaries are defined as the midpoint between each level. If the levels are 
    unevenly spaced the pressure levels themselves are then re-defined to be at the midpoint within each layer 
    (or as the edge if set_as_edges=True). The top and bottom layer are treated as follows:
         - Bottom layer - Either whole or 1/2 difference between NAME or satellite surface pressure is used to define layer extent and
                          boundaries depending on if set_as_edges is True or False.
         - Top layer - set layer extent and as matching previous level value OR set layer extent and boundaries based on a minumum 
                       pressure (set as 0 hPa), whichever dpressure is smaller.
     
    Args:
        ds (xarray.Dataset) : 
            Dataset containg multiple data points with a set of pressure levels for each.
        pressure_domain (str/None, optional) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            Must be specified if pressure_NAME is not included explicitly.
            Check $DATA_PATH/LPDM/surface_pressure folder to see which domains currently exist.
        p_column (str, optional) : 
            Name for the pressure_levels data (str). Default = "pressure_levels"
        set_as_edges (bool, optional) : 
            Whether to treat the pressure levels as points at the edge of the layers rather than within the layers .
            Default = False.
            Note: level 0 will ALWAYS be treated as being the surface pressure (i.e. at the boundary) regardless 
            of this parameter.
        include_bounds (bool, optional)   : 
            Whether to output the boundary values for the layers as well as the pressure_levels and dpressure.
            Default = False.
        use_name_pressure (bool, optional) : 
            If to use the NAME surface pressure rather than the GOSAT value included within the data. 
            Default = False.
        pressure_NAME (numpy.array, optional) : 
            This and the following parameters are only used when use_name_pressure=True.
            If pressure from NAME run has already been extracted, this can be specified with this parameter to 
            save computing time.
            If not specified, columns from ds and pressure_base_dir and pressure_domain will be used to extract matching 
            pressure values.
        columns (list, optional) : 
            Names of data variables or co-ords within input Dataset for the latitude, longitude and time values 
            (3 item list). Default = ["latitude","longitude","time"]
        pressure_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        max_days (int, optional) : 
            Maximum number of days from time within ds to search for the relevant pressure data. 
            Default = 31 (days).
        day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Default = True.
                            
    
    Returns:
        if include_bounds == True:
            (np.array,np.array,np.array):
                pressure_levels,dpressure,pressure_bounds
        else:
            (np.array,np.array): 
                pressure_levels,dpressure

    '''

    pressure_levels = ds[p_column].copy(deep=True).values
    delta_pressure = pressure_levels[:,:-1] - pressure_levels[:,1:]
    
    if use_name_pressure:
        if pressure_NAME is None:
            if pressure_domain is not None:
                pressure_NAME = name_pressure_match(ds,pressure_domain=pressure_domain,
                                                pressure_base_dir=pressure_base_dir,columns=columns,
                                                day_template=day_template,max_days=max_days)
            else:
                raise Exception("Pressure_domain must be specified if pressure values need to be \
                                extracted to use define_pressure_levels function.")
    else:
        pressure_NAME = None
    
    # Use midpoint between points to define the boundary
    # Scale the pressure points to be at the centre of that boundary
    pressure_levels,dpressure = midpoint_bounds(pressure_levels,delta_pressure,pressure_NAME,set_as_edges=set_as_edges)

    if include_bounds:
        pressure_bound = np.zeros((len(pressure_levels),len(pressure_levels[0])+1))
        pressure_bound[:,:-1] = pressure_levels + np.abs(dpressure/2.)
        pressure_bound[:,-1] = pressure_levels[:,-1] - np.abs(dpressure[:,-1]/2.)
        
        return pressure_levels,dpressure,pressure_bound
    else:
        return pressure_levels,dpressure 
    
                 
def name_pressure_match(ds,pressure_domain,columns=["latitude","longitude","time"],p_column="surface_pressure",
                        pressure_base_dir=name_pressure_directory,pressure_convert=1/100.,
                        max_days=31,day_template=True):
    '''
    The name_pressure_match function finds the corresponding NAME surface pressure values based on a set of
    lat,lon and time values extracted from the input dataset.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset containing pressure values to match the NAME surface pressure values to. Within the 
            dataset the pressure values should be related to latitude, longitude and time values.
        pressure_domain (str) :
            Domain over which surface pressure values have been extracted (can be distinct from 
            domain if pressure_domain contains area of domain).
            Check $DATA_PATH/LPDM/surface_pressure folder to see which domains currently exist.
        columns (list, optional) : 
            Names of data variables or co-ords within input Dataset for the latitude, longitude and time values.
            Default = ["latitude","longitude","time"]
        pressure_base_dir (str, optional) : 
            Base directory containing the NAME output files for the SurfacePressure run.
            Filename is assumed to be of the form "Pressure_C1_*.txt"
            See name_pressure_file() function for more details.
        pressure_convert (float, optional) : 
            If pressure values extracted from NAME are not in the required units, pressure_convert
            should be set to the scaling factor to convert these units.
            By default, we assume we want to convert to hPa from NAME input in Pa (pressure_covert=1/100.)
        max_days (int, optional) : 
            Number of days tolerance to allow when using time stamp to find relevant NAME pressure values. 
            Default = 31 (days).
        day_template (bool, optional) :
            Use nearest day as a template for the change of pressure over the course of the day and match
            to the nearest time on that day.
            E.g. if datetime is 2012-05-01 03:00:00, max_days is 31 and nearest day is 2012-01-01 then 
            use entry from 2012-01-01 03:00:00 (rather than 2012-02-01 00:00:00, which would be the 
            nearest entry).
            Note: Assumes NAME run for e.g. 2012-01-01 ends with datetime 2012-02-01 00:00:00.
            Default = True.
        
    Returns:
        numpy.array : 
            Array of pressure values. One for each lat,lon,time value.    
    '''
    
    #if len(lat) != len(lon) or len(lat) != len(time):
    #    print 'To extract matched name pressure values the same number of latitude, longitude and time values must be provided.'
    #    print 'Lat: {0}, Lon: {1}, Time: {2}'.format(len(lat),len(lon),len(time))
    #    return None
    
    lat = ds[columns[0]].values
    lon = ds[columns[1]].values
    time = ds[columns[2]].values
    
    # Can't apply a list tolerance with current version of xarray. May be able to add with newer versions.
    #lat_tolerance = 5.0
    #lon_tolerance = 5.0
    #time_tolerance = np.timedelta64(max_days,'D')
    #tolerance = [lat_tolerance,lon_tolerance,time_tolerance]
    
    start_date = min(time).astype('M8[D]').astype(str)
    end_date = max(time).astype('M8[D]').astype(str)
    
    pressure_dir = os.path.join(pressure_base_dir,pressure_domain)
    
    print('Extracting pressure_NAME from name_pressure function')
    pressure_NAME = name_pressure(pressure_dir,start_date=start_date,end_date=end_date,name=p_column,column_names=columns,max_days=max_days)
    #matched_pressure_NAME = pressure_NAME.sel(method="nearest",**{columns[0]:lat,columns[1]:lon,columns[2]:time}) # Returns grid of values at timexlatxlon
    #matched_pressure_NAME = matched_pressure_NAME[p_column].values
    #matched_pressure_NAME = np.diagonal(np.diagonal(matched_pressure_NAME)) # Only want values in this grid on the diagonal for our purposes (i.e. time1,lat1,lon1; time2,lat2,lon2, not time1,lat2,lon1 etc.).
    matched_pressure_NAME = np.zeros(len(lat))
    if day_template:
        end_datetime = np.max(pressure_NAME[columns[2]].values) # Extract maximum datetime value
    for i,la,lo,t in zip(list(range(len(lat))),lat,lon,time):
        if not day_template or t <= end_datetime:
            match = pressure_NAME.sel(method="nearest",**{columns[0]:la,columns[1]:lo,columns[2]:t})
        else:
            # Assumes offset is +1 day because we expect NAME run for e.g. 2012-01-01 to end with datetime 2012-02-01 00:00:00.
            max_offset_days = np.timedelta64((t - end_datetime),'D') + np.timedelta64(1,"D")
            pressure_NAME_offset = pressure_NAME.copy(deep=True)
            t_offset = pressure_NAME_offset[columns[2]]+max_offset_days 
            pressure_NAME_offset = pressure_NAME_offset.assign_coords(**{columns[2]:t_offset}) # Create copy of dataset with date offset to date we're looking for.
            #day_tolerance = tolerance[:1] + [np.timedelta64(1,'D')]
            match = pressure_NAME_offset.sel(method="nearest",**{columns[0]:la,columns[1]:lo,columns[2]:t})
        matched_pressure_NAME[i] = match[p_column].values
    
    matched_pressure_NAME = matched_pressure_NAME*pressure_convert
    
    return matched_pressure_NAME
        

def name_pressure(directory,start_date=None,end_date=None,name="surface_pressure",
                  column_names=["latitude","longitude","time"],
                  set_columns=["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"],max_days=31):
    '''
    The name_pressure function extracts the pressure values from a special NAME run to find surface pressure 
    values across a given domain.
    These values will often be spread over multiple files and this function can take several files from a 
    directory and create one dataset.
    
    Args:
        directory (str) : 
            Directory containing the NAME output files for the SurfacePressure run.
        start_date (str, optional) : 
            Start date range of interest for the pressure values.
        end_date (str, optional) : 
            End date range of interest for the pressure values (optional)
        name (str, optional) : 
            Name to use when creating new pressure data variables within output Dataset.
            Default = "surface_pressure"
        column_names (list, optional) : 
            Names for co-ordinate axes within dataset. Should be a 3-item list for latitude, longitude and time.
            Default = ["latitude","longitude","time"]
        set_columns (list, optional) : 
            Column names within input file which contain longitude, latitude and height (in that order).
            Function assumes these columns will be the first in the file and all other populated columns
            contain the surface pressure (or topography) values.
            Default = ["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"] based on current BackRuns script.
        max_days (int, optional) : 
            Number of days tolerance to allow when using time stamp to find relevant NAME pressure values. 
            Default = 31 (days).
    
    Returns:
        xarray.Dataset : 
            Surface pressure data covering at least start_date - end_date range (if specified)    
    '''
    
    #search_str = "Pressure_C1_*_{0}d.txt".format(max_days)
    search_str = "Pressure_C1_*.txt*"
    
    if start_date and end_date:
        print('Setting tolerance of {0} days back when searching for NAME pressure files.'.format(max_days))
        start_date = str(np.datetime64(start_date) - np.timedelta64(max_days-1,'D'))#.astype(str)
        end_date = str(np.datetime64(end_date) + np.timedelta64(1,'D'))#.astype(str)
    
    files = extract_files(directory,search_str=search_str,start=start_date,end=end_date)
    
    if len(files) == 0:
        raise Exception('No NAME pressure files found within directory: {0} for date range {1} - {2} (using search_str {3})'.format(directory,start_date,end_date,search_str))
        #return None
    
    for i,filename in enumerate(files):
        if i == 0:
            name_ds = name_pressure_file(filename,name=name,column_names=column_names,set_columns=set_columns)
        else:
            ds = name_pressure_file(filename,name=name,column_names=column_names,set_columns=set_columns)
            ds = ds_check_unique(ds,name_ds,dim_apply="time")
            name_ds = name_ds.merge(ds)
    
    return name_ds   

def ds_check_unique(ds1,ds2,dim_apply="time"):
    '''
    The ds_check_unique function checks for any repeats between two datasets on the dim_apply axis.
    Any repeats are removed from the first dataset (ds1) and returned.
    
    Args:
        ds1 (xarray.Dataset) : 
            First Dataset with dim_apply dimension
        ds2 (xarray.Dataset) : 
            Second Dataset with dim_apply dimension to compare to.
        dim_apply (str, optional) : 
            Dimension to check for repeats (str). Default = "time"
    
    Returns:
        xarray.Dataset : 
            ds1 with any matching values removed.    
    '''
    
    t1 = ds1[dim_apply]
    t2 = ds2[dim_apply]
    
    filt = ~np.in1d(t1,t2) # Compare values in ds2 to ds1 and check for repeats. Want to keep values that aren't repeats.
    ds = apply_filter(ds1,filter_array=filt,dim_apply=dim_apply)
    
    return ds

def apply_filter(ds,filter_array,dim_apply='time'):
    '''
    The apply_filter function applies a filter array of indices to a dataset.
    All data variables with the dim_apply value as their first dimension will be filtered.
    
    Args:
        ds (xarray.Dataset) : 
            Dataset to filter. 
        filter_array (np.array) : 
            Array of indices to keep within the dataset. All other indices will be removed.
        dim_apply (str, optional) : 
            First order dimension to apply the filter to.
            This was almost always be "time" (default) but is included here for generality.
            
    Returns:
        xarray.Dataset : 
            Copy of ds with data variables (data_vars) filtered.    
    '''
    
    data_vars = ds.data_vars
    ds = ds.copy(deep=True)

    # Apply filter to data variables
    if len(filter_array) > 0:
        for name in data_vars:
            try:
                dim_order = coord_order(ds,data_vars=[name])
                if dim_apply == dim_order["1"][0]:
                    ds[name] = ds[name][filter_array]
            except IndexError:
                raise IndexError("Indices from filter_array are out of bounds for data variable {0}.\n Filter array: {1}".format(name,filter_array))
            except KeyError:
                raise KeyError("Input for filter_array ({0}) should be an array of indices to keep. Current type: {1}".format(filter_array,type(filter_array)))
                
    else:
        # If filter is an empty array, we assume all elements should be removed.
        print('WARNING: Applying this filter has removed all elements.')
        for name in data_vars:
            dim_order = coord_order(ds,data_vars=[name])
            if dim_apply == dim_order["1"][0]:
                ds[name].values.fill(np.nan)

    ds = ds.dropna(dim_apply)
    
    return ds

def coord_order(ds,data_vars=[],check_consistent=False):
    '''
    The coord_order function determines the order of the coordinates for the data variables.
    
    For example for the xarray.Dataset:
        <xarray.Dataset>
            Dimensions:                (level: 3, time: 10)
            Coordinates:
              * time                   (time) datetime64[ns] 2010-01-01T01:30:00 ...
              * level                  (level) int64 0 1 2
            Data variables:
                lat                    (time) float64 nan nan 30.67 31.0 31.33 31.67 ...
                xch4_averaging_kernel  (time, level) float64 nan nan nan nan nan nan 0.0 ...
    
    The two coordinate axes are "time" and "level". This function would determine that "time" is only ever a first
    order coordinate (only used for axis=0) and that "level" is a second order coordinate (only used for axis=1).
    
    Args:
        ds (xarray.Dataset) : 
            Dataset
        data_vars (list, optional) : 
            Which data variables within ds to check.
            If data_vars if left blank ([] by default) all data variables will be examined.
        check_consistent (bool, optional) : 
            If True: checks that each coordinate is only ever used at one order (e.g. time is only
            ever used as a first order coordinate, level as a second order coordinate etc.)
    
    Returns:
        If check_consistent == False:
            dict : orders as a keys each with a list of relevant coordinate names
                    e.g. {1: ["time"], 2: ["level"]}
        If check_consistent == True:
            dict, bool : dictionary as described above, boolean flag for whether order is consistent
    '''
    if not data_vars:
        data_vars = ds.data_vars # Extract data varaible names from Dataset
    dim_order = {} 
    consistent_dims = True
    for name in data_vars:
        dims = ds[name].dims
        for i,d in enumerate(dims):
            order = str(i+1) # Set order as 1 for first order, 2 for second order...
            if order in list(dim_order.keys()): # Check if key already exists
                if d not in dim_order[order]: # Check that dimension is not already present for that key
                    dim_order[order].append(d)
            else:
                dim_order[order] = [d] # Create key if it doesn't exist
            
            if check_consistent:
                for key in list(dim_order.keys()):
                    if key != order: # Check all other keys except where we just added the dimension name
                        if d in dim_order[key]:
                            # If the dimension is being used for different orders by different data variables
                            # the dimensions are not consistent.
                            consistent_dims = False
    
    if check_consistent:
        return dim_order,consistent_dims
    else:
        return dim_order
    
def name_pressure_file(filename,name='surface_pressure',column_names=["latitude","longitude","time"],
                       set_columns=["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"]):
    '''
    The name_pressure_file function extracts the pressure values from one file of a special NAME run to find surface pressure 
    values across a given domain.
    
    NOTE: This function can also be used to extract Topography information from the same NAME run
    To do this set the filename to a topography file (e.g. "Topog_C1_T1_201201010000.txt") and the name to e.g. "topography".
    The set_columns should be the same between pressure and topography files.
    
    Output from:
    - BackRuns_SurfaceField.scr (with SurfaceFieldOutputRequests.txt)
    
    Args:
        filename (str) : 
            Pressure file produced by the appropriate BackRuns script including path information.
            Current script calls this: "Pressure_C1_'YYYYMMDD'_1d.txt" e.g. "Pressure_C1_20120101_15d.txt"
        name (str) : 
            Name to use when creating new pressure (or topography) data variable within output Dataset.
            Default = "surface_pressure".
        column_names (list) : 
            Names for co-ordinate axes within dataset. Should be a 3-item list for latitude, longitude and time.
            Default = ["latitude","longitude","time"]
        set_columns (list) : 
            Column names within input file which contain longitude, latitude and height (in that order).
            Function assumes these columns will be the first in the file and all other populated columns
            contain the surface pressure (or topography) values.
            Default = ["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"] based on current BackRuns script.
    Returns:
        xarray.Dataset : 
            Dataset containing extracted information with dimensions of latitude,longitude,time
    '''
    
    ## Extract data from input file
    lineskip = 35 # 35 lines to skip for a standard NAME file before reaching final header row
    set_columns = ["X (Lat-Long)","Y (Lat-Long)","Z (m agl)"]
    num_set_columns = len(set_columns)
    df = pd.read_csv(filename,skipinitialspace=True,skiprows=lineskip)
        
    ## Remove empty column from data, assume this has been dropeped from the end and compare to original length    
    initial_len = len(df.columns)
    df = df.dropna(axis=1) # Normally contains empty column at the end
    if len(df.columns) != initial_len:
        drop_last = True
    else:
        drop_last = False
    
    ## Extract longitdue and latitude from dataframe and then drop all label columns
    lon = df[set_columns[0]].values
    lat = df[set_columns[1]].values
    df = df.drop(set_columns,axis=1)

    ## Extract pressure values from dataframe and rearrange from row-based to column-based
    pressures = df[:].values
    pressures = np.swapaxes(pressures,1,0) # Swap from row-based to column-based
    
    ## Extract dates lines from file (skip all rows expect this line).
    date_row = 33 # Date row should be on the 33rd (34th?) line
    dates_str = pd.read_csv(filename,skipinitialspace=True,skiprows=date_row-1,nrows=1)
    if drop_last:
        dates_str = dates_str[:].values[0][num_set_columns:-1]
    else:
        dates_str = dates_str[:].values[0][num_set_columns:]
    
    ## Convert dates to datetime object via dt.datetime module; expect date format of: "02/01/2012 00:00 UTC"
    dates = [np.datetime64(dt.datetime.strptime(date.rstrip("UTC"),'%d/%m/%Y %H:%M ')) for date in dates_str]

#   ### Previous method for converting to xarray Dataset but was slow    
#    ## Define axes and elements to create new Multi-Index dataframe - can be converted simply to an xarray dataset
#    num_elements = len(df.index)*len(df.columns)
#    
#    lon_index = np.zeros(num_elements)
#    lat_index = np.zeros(num_elements)
#    dates_index = np.array([np.datetime64("1970-01-01 00:00")]*num_elements,dtype=np.datetime64) # Enter dummy datetime value to initialise array.
#
#    num_rows = len(df.index) # Number of rows for each surface pressure column
#    for i,column in enumerate(pressures):
#        dates_index[i*num_rows:(i+1)*num_rows] = np.array([dates[i]]*num_rows)
#        lon_index[i*num_rows:(i+1)*num_rows] = lon
#        lat_index[i*num_rows:(i+1)*num_rows] = lat
#    
#    pressures_stack = np.ravel(pressures) # Stack all pressure values into one column
#    
#    df_split = pd.DataFrame(pressures_stack,index=[dates_index,lon_index,lat_index]) # Create MultiIndex DataFrame
#
#    ## Create xarray dataset and convert into expected format
#    ds = df_split.to_xarray() # Cast to xarray.dataset
#    ds = ds.rename({"level_0":column_names[2],"level_1":column_names[1],"level_2":column_names[0],0:name})
#    ds = ds.transpose(*column_names) # Rearrange dimensions to match typical footprint order
#   ######
    
    # Lon and Lat columns are presented as a grid of numbers, with one lon value for a set of lat values, so should be safe to do this.
    lon_index = np.unique(lon)
    lat_index = np.unique(lat)
    pressures = np.reshape(pressures,(pressures.shape[0],len(lon_index),len(lat_index)))

    p_col_names = [column_names[2],column_names[1],column_names[0]]
    coords = OrderedDict([(p_col_names[0],np.array(dates)),(p_col_names[1],lon_index),(p_col_names[2],lat_index)])
    
    ds = xray.Dataset(data_vars={name:(p_col_names,pressures)},coords=coords)
    ds = ds.transpose(*column_names) # Rearrange dimensions to match typical footprint order

    if len(np.where(ds["time"] == "1970-01-01 00:00")[0]) != 0:
        print('WARNING: Some time values have not been propagated correctly. Default entries of "1970-01-01 00:00" entries are still present')
        return None
    
    return ds


def extract_files(directory,search_str=None,start=None,end=None,date_separator='',day=True):
    '''
    The extract_files function extracts filenames from a directory based on a search_str and/or based on a start
    and end dates (if filename contains details of a datestamp).
    
    Note:
        By default if using start and end dates files of the form "*YYYYMMDD*" are expected e.g. filename_20111001.csv.
        A date_separator can also be specified between YYYY/MM and MM/DD. 
        For example e.g. date_separator='-', will search for files of the form "*YYYY-MM-DD*".
    
    Args:
        directory (str) : 
            Directory to search (str)
        search_str (str/None, optional) : 
            String to use to search directory (using glob). If full filename is not specified this should 
            contain at least one wildcard character ("*").
            If no search_str specified, all files from directory will be returned and filtered by any start 
            and end criteria.
        start (str/None, optional) : 
            Start date of the form "YYYY-MM-DD".
        end (str/None, optional) : 
            End date of the form "YYYY-MM-DD". Date range is up to but not including this date.
            Either both or neither of start and end should be specified.
        date_separator (str) : 
            Date separator string between year, month and day in filename.
            By default this is "", meaning dates of the form "YYYYMMDD" will be searched for.
        day (bool) :
            Expect day to be included in any date within the filename.
            If day is True expect dates of the form e.g. 20120101 or 2012-01-01.
            If day is False expect dates of the form e.g. 201201 or 2012-01
    
    Returns:
        list : 
            Filenames as string s(full path information).  
    '''
    
    if start:
        if start.find('-') == -1:
            print('WARNING: Start date to extract files is not in correct format: should be YYYY-MM-DD. Unable to set date range.')
            #start=None
            #end=None
            return None
    if end:
        if end.find('-') == -1:
            print('WARNING: End date to extract files is not in correct format: should be YYYY-MM-DD.  Unable to set date range.')
            #start=None
            #end=None
            return None
    
    if (start and not end) or (end and not start):
        print('WARNING: Start and end must both be specified to extract date range of files. Unable to set date range.')
        #start=None
        #end=None
        return None
    
    search_str_short = search_str
    if search_str:
        search_str = os.path.join(directory,search_str)
    else:
        search_str = os.path.join(directory ,"*")

    filenames = glob.glob(search_str)
    filenames.sort()
 
    if start and end:
        #date_range = (np.arange(start,end,dtype="datetime64[D]").astype(str))
        date_range = np.arange(start,end,dtype="datetime64[D]")
        date_range = [str(date) for date in date_range]
        if not day:
           date_range = ['-'.join(date.split('-')[0:2]) for date in date_range] 
        date_range = [date.replace('-',date_separator) for date in date_range]
        print('Finding files in range: {0} - {1} in directory {2} using search string {3}'.format(start,end,directory,search_str_short))
        
        files = []
        for filename in filenames:
            try:
                if date_separator and day:
                    d_sep = "[{0}]".format(date_separator)
                    re_str = r"\d{4}"+d_sep+r"\d{2}"+d_sep+r"\d{2}"  # Creating a regular expression to find 8 digits separated by date_separator e.g. 2012-01-01
                elif date_separator and not day:
                    d_sep = "[{0}]".format(date_separator)
                    re_str = r"\d{4}"+d_sep+r"\d{2}"  # Creating a regular expression to find 6 digits separated by date_separator e.g. 2012-01
                elif day:
                    re_str = r"\d{8}" # Creating a regular expression to find an 8 digit number (should be the date) e.g. 20120101
                elif not day:
                    re_str = r"\d{6}" # Creating a regular expression to find an 6 digit number (should be the date) e.g. 201201
                d = re.search(re_str,filename) 
                d = d.group() # Extract value from regular expression compiler
            except AttributeError:
                pass
            else:
                if d in date_range:
                    files.append(filename)
    else:
        print('Finding files in directory {0} with search str {1}'.format(directory,search_str_short))
        files = filenames
    
    return files


def midpoint_bounds(pressure_levels,delta_pressure=None,pressure_NAME=None,set_as_edges=False,
                    calc_boundaries=True):
    '''
    The midpoint_bounds function uses midpoint between pressure points to define the boundaries. The pressure levels are
    then scaled to be at the centre of that boundary.
    
    Note:
        The assumption is made that first level represents the surface pressure (unless pressure_NAME is specified, where NAME 
        surface pressure value replaces p_0).
    
    Args:
        pressure_levels (numpy.array) : 
            Array containing pressure levels values for each data point (2D array)
        delta_pressure (numpy.array, optional) : 
            Difference between the pressure levels.
        pressure_NAME (numpy.array, optional) : 
            NAME surface pressure values which correspond to pressure_levels latitude, longitude and time coordinates.
            (1D numpy.array). These values can be extracted using the name_pressure_match() function.
            If this is specified, these NAME values will be used as the surface pressure values rather than the 
            first level for each data point within the pressure_levels array.
        set_as_edges (bool, optional) : 
            Whether to treat the pressure levels as points at the edge of the layers rather than within the layers.
            Default = False.
        calc_boundaries (bool, optional) : 
            Whether to apply the necessary special conditions to calculate the bottom and top layers. 
            Note: if this value is set to False, the pressure_levels and dpressure output arrays will not be the 
            same length as the input pressure_levels array.
            Default = True.
    
    Returns:
        (np.array,np.array) : 
            new_pressure_levels,dpressure
    '''

    min_pressure = 0.0
    
    if delta_pressure is not None:
        delta_pressure = pressure_levels[:,:-1] - pressure_levels[:,1:]

    if pressure_NAME is not None:
        use_name_pressure = True
    else:
        use_name_pressure = False
    
    if calc_boundaries:
        dpressure = np.zeros(pressure_levels.shape)
        new_pressure_levels = np.zeros(pressure_levels.shape)
        
        if set_as_edges:
            if use_name_pressure:
                delta_pressure[:,0] = pressure_NAME - pressure_levels[:,1]
                pressure_levels[:,0] = pressure_NAME
            
            dpressure[:,0] = delta_pressure[:,0]
            new_pressure_levels[:,0] = pressure_levels[:,0] - dpressure[:,0]/2. # Define new pressure level from p0 (surface)
            
            dpressure[:,1:-1] = delta_pressure[:,1:]

            dpressure[:,-1] = delta_pressure[:,-1] # Set last dpressure as matching previous
            bound_below_min = np.where(pressure_levels[:,-1] - dpressure[:,-1] < min_pressure)[0]
            dpressure[bound_below_min,-1] = pressure_levels[:,-1] - min_pressure # If pressure would be < min_pressure, set dpressure as distance from min_pressure (e.g. 0.0)

            new_pressure_levels[:,1:] = pressure_levels[:,1:] - dpressure[:,1:]/2.
        else:
            if use_name_pressure:
                delta_pressure[:,0] = pressure_NAME - pressure_levels[:,1]
                pressure_levels[:,0] = pressure_NAME
            
            dpressure[:,0] = delta_pressure[:,0]/2.
            new_pressure_levels[:,0] = pressure_levels[:,0] - dpressure[:,0]/2.
            
            dpressure[:,1:-1] = delta_pressure[:,:-1]/2. + delta_pressure[:,1:]/2.
            
            dpressure[:,-1] = delta_pressure[:,-1] # Set last dpressure as matching previous
            bound_below_min = np.where(pressure_levels[:,-1] - dpressure[:,-1]/2. < min_pressure)[0]
            if len(bound_below_min):
                dpressure[bound_below_min,-1] = pressure_levels[:,-1] + delta_pressure[:,-1]/2. - min_pressure # If pressure would be < min pressure (e.g. 0), set dpressure as distance from min_pressure
            
            new_pressure_levels[:,1:] = (pressure_levels[:,1:] + delta_pressure/2.) - dpressure[:,1:]/2.
    else:
        dpressure = delta_pressure[:,:-1]/2. + delta_pressure[:,1:]/2.
        new_pressure_levels = pressure_levels[:,1:-1] + delta_pressure[:,:-1]/2. - dpressure/2.

    return new_pressure_levels,dpressure
