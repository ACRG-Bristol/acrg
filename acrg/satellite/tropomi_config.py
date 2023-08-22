#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 13:28:05 2021

@author: rt17603
"""

import os
import sys
from pathlib import Path
import acrg.config.config as config

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH") 
else:
    from acrg.config.paths import Paths
    acrg_path = Paths.acrg
    
def tropomi_template_file():
    '''
    Define template file for inputs to tropomi.tropomi_process() function.
    '''
    reference_file = os.path.join(acrg_path,"acrg/config/templates/tropomi_process_template.ini")
    return reference_file

def tropomi_param_dict():
    '''
    This function defines the nested dictionary for the TROPOMI process 
    input parameters to be passed to the tropomi_process function.
    Currently this format is based on "tropomi_process_template.ini" and 
    this file should not be altered unless you wish to change the 
    expected parameters and types for all parameter files of this type.
    
    Returns
        param_type nested dictionary
    '''
   
    reference_file = tropomi_template_file()
    param_dict = config.generate_param_dict(reference_file)
    
    return param_dict  

def tropomi_essential_param():
    '''
    This function defines the parameters which *must* be specified 
    in the input parameter file.
    
    All other parameters will are allowed to be optional.
    
    Returns:
        list: parameter names (str)
    '''
    essential_param = ["site", "start_date", "end_date", 
                       "lat_bounds", "lon_bounds",
                       "coord_bin", "regrid_method",
                       "output_directory"]
    
    ##TODO: Add domain input as alternative to "lat_bounds",
    # "lon_bounds" and "coord_bin" ??
    # Would need to change essential parametes and incorporate
    # into tropomi_check() function.
    
    return essential_param

def tropomi_check(param):
    '''
    This function applies additional checks on parameters within the 
    param dictionary (extracted from a configuration file) to ensure 
    they are of the correct format.
    
    This includes:
        Check other necessary parameters are specified when certain conditions are set to true
            If 'write_name' is specified, ensure 'name_directory' is not set to "/path/to/output/directory/"
    
    Args:
        param (dict) : output of gosat_param function. Dictionary containing parameter and value details.
    
    Returns:
        None
        
        Raises Error if any of the above conditions are not met.
    '''

    path_check = "/path/to/output/directory/"
    
    if 'write_name' in param.keys() and 'name_directory' in param.keys():
        if param['write_name'] == True:
            if param['name_directory'] == path_check or not param['name_directory']:
                raise Exception('Please set name_directory parameter or remove option from input param file to use the default. Currently set to {}'.format(path_check))

    # TODO: Add this back in, if you make writing the netcdf file an 
    # optional output.
#    if 'write_nc' in param.keys() and 'output_directory' in param.keys():
#        if param['write_nc'] == True:
#            if param['output_directory'] == path_check:
#                raise Exception('Please set output_directory parameter or remove option from input param file to use the default. Currently set to {}'.format(path_check))

    # TODO: If use_surface_pressure is included check "pressure_domain" is set?
    
    # For any parameters which could be passed as empty strings
    # do an extra check to make sure these are defined if specified.
    param_not_empty = ["regrid_method", "time_increment",
                       "start_date", "end_date",
                       "input_directory","output_directory"]
    for p in param_not_empty:
        if p in param:
            if not param[p]:
                raise Exception(f"If '{p}' is included in configuration file, please make sure this is defined.")
    
    # paths = ["input_directory","output_directory","name_directory"]
    # for path in paths:
    #     if path in param:
    #         param[path] = Path(param[path])
    
    #regrid_method
    #time_increment
    #output_directory
    
    return param
    

    
def tropomi_param(config_file,**command_line):
    '''
    Extract parameters from input configuration file and associated MCMC function. 
    Checks the mcmc_type to extract the required parameters.
    
    Args:
        config_file (str):
            Configuration file name. Should be an .ini file.
        command_line :
            Any additional command line arguments to be added to the param dictionary or to superceed 
            values contained within the config file.
    
    Returns:
        collections.OrderedDict :
            Ordered dictionary of parameter names and values which can be passed to tropomi_process function
    '''

    param_dict = tropomi_param_dict()
    expected_param = tropomi_essential_param()

    # If an expected parameter has been passed from the command line, this does not need to be within the config file
    for key,value in command_line.items():
        if key in expected_param and value is not None:
            expected_param.remove(key)

    param = config.extract_params(config_file,
                                  expected_param=expected_param,
                                  param_type=param_dict,
                                  exclude_not_found=True)

    # Command line values added to param (or superceed inputs from the config file)
    for key,value in command_line.items():
        if value is not None:
            param[key] = value

    # Apply additional parameter checks
    param = tropomi_check(param)
    
    # If configuration file does not include values for the required parameters - produce an error
    for ep in expected_param:
        if not param[ep]:
            raise Exception(f"Required parameter '{ep}' has not been defined")

    return param 
