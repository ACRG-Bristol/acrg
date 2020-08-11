#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 12:20:35 2018

This module describes the specific configuration set up for the running the gosat process code.
This builds on the more general functionality of configuration file functions acrg_config.config 
module. 
Please refer to this module for further details.

gosat_process
+++++++++++++

For the gosat_process() inputs an input parameter dictionary is pre-defined (gosat_param_dict()) 
and separates the parameters by section.
This dictionary allows us to define which inputs we expect and the associated types.
Unless specified as optional parameters an error is returned if those parameters aren't present.

See templates/gosat_process_template.ini for the template configuration file on which the parameter
dictionary is based.

Extracted parameters are returned as an OrderedDict object (from collections module).

Updating template file
++++++++++++++++++++++

If you wish to update the format of the gosat configuration input this can be done by updating the 
templates/gosat_process_template.ini file. However, ** this will impact all files of this type ** 
and should only be done if this is the intended outcome.

Note: you may also need to update the gosat_essential_param() or gosat_check() parameters to 
change or add additional checks where necessary.

How to run
+++++++++++

The main function to use for reading in parameters from a config file are:
    
    * gosat_param(config_file) - Extract parameters for input into gosat_process() based on a pre-defined 
    param_type dictionary (see above)

@author: rt17603
"""

import acrg_config.config as config
import os
import sys

if sys.version_info[0] == 2: # If major python version is 2, can't use paths module
    acrg_path = os.getenv("ACRG_PATH") 
else:
    from acrg_config.paths import paths
    acrg_path = paths.acrg

def gosat_param_dict():
    '''
    The gosat_param_dict function defines the nested dictionary for the GOSAT process input parameters to be passed to
    the gosat_process function.
    Currently this format is based on "gosat_process_template.ini" and this file should not be altered unless you wish to
    change the expected parameters and types for all parameter files of this type.
    
    Returns
        param_type nested dictionary
    '''
   
#    param_dict = {"MEASUREMENTS":
#                     {"input_directory":str,
#                      "species":str},
#                  "MEASUREMENTS.SELECTION":
#                      {"site":str,
#                       "lat_bounds":list,
#                       "lon_bounds":list,
#                       "domain":str,
#                       "coord_bin":list,
#                       "start":str,
#                       "end":str},
#                  "MEASUREMENTS.FILTER":
#                      {"quality_filt":bool,
#                       "bad_pressure_filt":bool,
#                       "mode":str},
#                  "MEASUREMENTS.NAME_SP_FILT":   
#                     {"name_sp_filt":bool,
#                      "name_filters":list,
#                      "cutoff":float,
#                      "layer_range":list},
#                  "NAME.SURFACE_PRESSURE":
#                      {"use_name_pressure":bool,
#                       "pressure_dir":str},
#                  "NAME.OUTPUT":
#                      {"write_name":bool,
#                       "name_file_per_day":bool,
#                       "name_directory":str},
#                  "NC.OUTPUT":
#                      {"write_nc":bool,
#                       "output_directory":str},
#                  "OUTPUT":
#                      {"overwrite":bool}
#                    }
    
    reference_file = os.path.join(acrg_path,"acrg_config/templates/gosat_process_template.ini")
    param_dict = config.generate_param_dict(reference_file)
    
    return param_dict  

def gosat_essential_param():
    '''
    The gosat_essential_param function defines the parameters which *must* be specified in the 
    gosat input parameter file.
    
    All other parameters will be set as optional.
    
    Returns:
        list: parameter names (str)
    '''
    essential_param = ["site"]
    return essential_param

def gosat_check(param):
    '''
    The gosat_check function applies additional checks on parameters within the param dictionary (extracted from a 
    configuration file) to ensure they are of the correct format.
    
    This includes:
        Check other necessary parameters are specified when certain conditions are set to true
            If 'write_name' is specified, ensure 'name_directory' is not set to "/path/to/output/directory/"
            If 'write_nc' is specified, ensure 'output_directory' is not set to "/path/to/output/directory/"
    
    Args:
        param (dict) : output of gosat_param function. Dictionary containing parameter and value details.
    
    Returns:
        None
        
        Raises Error if any of the above conditions are not met.
    '''

    path_check = "/path/to/output/directory/"
    if 'write_name' in param.keys() and 'name_directory' in param.keys():
        if param['write_name'] == True:
            if param['name_directory'] == path_check:
                raise Exception('Please set name_directory parameter or remove option from input param file to use the default. Currently set to {}'.format(path_check))
    
    if 'write_nc' in param.keys() and 'output_directory' in param.keys():
        if param['write_nc'] == True:
            if param['output_directory'] == path_check:
                raise Exception('Please set output_directory parameter or remove option from input param file to use the default. Currently set to {}'.format(path_check))
    
#    if 'write_name' in param.keys():
#        if param['write_name'] == True:
#            if 'name_directory' not in param.keys():
#                raise Exception('Output directory for NAME files must be specified when write_name=True')
#            elif not param['name_directory']:
#                raise Exception('Output directory for NAME files must contain a value when write_name=True')
#    if 'write_nc' in param.keys():
#        if param['write_nc'] == True:
#            if 'output_directory' not in param.keys():
#                raise Exception('Output directory for netCDF files must be specified when write_nc=True')
#            elif not param['output_directory']:
#                raise Exception('Output directory for netCDF files must contain a value when write_nc=True')
    
def gosat_param(config_file):
    '''
    The gosat_param function reads the parameters for input into the acrg_gosat.gosat_process(...) function.
    If param = gosat_param(config_file), the acrg_gosat.gosat_process function can be called as:
        acrg_gosat.gosat_process(**param)
    
    Note: Parameters which must always be specified in the input configuration file are defined by 
    gosat_essential_param() function.
    Note: Additional checks on the format and inclusion of some parameters will be performed by the gosat_check() 
    function.
    
    Args:
        config_file (str) : Configuration filename. Should match format of gosat_process_template.ini file.
    
    Returns:
        collections.OreredDict : dictionary of parameters from config file.
    '''
    
    param_dict = gosat_param_dict()
    
    #all_params = config.all_parameters_in_param_type(param_dict)
    essential_param = gosat_essential_param()
    #optional_param = all_params[:]
    #for p in essential_param:
    #    optional_param.remove(p)
    
    param = config.all_param(config_file,
                             #optional_param=optional_param,
                             expected_param=essential_param,
                             param_type=param_dict,
                             exclude_not_found=True)
    gosat_check(param)
    
    return param
