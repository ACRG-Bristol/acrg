#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 12:41:21 2018

This module describes the specific configuration set up for the running the tdmcmc code.
This builds on the more general functionality of configuration file functions acrg_config.config module. 
Please refer to this module for further details.

MCMC
++++

For the MCMC inputs an input parameter dictionary is pre-defined (mcmc_param_type()) and separates the parameters by classifications.
This dictionary allows us to define which inputs we expect and the associated types.
Unless specified as optional parameters an error is returned if those parameters aren't present.

By default this includes three classifications:
    - MEASUREMENTS - defines expected parameters from [MEASUREMENTS] section
    - MCMC         - defines expected parameters from all [MCMC.***] sections collected together
    - TDMCMC       - defines expected parameters from all [TDMCMC.***] sections collected together
See mcmc_param_type() function for full list of pre-defined parameters and types

Extracted parameters are returned as an OrderedDict object (from collections module).

How to run
++++++++++

The main function to use for reading in parameters from a config file are:
    
    * all_mcmc_param(config_file,...) - Extract MCMC parameters specifically based on a pre-defined param_type 
    dictionary (see above)

Note: At present the param_type dict specified by the gosat_param_dict() function is created from an in-built 
nested dictionary. However, this can be modified to rely on a template configuration file instead. See
acrg_config.generate_param_dict() function for more details. 

@author: rt17603
"""

import os
import acrg_agage as agage
import numpy as np
import acrg_config.config as config
from collections import OrderedDict

acrg_path = os.getenv("ACRG_PATH")

## Functions below are specifically related to the MCMC code which build on the functions within acrg_config module

def mcmc_param_type():
    '''
    The mcmc_param_type function defines the names of expected input parameters from the config file and 
    the required Python object types.
    
    The mcmc_param_type output is currently based on an input template configuration file:
        $ACRG_PATH/acrg_config/templates/tdmcmc_template.ini
    
    WARNING: The description that follows may be subject to change, always check the template config file.
    The section headings are split into three section groups for the input parameters:
        'MEASUREMENTS' - details related to the measurements made
        'MCMC'         - parameters for running the mcmc model
        'TDMCMC'       - parameters for running the tdmcmc model (in addition to mcmc parameters needed)
    
    Returns:
        nested OrderedDict of sections and parameters within the configuration file

    e.g.    
    OrderedDict([('MEASUREMENTS',
              OrderedDict([('sites', list),
                           ('species', str),
                           ('start_date', str),
                           ('end_date', str),
                           ('domain', str),
                           ('network', str)])),
             ('MCMC.MEASUREMENTS',
              OrderedDict([('meas_period', list), ('av_period', list)])),
             ('MCMC.ITERATIONS',
              OrderedDict([('nIt', int), ('burn_in', int), ('nsub', int)])),
             ('MCMC.BASIS_FUNCTIONS',
              OrderedDict([('fp_basis_case', str), ('bc_basis_case', str)])),
             ('TDMCMC.SET_UP',
              OrderedDict([('reversible_jump', bool),
                           ('kmin', int),
                           ('kmax', int),
                           ('k_ap', int)])),
             ('MCMC.PDFS',
              OrderedDict([('x_pdf0', int),
                           ('pdf_param1_pdf0', int),
                           ('pdf_param2_pdf0', int),
                           ('sigma_model_pdf', int)])),
             ('MCMC.HYPERPARAMETERS',
              OrderedDict([('pdf_param10', float),
                           ('pdf_param20', float),
                           ('pdf_p1_hparam10', float),
                           ('pdf_p1_hparam20', int),
                           ('pdf_p2_hparam10', float),
                           ('pdf_p2_hparam20', int)])),
             ('MCMC.MODEL_UNCERTAINTY',
              OrderedDict([('sigma_model_ap', float),
                           ('sigma_model_hparams', numpy.ndarray),
                           ('bl_period', int),
                           ('bl_split', bool),
                           ('levels', int)])),
             ('MCMC.STEPSIZE',
              OrderedDict([('stepsize', float),
                           ('stepsize_pdf_p1', float),
                           ('stepsize_pdf_p2', float),
                           ('stepsize_sigma_y', float),
                           ('stepsize_clon', float),
                           ('stepsize_clat', float),
                           ('stepsize_bd', int)])),
             ('MCMC.COVARIANCE_MATRIX',
              OrderedDict([('inv_type', str),
                           ('tau_ap', float),
                           ('tau_hparams', numpy.ndarray),
                           ('tau_pdf', int),
                           ('stepsize_tau', float)])),
             ('MCMC.FILTERS', OrderedDict([('filters', list)])),
             ('MCMC.PARALLEL_TEMPERING',
              OrderedDict([('parallel_tempering', bool), ('nbeta', int)])),
             ('MCMC.OUTPUT',
              OrderedDict([('output_dir', str), ('unique_copy', bool)]))])

    '''
    
#    measurements = OrderedDict([('sites',list),
#                                ('species',str),
#                                ('start_date',str),
#                                ('end_date',str),
#                                ('domain',str),
#                                ('network',str)])
#    
#    mcmc = OrderedDict([('meas_period',list),
#                        ('av_period',list),
#                        ('nIt',int),
#                        ('burn_in',int),
#                        ('nsub',int),
#                        ('fp_basis_case',str),
#                        ('bc_basis_case',str),
#                        ('x_pdf0',int),
#                        ('pdf_param1_pdf0',int),
#                        ('pdf_param2_pdf0',int),
#                        ('sigma_model_pdf',int),
#                        ('pdf_param10',float),
#                        ('pdf_param20',float),
#                        ('pdf_p1_hparam10',float),
#                        ('pdf_p1_hparam20',float),
#                        ('pdf_p2_hparam10',float),
#                        ('pdf_p2_hparam20',float),
#                        ('sigma_model_ap',float),
#                        ('sigma_model_hparams',np.ndarray),
#                        ('bl_period',int),
#                        ('bl_split',bool),
#                        ('levels',list),
#                        ('stepsize',float),
#                        ('stepsize_pdf_p1',float),
#                        ('stepsize_pdf_p2',float),
#                        ('stepsize_sigma_y',float),
#                        ('stepsize_clon',float),
#                        ('stepsize_clat',float),
#                        ('stepsize_bd',int),
#                        ('inv_type',str),
#                        ('tau_ap',float),
#                        ('tau_hparams',np.ndarray),
#                        ('tau_pdf',int),
#                        ('stepsize_tau',float), 
#                        ('filters',list),
#                        ('parallel_tempering',bool),
#                        ('nbeta',int),
#                        ('output_dir',str),
#                        ('unique_copy',bool)])
#    
#    tdmcmc = OrderedDict([('reversible_jump',bool),
#                       ('kmin',int),
#                       ('kmax',int),
#                       ('k_ap',int)])
#    
#    param_dict = OrderedDict([('MEASUREMENTS', measurements),
#                              ('MCMC', mcmc),
#                              ('TDMCMC', tdmcmc)])
    reference_file = os.path.join(acrg_path,"acrg_config/templates/tdmcmc_template.ini")
    param_dict = config.generate_param_dict(reference_file)
    
    return param_dict

#def get_meas_params():
#    '''
#    The get_meas_param function returns all parameter names associated the the 'MEASUREMENTS' group
#    Returns:
#        OrderedDict: parameter names, str: group name
#    '''
#    key = 'MEASUREMENTS'
#    return mcmc_param_type()[key].keys(),key
#
#
#def get_mcmc_params():
#    '''
#    The get_meas_param function returns all parameter names associated the the 'MCMC' group
#    Returns:
#        OrderedDict: parameter names, str: group name
#    '''
#    key = 'MCMC'
#    return mcmc_param_type()[key].keys(),key
#
#
#def get_tdmcmc_params():
#    '''
#    The get_meas_param function returns all parameter names associated the the 'TDMCMC' group
#    Returns:
#        OrderedDict: parameter names, str: group name
#    '''
#    key = 'TDMCMC'
#    return mcmc_param_type()[key].keys(),key
  

def optional_parameters(section_group=None):
    '''
    The optional_parameters function defines the set of parameters whcih do not need to be specified within the 
    configuration file. (All other parameters within dictionary from mcmc_param_type() MUST be included).
    This includes some parameters which can be calculated independently without an output or have a default value.

    TODO: Add any additional optional parameters to this function as the configuration set-up changes
    and expands.
    
    Args:
        section_group (str/None, optional):
            Allows selection of optional parameters for a particular section_group.
            Should be one of: "MEASUREMENTS","MCMC","TDMCMC"
            If None is specified, optional parameters for all section_groups will be returned
    
    Returns:
        list:
            Optional parameters for MCMC code
    '''
    meas_params = ["network","start_date","end_date"]
    mcmc_params = ["unique_copy"]
    tdmcmc_params = []
    
    optional_param = []
    if section_group is None or section_group == "MEASUREMENTS":
        optional_param += meas_params
    if section_group is None or section_group == "MCMC":
        optional_param += mcmc_params
    if section_group is None or section_group == "TDMCMC":
        optional_param += tdmcmc_params
    
    return optional_param

def add_defaults(param,section_group=None):
    '''
    The add_defaults function adds any default values if they have not been explictly specified within
    the input configuation file.
    
    The current defaults that can be added are:
        "network" (within MEASUREMENTS section group) - extracted from acrg_site_info.json file
        "unique_copy" (within MCMC section group)     - set to False as a default
    
    TODO: Add additional defaults as the configuration set-up changes and expands.
    
    Args:
        param (dict) :
            Dictionary of parameter names and values.
            Output from all_mcmc_param, measurements_param, mcmc_param, tdmcmc_param functions.
        section_group (str/None, optional):
            Only apply defaults for the particular section_group.
            Should be one of: "MEASUREMENTS","MCMC","TDMCMC"
            If None is specified, defaults will be added for all section_groups
    
    Returns:
        dict:
            param dictionary with default values added where necessary
    
    '''
    if section_group is None or section_group == "MEASUREMENTS":
        if ("network" not in param.keys()) or (not param["network"]):
            site1 = param['sites'][0]
            param["network"] = agage.site_info[site1]["network"]
            print 'Extracting network for first site from json file'
    
    if section_group is None or section_group == "MCMC":
        if ("unique_copy" not in param.keys()) or (param["unique_copy"] == None):
            param["unique_copy"] = False

    return param

def measurements_param(config_file,optional_param=[]):
    '''
    The measurements_param function extracts all parameters relevant to measurement details (see mcmc_param_type for full list)

    In addition to optional_param values, parameters from optional_parameters() function for section_group="MEASUREMENTS" 
    are also treated as optional when reading the configuration file.
     
    Args:
        config_file (str) : 
            Filename for input configuration file
        optional_param (list, optional) : 
            Additional parameters within configuration file which are optional. If the parameter cannot be found in input file, 
            value will set to None.
    
    Returns:
        OrderedDict: 
            parameter names and values
        
        If any measurement parameter cannot be found (not specified as an optional param):
            Exception raised and program exited
    '''
    
    meas_group = 'MEASUREMENTS'
    
    optional_param += optional_parameters(section_group=meas_group)
    param_type = mcmc_param_type()

    param = config.extract_params(config_file,section_group=meas_group,optional_param=optional_param,param_type=param_type)
    
    param = add_defaults(param,section_group=meas_group)
    
    return param

   
def mcmc_param(config_file,optional_param=[]):
    '''
    The mcmc_param function extracts all parameters for the MCMC run (see mcmc_param_type for full list)

    In addition to optional_param values, parameters from optional_parameters() function for section_group="MCMC" are 
    also treated as optional when reading the configuration file.
     
    Args:
        config_file (str) : 
            Filename for input configuration file
        optional_param (list, optional) : 
            Additional parameters within configuration file which are optional. If the parameter cannot be found in input file, 
            value will set to None.
    
    Returns:
        OrderedDict: 
            parameter names and values   
        
        If any MCMC parameter cannot be found (not specified as an optional param):
            Exception raised and program exited
    '''
    
    mcmc_group = "MCMC"

    optional_param += optional_parameters(section_group=mcmc_group)
    param_type = mcmc_param_type()
    
    param = config.extract_params(config_file,section_group=mcmc_group,optional_param=optional_param,param_type=param_type)
    
    param = add_defaults(param,section_group=mcmc_group)
    
    return param

def tdmcmc_param(config_file,optional_param=[]):
    '''
    The tdmcmc_param function extracts all parameters for the MCMC run (see mcmc_param_type for full list)

    In addition to optional_param values, parameters from optional_parameters() function for section_group="TDMCMC" are 
    also treated as optional when reading the configuration file.
     
    Args:
        config_file (str) : 
            Filename for input configuration file
        optional_param (list, optional) : 
            Additional parameters within configuration file which are optional. If the parameter cannot be found in input file, 
            value will set to None.
    
    Returns:
        OrderedDict: 
            parameter names and values   
        
        If any MCMC parameter cannot be found (not specified as an optional param):
            Exception raised and program exited
    '''
    
    tdmcmc_group = "TDMCMC"

    optional_param += optional_parameters(section_group=tdmcmc_group)
    param_type = mcmc_param_type()
    
    param = config.extract_params(config_file,section_group=tdmcmc_group,optional_param=optional_param,param_type=param_type)
    
    param = add_defaults(param,section_group=tdmcmc_group)
    
    return param

def all_mcmc_param(config_file,optional_param=[]):
    '''
    The all_mcmc_param function extracts all parameters related to running the tdmcmc code as defined in mcmc_param_type()
    
    In addition to optional_param values, parameters from optional_parameters() function are also treated as optional 
    when reading the configuration file.
    
    Args:
        config_file (str) : 
            Filename for input configuration file
        optional_param (list, optional) : 
            Additional parameters within configuration file which are optional. If the parameter cannot be found in input file, 
            value will set to None.

    Returns:
        OrderedDict: 
            parameter names and values 
        
        If any parameter defined in mcmc_param_type() cannot be found (not specified as an optional param):
            Exception raised and program exited    
    '''
    
    param_type = mcmc_param_type()
    optional_param += optional_parameters()
    param = config.extract_params(config_file,optional_param=optional_param,param_type=param_type)
    param = add_defaults(param)
    
#    meas_parameters = measurements_param(config_file,optional_param)
#    mcmc_parameters = mcmc_param(config_file,optional_param)
#    tdmcmc_parameters = tdmcmc_param(config_file,optional_param)
    
#    param = OrderedDict({})
#    param.update(meas_parameters)
#    param.update(mcmc_parameters)
#    param.update(tdmcmc_parameters)
    
#    # Checking if any additional keys are present except those explictly used above
#    known_keys = [get_meas_params()[1],get_mcmc_params()[1],get_tdmcmc_params()[1]]
#    param_type = mcmc_param_type()
#    
#    for key in param_type:
#        if key not in known_keys:
#            # Extract the extra parameters but print a warning as these values should really be incorporated into code
#            print 'WARNING: Additional unknown key {0} extracted from mcmc_param_type. May be worth adding additional functions for this?'.format(key)
#            extra_parameters = config.extract_params(config_file,section_group=key,optional_param=optional_param,param_type=param_type)
#            param.update(extra_parameters)
    
    return param