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

import acrg_agage as agage
import numpy as np
import acrg_config.config as config
from collections import OrderedDict

## Functions below are specifically related to the MCMC code which build on the functions above

def mcmc_param_type():
    '''
    The mcmc_param_type function specifies the names of expected input parameters from the config file and the required Python object types.
    This includes three section groups for the input parameter:
        'MEASUREMENTS' - details related to the measurements made
        'MCMC'         - parameters for running the mcmc model
        'TDMCMC'       - parameters for running the tdmcmc model (in addition to mcmc parameters needed)
    
    
    Returns:
        OrderedDict of section groups:
            Each section group contains a OrderedDict of parameter names and associated type:
            'MEASUREMENTS' ('sites',list),('species',str),('start_date',str),('end_date',str),('domain',str),
                           ('network',str)
            'MCMC'         ('meas_period',list),('av_period',list),('nIt',int),('burn_in',int),
                           ('nsub',int),('fp_basis_case',str),('bc_basis_case',str),('x_pdf0',int),
                           ('pdf_param1_pdf0',int),('pdf_param2_pdf0',int),('sigma_model_pdf',int),
                           ('pdf_param10',float),('pdf_param20',float),('pdf_p1_hparam10',float),
                           ('pdf_p1_hparam20',float),('pdf_p2_hparam10',float),('pdf_p2_hparam20',float),
                           ('sigma_model_ap',float),('sigma_model_hparams',list),('bl_period',int),
                           ('bl_split',bool),('levels',int),('stepsize',float),('stepsize_pdf_p1',float),
                           ('stepsize_pdf_p2',float),('stepsize_sigma_y',float),('stepsize_clon',float),
                           ('stepsize_clat',float),('stepsize_bd',int),('inv_type',str),('tau_ap',float),
                           ('tau_hparams',list),('tau_pdf',int),('stepsize_tau',float), ('filters',list),
                           ('parallel_tempering',bool),('nbeta',int),('output_dir',str),('unique_copy',boolean)
            'TDMCMC'       ('reversible_jump',bool),('kmin',int),('kmax',int),('k_ap',int)
    
    '''
    
    measurements = OrderedDict([('sites',list),
                                ('species',str),
                                ('start_date',str),
                                ('end_date',str),
                                ('domain',str),
                                ('network',str)])
    
    mcmc = OrderedDict([('meas_period',list),
                        ('av_period',list),
                        ('nIt',int),
                        ('burn_in',int),
                        ('nsub',int),
                        ('fp_basis_case',str),
                        ('bc_basis_case',str),
                        ('x_pdf0',int),
                        ('pdf_param1_pdf0',int),
                        ('pdf_param2_pdf0',int),
                        ('sigma_model_pdf',int),
                        ('pdf_param10',float),
                        ('pdf_param20',float),
                        ('pdf_p1_hparam10',float),
                        ('pdf_p1_hparam20',float),
                        ('pdf_p2_hparam10',float),
                        ('pdf_p2_hparam20',float),
                        ('sigma_model_ap',float),
                        ('sigma_model_hparams',np.ndarray),
                        ('bl_period',int),
                        ('bl_split',bool),
                        ('levels',list),
                        ('stepsize',float),
                        ('stepsize_pdf_p1',float),
                        ('stepsize_pdf_p2',float),
                        ('stepsize_sigma_y',float),
                        ('stepsize_clon',float),
                        ('stepsize_clat',float),
                        ('stepsize_bd',int),
                        ('inv_type',str),
                        ('tau_ap',float),
                        ('tau_hparams',np.ndarray),
                        ('tau_pdf',int),
                        ('stepsize_tau',float), 
                        ('filters',list),
                        ('parallel_tempering',bool),
                        ('nbeta',int),
                        ('output_dir',str),
                        ('unique_copy',bool)])
    
    tdmcmc = OrderedDict([('reversible_jump',bool),
                       ('kmin',int),
                       ('kmax',int),
                       ('k_ap',int)])
    
    param_dict = OrderedDict([('MEASUREMENTS', measurements),
                              ('MCMC', mcmc),
                              ('TDMCMC', tdmcmc)])
                          
    
    return param_dict

def get_meas_params():
    '''
    The get_meas_param function returns all parameter names associated the the 'MEASUREMENTS' group
    Returns:
        OrderedDict: parameter names, str: group name
    '''
    key = 'MEASUREMENTS'
    return mcmc_param_type()[key].keys(),key


def get_mcmc_params():
    '''
    The get_meas_param function returns all parameter names associated the the 'MCMC' group
    Returns:
        OrderedDict: parameter names, str: group name
    '''
    key = 'MCMC'
    return mcmc_param_type()[key].keys(),key


def get_tdmcmc_params():
    '''
    The get_meas_param function returns all parameter names associated the the 'TDMCMC' group
    Returns:
        OrderedDict: parameter names, str: group name
    '''
    key = 'TDMCMC'
    return mcmc_param_type()[key].keys(),key
  

def measurements_param(config_file,optional_param=[]):
    '''
    The measurements_param function extracts all parameters relevant to measurement details (see mcmc_param_type for full list)
    
    Args:
        config_file (str)     : filename for input configuration file
        optional_param (list) : parameters which are optional. If the param cannot be found value will be set to None.
    
    Returns:
        OrderedDict: parameter names and values
        
        If any measurement parameter cannot be found (not specified as an optional param):
            Exception raised and program exited
    '''
    
    meas_names,meas_section = get_meas_params()
    can_calculate = ['network'] # Allow values which can be calculated to be optional parameters (currently only 'network')
    optional_param += ['start_date','end_date']
    param_type = mcmc_param_type()
    
    param = config.extract_params(config_file,meas_section,names=meas_names,optional_param=can_calculate+optional_param,param_type=param_type)
    
    for key in can_calculate:
        if key == 'network': # If network is not specified, extract from site_info file
            try:
                param['network']
            except KeyError:
                site1 = param['sites'][0]
                param[key] = agage.site_info[site1]["network"]
                print 'Extracting network for first site from json file'
        else:
            pass # Can add in extra checks here if value can be extracted from another source
    
    #sites = param['sites']
    #species = param['species']
    #domain = param['domain']
    #network = param['network']
    #start_date = param['start_date']
    #end_date = param['end_date']
    
    #return sites,species,domain,network,start_date,end_date
    
    return param

   
def mcmc_param(config_file,optional_param=[]):
    '''
    The mcmc_param function extracts all parameters for the MCMC run (see mcmc_param_type for full list)
     
    Args:
        config_file (str)     : filename for input configuration file
        optional_param (list) : parameters which are optional. If the param cannot be found value will be set to None
    
    Returns:
        OrderedDict: parameter names and values   
        
        If any MCMC parameter cannot be found (not specified as an optional param):
            Exception raised and program exited
    '''
    
    names,section_group = get_mcmc_params()
    optional_param += ['unique_copy']
    param_type = mcmc_param_type()
    
    param = config.extract_params(config_file,section_group=section_group,names=names,optional_param=optional_param,param_type=param_type)
    
    if param['unique_copy'] is None:
        param['unique_copy'] = False
    
    return param

def tdmcmc_param(config_file,optional_param=[]):
    '''
    The tdmcmc_param function extracts all parameters for the MCMC run (see mcmc_param_type for full list)
     
    Args:
        config_file (str)     : filename for input configuration file
        optional_param (list) : parameters which are optional. If the param cannot be found value will be set to None
    
    Returns:
        OrderedDict: parameter names and values   
        
        If any MCMC parameter cannot be found (not specified as an optional param):
            Exception raised and program exited
    '''
    
    names,section_group = get_tdmcmc_params()
    param_type = mcmc_param_type()
    
    param = config.extract_params(config_file,section_group=section_group,names=names,optional_param=optional_param,param_type=param_type)
    
    return param

def all_mcmc_param(config_file,optional_param=[]):
    '''
    The all_mcmc_param function extracts all parameters related to running the tdmcmc code as defined in mcmc_param_type()
    
    Args:
        config_file    : filename for input configuration file
        optional_param : parameters which are optional. If the param cannot be found in input file, value will be set to None
    
    Returns:
        OrderedDict: parameter names and values 
        
        If any parameter defined in mcmc_param_type() cannot be found (not specified as an optional param):
            Exception raised and program exited    
    '''
    
    meas_parameters = measurements_param(config_file,optional_param)
    mcmc_parameters = mcmc_param(config_file,optional_param)
    tdmcmc_parameters = tdmcmc_param(config_file,optional_param)
    
    param = OrderedDict({})
    param.update(meas_parameters)
    param.update(mcmc_parameters)
    param.update(tdmcmc_parameters)
    
    # Checking if any additional keys are present except those explictly used above
    known_keys = [get_meas_params()[1],get_mcmc_params()[1],get_tdmcmc_params()[1]]
    param_type = mcmc_param_type()
    
    for key in param_type:
        if key not in known_keys:
            # Extract the extra parameters but print a warning as these values should really be incorporated into code
            print 'WARNING: Additional unknown key {0} extracted from mcmc_param_type. May be worth adding additional functions for this?'.format(key)
            extra_parameters = config.extract_params(config_file,section_group=key,optional_param=optional_param,param_type=param_type)
            param.update(extra_parameters)
    
    return param