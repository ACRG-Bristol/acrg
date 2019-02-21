#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 14:19:35 2017

Test suite for acrg_config module.

To run this test suite only from within the tests/ directory use the syntax
>> pytest test_acrg_config.py

Some tests have an additional decorator:
    @pytest.mark.basic
which marks these tests as being part of a "basic" subset to be run if a full test run is not required.

To run only the basic tests use the syntax
>> pytest test_acrg_config.py -m basic

@author: rt17603
"""

from builtins import zip
from builtins import str
import numpy as np
import os
import pytest

import acrg_config.config as configread
from acrg_config.config import extract_params,all_param
import acrg_tdmcmc.tdmcmc_config as config_tdmcmc
from acrg_tdmcmc.tdmcmc_config import all_mcmc_param

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")
 
test_config_path = os.path.join(acrg_path,"tests/files/config")

@pytest.fixture(scope="module")
def tdmcmc_config():
    ''' Define tdmcmc config file input '''
    filename = os.path.join(test_config_path,'tdmcmc_input_1.ini')
    return filename

@pytest.fixture(scope="module")
def tdmcmc_config_option():
    ''' Define tdmcmc config file input with optional parameters not included '''
    filename = os.path.join(test_config_path,'tdmcmc_input_2.ini')
    return filename

@pytest.fixture(scope="module")
def tdmcmc_config_missing():
    ''' Define tdmcmc config file input with optional parameters and an expected parameter not included '''
    filename = os.path.join(test_config_path,'tdmcmc_input_3.ini')
    return filename

@pytest.fixture(scope="module")
def tdmcmc_config_extra():
    ''' Define tdmcmc config file input with extra parameters included '''
    filename = os.path.join(test_config_path,'tdmcmc_input_4.ini')
    return filename

@pytest.fixture(scope="module")
def tdmcmc_config_missing_section():
    ''' Define tdmcmc config file input with optional section not included '''
    filename = os.path.join(test_config_path,'tdmcmc_input_5.ini')
    return filename

@pytest.fixture(scope="module")
def mcmc_param_type():
    ''' Define MCMC parameter type nested dictionary for tdmcmc config input '''
    alt_filename = os.path.join(acrg_path,"tests/files/config/tdmcmc_template.ini")
    param_type = config_tdmcmc.mcmc_param_type(alt_filename)
    return param_type

@pytest.fixture(scope="module")
def example_config():
    ''' Define example config file input '''
    filename = os.path.join(test_config_path,'example_input.ini')
    return filename

@pytest.fixture(scope="module")
def example_param_type():
    ''' Define parameter type nested dictionary for example config file '''
    section_1 = {'parameter_1':float,
                 'parameter_2':int,
                 'parameter_3':list}
    
    section_2 = {'parameter_4':str,
                 'parameter_5':np.ndarray}
    
    section_3 = {'parameter_6':bool}
    
    param_type_1 = {'SECTION':section_1,
                    'NAME.OTHER_SECTION':section_2,
                    'NAME.THIRD_SECTION':section_3}
    
    return param_type_1

@pytest.mark.basic    
def test_extract_all(tdmcmc_config):
    ''' Test parameters can be extracted from input file with no other inputs '''
    x = extract_params(tdmcmc_config)
    assert x

@pytest.mark.basic
def test_extract_section_1(tdmcmc_config):
    ''' Test a simple form section can be specified '''
    x = extract_params(tdmcmc_config,section='MEASUREMENTS',param_type=None)
    assert x

@pytest.mark.basic
def test_extract_section_2(tdmcmc_config):
    ''' Test a more complex form section can be specified '''
    x = extract_params(tdmcmc_config,section='MCMC.MEASUREMENTS',param_type=None)
    assert x        

@pytest.mark.basic
def test_extract_section_group_1(tdmcmc_config):
    ''' Test a section group can be specified '''
    x = extract_params(tdmcmc_config,section_group='MEASUREMENTS',param_type=None)
    assert x

def test_extract_section_group_2(tdmcmc_config):
    ''' Test a different section group can be specified '''
    x = extract_params(tdmcmc_config,section_group='MCMC',param_type=None)
    assert x

@pytest.mark.basic
def test_extract_names(tdmcmc_config):
    ''' Test specific parameter names can be specified and extracted '''
    x = extract_params(tdmcmc_config,names=['sites'],param_type=None)
    assert x

def test_extract_names_section(tdmcmc_config):
    ''' Test specific parameter names can be found within the specified section '''
    x = extract_params(tdmcmc_config,section='MEASUREMENTS',names=['sites','species','domain'],param_type=None)
    assert x

def test_extract_names_section_group(tdmcmc_config):
    ''' Test parameter names can be found within a section covered by the specified section group '''
    x = extract_params(tdmcmc_config,section_group='MCMC',names=['meas_period','nIt'],param_type=None)
    assert x

def test_extract_names_WrongSectionForName(tdmcmc_config):
    ''' Test an empty object is returned if parameter names cannot be found within the input section.'''
    x = extract_params(tdmcmc_config,section='MEASUREMENTS',names=['meas_period','nIt'],param_type=None)
    assert not x

def test_extract_names_WrongSectionGroupForName(tdmcmc_config):
    ''' Test an empty object is returned if parameter names cannot be found within the input sections within the section group. '''
    x = extract_params(tdmcmc_config,section_group='MCMC',names=['sites','species','domain'],param_type=None)
    assert not x

def test_extract_IncorrectSection(tdmcmc_config):
    ''' Test an error is raised if the specified section cannot be found '''
    with pytest.raises(KeyError) as e_info:
        extract_params(tdmcmc_config,section='FLIBBLE',param_type=None)
    #assert not x

def test_extract_IncorrectSectionGroup(tdmcmc_config):
    ''' Test an error is raised when no sections are found for the specified section group. '''
    with pytest.raises(KeyError) as e_info:
        extract_params(tdmcmc_config,section_group='FLIBBLE',param_type=None)
    #assert not x

#@pytest.mark.basic
#def test_extract_all_param_type(tdmcmc_config,mcmc_param_type):
#    ''' Test correct types are created when extracting using a param_type dictionary '''
#    x = extract_params(tdmcmc_config,param_type=mcmc_param_type)
#    
#    parameters = []
#    types = []
#    for section_group in mcmc_param_type.values():
#        parameters.extend(section_group.keys())
#        types.extend(section_group.values())
#    
#    for param,t in zip(parameters,types):
#        assert isinstance(x[param],t)

@pytest.mark.basic
def test_extract_all_param_type(example_config):
    ''' Test correct types are created when extracting using a param_type dictionary '''
    
    param_type = {"SECTION":
                      {"parameter_1":float,
                       "parameter_2":int,
                       "parameter_3":list},
                  "NAME.OTHER_SECTION":
                      {"parameter_4":float,
                       "parameter_5":list},
                  "NAME.THIRD_SECTION":
                      {"parameter_6":bool}
                  }
    
    x = extract_params(example_config,param_type=param_type)
    
    parameters = []
    types = []
    for section in list(param_type.values()):
        parameters.extend(list(section.keys()))
        types.extend(list(section.values()))
    
    for param,t in zip(parameters,types):
        assert isinstance(x[param],t)

def check_types(param_output,param_type,section=None,section_group=None):
    ''' Function to check correct parameters and types are being returned '''
    if section:
        spec = 'section'
        key = section
        section_group = section.split('.')[0]
        param_spec = configread.find_param_key(param_type,section)[1]
    elif section_group:
        spec = 'section_group'
        key = section_group
        param_spec = configread.find_param_key(param_type,section_group)[1]

    if param_spec == spec:
        parameters = list(param_type[key].keys())
        types = list(param_type[key].values())
        for param,t in zip(parameters,types):
            assert isinstance(param_output[param],t)
    elif (param_spec == 'section_group') and (spec == 'section'):
        all_parameters = list(param_type[section_group].keys())
        all_types = list(param_type[section_group].values())
        parameters = list(param_output.keys())
        for param,t in zip(all_parameters,all_types):
            if param in parameters:
                assert isinstance(param_output[param],t)
    else:
        assert param_output

@pytest.mark.basic
def test_extract_section_1_param_type(tdmcmc_config,mcmc_param_type):
    ''' Test a parameter type dictionary can be used when extracting a specific simple form section '''
    x = extract_params(tdmcmc_config,section='MEASUREMENTS',param_type=mcmc_param_type)
    check_types(x,mcmc_param_type,section='MEASUREMENTS')

def test_extract_section_2_param_type(tdmcmc_config,mcmc_param_type):
    ''' Test a parameter type dictionary can be used to extracting a specific more complex form section '''
    x = extract_params(tdmcmc_config,section='MCMC.MEASUREMENTS',param_type=mcmc_param_type)
    check_types(x,mcmc_param_type,section='MCMC.MEASUREMENTS')        

def test_extract_section_group_1_param_type(tdmcmc_config,mcmc_param_type):
    ''' Test a parameter type dictionary can be used to extracting a sections covered by a specific section group '''
    x = extract_params(tdmcmc_config,section_group='MEASUREMENTS',param_type=mcmc_param_type)
    check_types(x,mcmc_param_type,section_group='MEASUREMENTS')

def test_extract_section_group_2_param_type(tdmcmc_config,mcmc_param_type):
    ''' Test a parameter type dictionary can be used to extracting a sections covered by a specific section group '''
    x = extract_params(tdmcmc_config,section_group='MCMC',param_type=mcmc_param_type)
    check_types(x,mcmc_param_type,section_group='MCMC')

@pytest.mark.basic
def test_extract_param_optional(tdmcmc_config_option,mcmc_param_type):
    ''' Test optional parameters can be specified for the whole config file '''
    x = extract_params(tdmcmc_config_option,optional_param=['network','unique_copy','emissions_name'],param_type=mcmc_param_type)
    assert x

def test_extract_param_section_optional(tdmcmc_config_option,mcmc_param_type):
    ''' Test optional parameters can be specified for one section '''    
    x = extract_params(tdmcmc_config_option,section_group='MEASUREMENTS',optional_param=['network','emissions_name'],param_type=mcmc_param_type)
    assert x

def test_extract_WrongSectionForName_param_type(tdmcmc_config,mcmc_param_type):
    ''' Test an error is raised if parameter names cannot be found within the input section when param type dictionary is specified. '''
    with pytest.raises(KeyError) as e_info:
        extract_params(tdmcmc_config,section='MEASUREMENTS',names=['meas_period','nIt'],param_type=mcmc_param_type)

def test_extract_WrongSectionGroupForName_param_type(tdmcmc_config,mcmc_param_type):
    ''' Test an error is raised if parameter names cannot be found within a section covered by the the input section group
    when param type dictionary is specified. '''
    with pytest.raises(KeyError) as e_info:
        extract_params(tdmcmc_config,section_group='MCMC',names=['sites','species','domain'],param_type=mcmc_param_type)

def test_extract_NameNotInFile(tdmcmc_config,mcmc_param_type):
    ''' Test an error is raised when the input parameter names cannot be found within the file. '''
    with pytest.raises(KeyError) as e_info:
        extract_params(tdmcmc_config,names=['FLIBBLE'],param_type=mcmc_param_type)

def test_extract_param_keep_empty(tdmcmc_config_option,mcmc_param_type):
    ''' Test functionality of exclude_not_found=False option. 
    Check empty values are included for optional parameters when not found in configuration file '''
    optional_param = ['network','unique_copy','emissions_name']
    x = extract_params(tdmcmc_config_option,optional_param=optional_param,exclude_not_found=False,param_type=mcmc_param_type)
    assert optional_param[0] in list(x.keys())
    assert optional_param[1] in list(x.keys())

def test_extract_param_remove_empty(tdmcmc_config_option,mcmc_param_type):
    ''' Test functionality of exclude_not_found=True option.
    Check optional parameters are not included in output when not found in configuration file '''
    optional_param = ['network','unique_copy','emissions_name']
    x = extract_params(tdmcmc_config_option,optional_param=optional_param,exclude_not_found=True,param_type=mcmc_param_type)
    assert optional_param[0] not in list(x.keys())
    assert optional_param[1] not in list(x.keys())

@pytest.mark.basic
def test_all_params(tdmcmc_config_option):
    ''' Test all parameters can be extracted from a file '''
    x = all_param(tdmcmc_config_option)
    assert x

def test_all_params_param_type_1(tdmcmc_config,mcmc_param_type):
    ''' Test all_param function can extract all parameters can be extracted from a file when the parameter type dictionary 
    is specified. '''
    x = all_param(tdmcmc_config,param_type=mcmc_param_type)
    assert x

def test_all_params_param_type_2(example_config,example_param_type):
    ''' Test a different param type dictionary (than mcmc_param_type) to make sure values can be extracted from file. '''
    x = all_param(example_config,param_type=example_param_type)
        
    parameters = []
    types = []
    for section_group in list(example_param_type.values()):
        parameters.extend(list(section_group.keys()))
        types.extend(list(section_group.values()))
    
    for param,t in zip(parameters,types):
        assert isinstance(x[param],t)

@pytest.mark.basic    
def test_all_mcmc_param(tdmcmc_config_option):
    ''' Test all_mcmc_param function can be used. '''
    x = all_mcmc_param(tdmcmc_config_option)
    assert x

def test_all_mcmc_missing(tdmcmc_config_missing):
    ''' Test KeyError is raised if non-optional parameter from template is missing from file '''
    with pytest.raises(KeyError) as e_info:
        all_mcmc_param(tdmcmc_config_missing)

def test_extra_param_all(tdmcmc_config_extra,mcmc_param_type):
    '''
    Test parameters which are within the input file but not in the configuration file will still
    be extracted and the correct type determined.
    '''
    param_names = ["extra_param1","extra_param2"]
    param_types = [int,dict]
    x = all_param(tdmcmc_config_extra,param_type=mcmc_param_type)
    for name,t in zip(param_names,param_types):
        assert isinstance(x[name],t)

def test_extra_param_section(tdmcmc_config_extra,mcmc_param_type):
    '''
    Test parameters which are within the input file but not in the configuration file will still
    be extracted and the correct type determined when a section is specified.
    '''
    param_name = "extra_param1"
    param_type = int
    section = "MEASUREMENTS"
    x = extract_params(tdmcmc_config_extra,param_type=mcmc_param_type,section=section)
    assert isinstance(x[param_name],param_type)

def test_extra_param_section_group(tdmcmc_config_extra,mcmc_param_type):
    '''
    Test parameters which are within the input file but not in the configuration file will still
    be extracted and the correct type determined when a section group is specified.
    '''
    param_name = "extra_param2"
    param_type = dict
    section_group = "MCMC"
    x = extract_params(tdmcmc_config_extra,param_type=mcmc_param_type,section_group=section_group)
    assert isinstance(x[param_name],param_type)
    
def test_extra_mcmc(tdmcmc_config_extra):
    ''' Test all_mcmc_param function can also extract extra parameters. '''
    param_names = ["extra_param1","extra_param2"]
    param_types = [int,dict]

    x = all_mcmc_param(tdmcmc_config_extra)
    for name,t in zip(param_names,param_types):
        assert isinstance(x[name],t)


def config_missing_section_param(param_type,section=None,section_group=None):
    '''
    Find optional parameter names and remaining expected names for optional section otr
    optional_section_group.
    '''
    if section:
        keys = configread.find_param_key(param_type,section=section)[0]
    elif section_group:
        keys = configread.find_param_key(param_type,section_group=section_group)[0]
        
    optional_param = []
    for key in keys:
        optional_param.extend(param_type[key])

    names = configread.all_parameters_in_param_type(param_type)    
    for op in optional_param:
        names.remove(op)
    
    return optional_param,names

def test_missing_section(tdmcmc_config_missing_section,mcmc_param_type):
    '''
    Test optional section can be specified. Also checking all optional parameters
    are not included in output (matching set up o config file)
    '''
    optional_section = "MCMC.MEASUREMENTS"
    
    x = extract_params(tdmcmc_config_missing_section,param_type=mcmc_param_type,
                       optional_section=optional_section,exclude_not_found=True)
    
    optional_param,expected_names=config_missing_section_param(mcmc_param_type,section=optional_section)
    
    for op in optional_param:
        assert op not in x
    
    for name in expected_names:
        assert name in x
    
def test_missing_section_group(tdmcmc_config_missing_section,mcmc_param_type):
    ''' Test optional section group can be specified '''
    optional_section_group = "MCMC"

    x = extract_params(tdmcmc_config_missing_section,param_type=mcmc_param_type,
                       optional_section_group=optional_section_group,exclude_not_found=True)
    
    optional_param,expected_names=config_missing_section_param(mcmc_param_type,section_group=optional_section_group)
    
    for name in expected_names:
        assert name in x
