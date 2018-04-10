#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 16:39:01 2017

This module allows configuration files in the INI format to be read and used.

Example of a section in the configuration file:
------------------------------------------------
[MEASUREMENTS]
# Measurement details

sites = ["GSN"]           ; Sites to read the data from as a list
species = "chcl3"
start_date = "2015-01-01" ; Default start date used if none specified on the command line
end_date = "2015-02-01"   ; Default start date used if none specified on the command line
domain = "EASTASIA"
network = "AGAGE"
------------------------------------------------

 - Sections are included in square brackets
 - Parameter name and value pairs are separated by an equals sign
     - Values can be specified with the same syntax as when creating a python object e.g. '' for string, [] for lists (and also for np.array - will be converted if their type has been set as an array)
 - ; and # symbols can be used to create new line and inline comments

Section headings can be of the form [NAME] or [GROUP.NAME]. This allows paramaters to be separated into several
section headings in the configuration file for clarity but grouped into one overall classification when inputted based on the GROUP
name.

param_type dictionary
+++++++++++++++++++++

To specify inputs and the types they should be cast to (e.g. str, array, boolean etc) a nested dictionary should be created.
This can be passed into several functions as the param_type argument.

This should be of the form of one of the following:
    {'SECTION_GROUP1':{'param1':str,'param2':float},'SECTION_GROUP2':{'param3':list,'param4':np.array}}
    {'SECTION1':{'param1':str},'SECTION2':'{param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}}
    OrderedDict(['SECTION_GROUP1':OrderedDict([('param1':str),('param2':float)]),'SECTION_GROUP2':OrderedDict([('param3':list),('param4':np.array)]))
    OrderedDict(['SECTION1':{'param1':str},'SECTION2':{'param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}])

This can either be created directly or a template configuration file can be created and a param_type dictionary created 
from this using the generate_param_dict() function.
These template files should be kept within the acrg_config/templates/ directory and, after creation, should not be altered
unless you wish to change the format for all config files of this type.

Note: if param_type is not defined, the code will attempt to cast the inputs to the most sensible type.
This should be fine in most cases but could cause issues.
This also means the inputs will not be checked for any missing values.


How to run
++++++++++

The main functions to use for reading in parameters from a config file are:
    
    * all_param(config_file,...)      - Extract all parameters from a configuration file.
    * extract_params(config_file,...) - Extract specific parameters from a file either based on parameter 
    names, sections or groups.

A param_type dictionary can be defined both to fix expected inputs and to explictly specify the parameter types.

@author: rt17603
"""

import configparser
from collections import OrderedDict
import numpy as np

def open_config(config_file):
    '''
    The open_config function is used to open configuration files in the ini format.
    
    Args:
        config_file : filename for input configuration file (str).
    
    Returns
        configparser.ConfigParser object
    '''
    config = configparser.ConfigParser(inline_comment_prefixes=(';','#'))
    config.optionxform=str # Keeps case when inputting option names
    
    with open(config_file) as fp:
        config.readfp(fp)
    
    return config

def generate_param_dict(config_file):
    '''
    The generate_param_dict function creates a param_type nested dictionary from an input configuration file.
    This could be used on some fixed example config file to generate the parameter type dictionary which can then be applied
    to other configuration files of the same type.
    
    Args:
        config_file : filename for input configuration file (str).
    
    Returns:
        nested OrderedDict : parameter type dictionary from configuration file input
    '''
    config = open_config(config_file)
    
    sections = config.sections() # Extract all section names from the config file
    
    param_type = OrderedDict([])
    
    for section in sections:
        section_param = config[section].keys()
        if section_param:
            print section_param
            types = [type(convert(value)) for value in config[section].values()]
            print types
            section_types = OrderedDict([(key,value) for key,value in zip(section_param,types)])
            param_type[section] = section_types
    
    return param_type


def str_check(string,error=True):
    '''
    The str_check function is used as part of checking the input from a configuration file.
    This function ensures the input remains as a string and removes any " or ' characters
    
    Args:
        string : value input from config file
    
    Returns
        string (formatted)
    '''
    
    # Remove any ' or " symbols surrounding the input string from the config_file

    string = string.strip() # Strip any whitespace just in case
    if (string[0] == "'" and string[-1] == "'"):
        string = string[1:-1]
    elif (string[0] == '"' and string[-1] == '"'):
        string = string[1:-1]
    
    try:
        out = str(string)
        out.encode('ascii','ignore')
    except (TypeError,SyntaxError):
        if error:
            print "Could not convert input parameter '{0}' to str.".format(out)
        return None
    
    return out


def eval_check(string,error=True):
    '''
    The eval_check function evaluates the input string from a configuration file to a python object.
    For example: '1' would evalute to an int object 1
               : '1.' or '1.0' would evaluate to float object 1.0
               : "[1,2,3]" would evalute to a list object [1,2,3]
               : "1,2,3" would evaluate to a tuple object (1,2,3)
               See eval() documentation for full list
    Args:
        string : value input from config file
        error  : return error message if unable to evaluate (bool).
    
    Returns:
        Python object
        
        If unable to convert string to python object:
            None
    '''
    try:
        out = eval(string)
    except (NameError,SyntaxError):
        try:
            out = eval("'"+string+"'") # An input string without quotes cannot be evaluated so try adding quotes
        except (NameError,SyntaxError):
            if error:
                print "Could evaluate input '{0}' to any type.".format(string)
            return None
    
    return out


def list_check(string,force_convert=True,error=True):
    '''
    The list_check function converts input string to a list.
    
    Args:
        string (str)         : value input from config file
        force_convert (bool) : whether the conversion should be forced by creating a list object
        error (bool)         : return error message if unable to convert.
    
    Returns:
        List
        
        If unable to convert to list:
            None    
    '''
    
    out = eval_check(string) # Try evaluating input
    
    if not isinstance(out, (list, str)): # If not already a list
        try:
            out = list(out)
        except TypeError:
            if force_convert:
                out = [out]
            else:
                if error:
                    print "Could not convert input parameter '{0}' to list.".format(out)
                return None
    elif isinstance(out, str):
        out = [out]
    
    return out


#def float_check(string,error=True):
#    
#    out = eval_check(string)
#    
#    try:
#        out = float(out)
#    except ValueError:
#        if error:
#            print 'Could not convert input parameter to float: {0}'.format(out)
#        return None
#    
#    return out    

#def int_check(string,error=True):
#    
#    out = eval_check(string)
#    
#    try:
#        out = int(out)
#    except ValueError:
#        if error:
#            print 'Could not convert input parameter to int: {0}'.format(out)
#        return None
#    
#    return out
    

def convert(string,value_type=None):
    '''
    The convert function converts the input string to the specified value_type.
    If no value_type is stated, the function attempts to discern the appropriate type.
    
    Args:
        string     : value input from config file
        value_type : object type. Values accepted: str,list
    
    Returns:
        Python object of specified type
        
        If unable to convert to value_type:
            None
        If unable to find a suitable type when value_type is not specified:
            None
    '''
    
    if value_type is str:
        out = str_check(string)
    elif value_type is list:
        out = list_check(string)
    elif value_type is np.ndarray:
        out = list_check(string)
        out = np.array(out)
#    elif value_type is float:
#        out = float_check(string)
#    elif value_type is int:
#        out = int_check(string)
    else:
        out = eval_check(string,error=False)
        if out is None:
            out = str_check(out,error=False)
    
    return out       


#def check_params(param,param_type,keys=[],optional_param=[],raise_exception=True):
#    '''
#    The check_params function
#    '''
#    names = param.keys()
#    
#    if keys:
#        check_keys = keys
#    else:
#        check_keys = param_type.keys()
#    
#    check_names = []
#    for key in check_keys:
#        try:
#            check_names.extend(param_type[key].keys())
#        except KeyError:
#            print "Key '{0}' cannot be found in param_type passed to check_params function".format(key)
#            return None
#    
#    for name in check_names:
#        if name not in names:
#            if name not in optional_param:
#                if raise_exception:
#                    raise Exception("Parameter '{0}' must be specified in configuration file".format(name))
#                else:
#                    return False
#    
#    return True

def all_parameters_in_param_type(param_type):
    '''
    The all_parameters_in_param_type function extracts all parameters (regardless of section/section_group) for a given 
    param_type nested dictionary.
    
    Args:
        param_type : nested dictionary of expected parameter names and types.
                     Key for each parameter dictionary can be the section heading or the overall group (e.g. for [MCMC.MEASUREMENTS], section group should be 'MCMC').
                     If param_type is explictly specified should be of form of one of the following:
                          {'SECTION_GROUP1':{'param1':str,'param2':float},'SECTION_GROUP2':{'param3':list,'param4':np.array}}
                          {'SECTION1':{'param1':str},'SECTION2':'{param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}}
                           OrderedDict(['SECTION_GROUP1':OrderedDict([('param1':str),('param2':float)]),'SECTION_GROUP2':OrderedDict([('param3':list),('param4':np.array)]))
                           OrderedDict(['SECTION1':{'param1':str},'SECTION2':{'param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}])

    Returns:
        list : list of all parameters in param_type dictionary
    '''
    types = param_type # Get dictionary containing parameter names and types
    keys = types.keys() # Extract all section/section_group headings
    parameters = []
    for key in keys:
        parameters.extend(types[key].keys()) # For each section/section_group extract all the associated parameters
    
    return parameters

def find_param_key(param_type,section=None,section_group=None):
    '''
    The find_param_key function checks whether the keys within the param_type dictionary are for sections (e.g. 'MCMC.MEASUREMENTS') 
    or groups (e.g. 'MCMC') and returns the relevant key.
    One of section or section_group should be specified.
    Returned key_type is one of 'section' or 'section_group'.
    
    Args:
        param_type          : nested dictionary of expected parameter names and types.
                              Key for each parameter dictionary can be the section heading or the overall group (e.g. for [MCMC.MEASUREMENTS], section group should be 'MCMC').
                              If param_type is explicly specified should be of form of one of the following:
                                  {'SECTION_GROUP1':{'param1':str,'param2':float},'SECTION_GROUP2':{'param3':list,'param4':np.array}}
                                  {'SECTION1':{'param1':str},'SECTION2':'{param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}}
                                  OrderedDict(['SECTION_GROUP1':OrderedDict([('param1':str),('param2':float)]),'SECTION_GROUP2':OrderedDict([('param3':list),('param4':np.array)]))
                                  OrderedDict(['SECTION1':{'param1':str},'SECTION2':{'param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}])
        section (str)       : name of section in config file for the parameter
        section_group (str) : name of group in config file for the parameter
    
    Returns:
        str,str: key, key_type                      
    
    '''
    types = param_type # Get dictionary containing parameter names and types
    keys = types.keys() # Extract all section/classification keys
    
    # Find oarameter class if not specified (should be defined as first part of section split by '.' e.g. MCMC.MEASUREMENTS, section_group='MCMC')
    if not section_group:
        section_group = section.split('.')[0]
        
    if section in keys:
        key = section
        key_type = 'section'
    elif section_group in keys:
        key = section_group
        key_type = 'section_group'
    elif section_group in [k.split('.')[0] for k in keys]:
        keys_starter = [k.split('.')[0] for k in keys]
        key = keys[keys_starter.index(section_group)]
        key_type = 'section'
    else:
        key = None
        key_type = None
        #raise Exception('Section/Classification {0}/{1} does not match to any key in input param_type'.format(section_group,section))
        #print 'Param class cannot be found i for section of parameters not defined. Using {0} as default'.format(section_groups[0])
        #section_group = section_groups[0]
    
    return key,key_type
    


def get_value(name,config,section,param_type=None):
    '''
    The get_value function extracts the value of a parameter from the configuration file.
    This value is then converted to the type specified within param_type (default from mcmc_param_type() function).
    
    Args:
        name (str)            : name of the parameter to extract
        config (ConfigParser) : ConfigParser object created using the configparser class. Data from config file should have been read in
        section (str)         : name of section in config file for the parameter
        param_type (dict)     : nested dictionary of parameter classes and expected parameter names and types.
                                Key for each parameter dictionary can be the section heading or the overall classification (e.g. for [MCMC.MEASUREMENTS], classification should be 'MCMC').
                                If param_type is explicly specified should be of form:
                                    {'SECTION_GROUP1':{'param1':str,'param2':float},'SECTION_GROUP2':{'param3':list,'param4':np.array}}
                                    {'SECTION1':{'param1':str},'SECTION2':'{param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}}
                                    OrderedDict(['SECTION_GROUP1':OrderedDict([('param1':str),('param2':float)]),'SECTION_GROUP2':OrderedDict([('param3':list),('param4':np.array)]))
                                    OrderedDict(['SECTION1':{'param1':str},'SECTION2':{'param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}])
                                
    Returns:
        value
        
        If param_type is specified
            If neither section nor section_group can be identified within param_type dictionary:
                Exception raised and program exited
            If parameter name cannot be identified within param_type dictionary:
                Exception raised and program exited
    '''
    
    if param_type:
        key,key_type = find_param_key(param_type,section)
        types = param_type
        try:
            value_type = types[key][name] # Find specified type of object for input parameter
        except KeyError:
            #raise Exception('Input parameter {0} in section {1} not expected (not found in param_type dictionary [{2}][{0}])'.format(name,section,key))
            print "Type for input name '{0}' is not specified.".format(name)
            value_type = None
    else:
        #print "Type for input name '{0}' is not specified.".format(name)
        #value_type = str
        value_type = None
    
    # For int, float and bool object functions within config module exist to cast directly to these types
    if value_type == int:
        value = config.getint(section,option=name)
    elif value_type == float:
        value = config.getfloat(section,option=name)
    elif value_type == bool:
        value = config.getboolean(section,option=name)
    else: # For lists, np.ndarrays and str objects we extract as a string and manually format
        value = config.get(section,option=name)
        value = convert(value,value_type=value_type)
    
    return value


def extract_params(config_file,section=None,section_group=None,names=[],optional_param=[],remove_not_found=False,
                   param_type=None):
    '''
    The extract_params function extracts parameter names and values from a configuration file.
    The parameters which are extracted is dependent on whether the section, section_group and/or names variables are specified.
    A param_type dictionary can be defined to ensure variables are cast to the correct types.
    
    Args:
        config_file             : filename for input configuration file.
        section (str)           : extract parameters from section name.
        section_group (str)     : extract parameters from all sections with this group.
                                  If section and section_group are both specified - section takes precedence.
        names (list)            : which parameter names to extract (within section or section_group)
        optional_param (list)   : parameters which are optional. If the param cannot be found value will be set to None
        param_type (dict)       : nested dictionary of sections or groups and expected parameter names and types.
                                  If param_type is explicly specified should be of form:
                                      {'SECTION_GROUP1':{'param1':str,'param2':float},'SECTION_GROUP2':{'param3':list,'param4':np.array}}
                                      {'SECTION1':{'param1':str},'SECTION2':'{param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}}
                                      OrderedDict(('SECTION_GROUP1':OrderedDict(('param1':str),('param2':float)),'SECTION_GROUP2':OrderedDict(('param3':list),('param4':np.array)))
                                      OrderedDict(('SECTION1':{'param1':str},'SECTION2':{'param2':float}'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}))        
        remove_not_found (bool) : whether to remove parameters which are not found in the input file or include them as None.
        
    Returns:
        OrderedDict : parameter names and values
        
        If parameter cannot be found (not specified as an optional param)
            Exception raised and program exited
    '''
    
    # Open config file with configparser
    config = open_config(config_file)
    
    all_sections = config.sections() # Extract all section names from the config file
    
    if section:
        if section in all_sections:
            select_sections = [section] # Only look within selected section
        else:
            raise KeyError('Specified section {0} could not be found in configuration file: {1}'.format(section,config_file[0]))
            #select_sections = []
    elif section_group:
        sections = all_sections
        select_sections = [s for s in sections if s.split('.')[0].lower() == section_group.lower()] # Find all sections covered by section_group (section_group.name)
        if not select_sections:
            raise KeyError('No sections could be found for specified section_group {0} in configuration file: {1}'.format(section_group,config_file[0]))
    else:
        select_sections = all_sections # Find all sections
    
    # Finding all parameter names in the sections we want to extract
    extracted_names = []
    match_section = []
    for sect in select_sections:
        k = config[sect].keys()
        s = [sect]*len(k)
        extracted_names.extend(k) # List of the parameter names within the input file, within specified sections
        match_section.extend(s)   # Associated list with the section heading for each parameter
       
    if not names:
        if param_type:
            if section_group:
                key_value,key_type = find_param_key(section_group=section_group,param_type=param_type)
                keys = [key_value]
            elif section:
                key_value,key_type = find_param_key(section=section,param_type=param_type)
                keys = [key_value]
            else:
                keys = param_type.keys()
            #print 'Keys to extract input names from param_type: {0}'.format(keys)
            names = []
            if (section_group and key_type == 'section_group') or (section and key_type == 'section') or (section_group and key_type == 'section'):
               for key in keys:
                   names.extend(param_type[key].keys())
            elif (section and key_type == 'section_group'):
                print 'Cannot create list of necessary input parameters created based on param_type input. Please check all inputs are included manually.'
                names = extracted_names # Set to just match names extracted from file           
            else:
                names = all_parameters_in_param_type(param_type) # Extract all parameter names from param_type dictionary
        else:
            names = extracted_names # Set to just match names extracted from file
    
    #print 'Parameter names to extract from configuration file: {0}'.format(names)
    #print 'Relevant parameter names within configuration file: {0}'.format(extracted_names)
    
    param = OrderedDict({})
    
    #if param_type:
    for name in names:
        if name in extracted_names:
            try:
                index = extracted_names.index(name)
            except ValueError:
                print "Parameter '{0}' not found in configuration file (check specified section {1} or section_group {2} is correct).".format(name,section,section_group)
            else:
                param[name] = get_value(name,config,match_section[index],param_type)
        else:
            if name in optional_param:
                param[name] = None
            elif not param_type:
                print "Parameter '{0}' not found in configuration file (check specified section {1} or section_group {2} is correct).".format(name,section,section_group)
            else:
                if section:
                    raise KeyError("Parameter '{0}' not found in input configuration file in section '{1}'".format(name,section))
                elif section_group:
                    raise KeyError("Parameter '{0}' not found in input configuration file within section_group '{1}'".format(name,section_group))
                else:
                    raise KeyError("Parameter '{0}' not found in input configuration file.".format(name))
    #else:
    #    for name in names:
    #        try:
    #            index = extracted_names.index(name)
    #        except ValueError:
    #            print "Parameter '{0}' not found in configuration file (check specified section {1} or section_group {2} is correct).".format(name,section,section_group)
    #        else:
    #            section = match_section[index]
    #            param[name] = get_value(name,config,section,param_type)
    
    if remove_not_found:
        param = OrderedDict([(key,value) for key,value in param.iteritems() if value != None])
    
    return param

def all_param(config_file,optional_param=[],param_type=None,remove_not_found=False):
    '''
    The all_param function extracts all parameters from a config file.
    If param_type specified will cast to the specified types, otherwise will attempt to discern the parameter types from the form of the values.
    
    Args:
        config_file             : filename for input configuration file
        optional_param          : parameters which are optional. If the param cannot be found in input file, value will be set to None
        param_type              : nested dictionary of sections or groups and expected parameter names and types.
                                  If param_type is specified should be of form:
                                   {'SECTION_GROUP1':{'param1':str,'param2':float},'SECTION_GROUP2':{'param3':list,'param4':np.array}}
                                   {'SECTION1':{'param1':str},'SECTION2':'{param2':float},'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}}
                                    OrderedDict(('SECTION_GROUP1':OrderedDict(('param1':str),('param2':float)),'SECTION_GROUP2':OrderedDict(('param3':list),('param4':np.array)))
                                    OrderedDict(('SECTION1':{'param1':str},'SECTION2':{'param2':float}'SECTION3':{'param3':list},'SECTION4':{'param4':np.array}))        
        remove_not_found (bool) : whether to remove parameters which are not found in the input file or include them as None.    
    Returns:
        OrderedDict: parameter names and values 
    '''
    
    param = OrderedDict({})
    param = extract_params(config_file,optional_param=optional_param,param_type=param_type,remove_not_found=remove_not_found)
    
    return param


    