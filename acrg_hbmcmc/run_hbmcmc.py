'''
Wrapper script to read in parameters from a configuration file and run the
acrg_hbmcmc.hbmcmc.fixedbasisMCMC() function.
'''

import os
import argparse
#import acrg_hbmcmc.hbmcmc as mcmc
import acrg_config.config as config
from acrg_config.paths import paths

def fixed_basis_expected_param(section_group=None):
    '''
    Define required parameters for acrg_hbmcmc.hcmcmc.fixedbasisMCMC() function.
    
    Expected parameters currently include:
        species, sites, meas_period, start_date, end_date, domain, outputpath, outputname
    
    Args:
        section_group (str/None, optional) :
            Select expected parameters for a particular section group. At the moment this includes:
                "INPUT"
                "MCMC"
            None - returns expected parameters for all groups
    
    Returns:
        list : required parameter names
    '''
    sg_expected_param = {"INPUT":["species","sites","meas_period","start_date","end_date","domain"],
                          "MCMC":["outputpath","outputname"]}
	
    if section_group is None:
        section_group = list(sg_expected_param.keys())

    expected_param = []
    for sg in section_group:
        expected_param.extend(sg_expected_param[sg])

    return expected_param

def extract_mcmc_type(config_file,default="fixed_basis"):
    '''
    Find value which describes the MCMC function to use.
    Checks the input configuation file the "mcmc_type" keyword within the "MCMC.TYPE" section.
    If not present, the default is used.
    
    Args:
        config_file (str) :
            Configuration file name. Should be an .ini file.
    
    Returns:
        str :
            Keyword for MCMC function to use
    '''
    mcmc_type_section = "MCMC.TYPE"
    mcmc_type_keyword = "mcmc_type"
    param_mcmc_type = config.extract_params(config_file,section=mcmc_type_section)
    if mcmc_type_keyword in param_mcmc_type:
        mcmc_type = param_mcmc_type[mcmc_type_keyword]
    else:
        mcmc_type = default
       
    return mcmc_type
    
def define_mcmc_function(mcmc_type):
    '''
    Links mcmc_type name to function.
    
    Current options:
        mcmc_type (str) :
            "fixed_basis" : acrg_hbmcmc.hbmcmc.fixedbasisMCMC(...) function
    
    Returns:
        Function
    '''
    function_dict = {"fixed_basis":"placeholder function name"}#mcmc.fixedbasisMCMC}
    
    return function_dict[mcmc_type]

def hbmcmc_extract_param(config_file,mcmc_type="fixed_basis",print_param=True):
    '''
    Extract parameters from input configuration file and associated MCMC function. 
    Checks the mcmc_type to extract the required parameters.
    
    Args:
        config_file (str):
            Configuration file name. Should be an .ini file.
        mcmc_type (str, optional) :
            Keyword for MCMC function to use.
            Default = "fixed_basis" (only option at present)
        print_param (bool, optional) :
            Print out extracted parameter names.
            Default = True
    
    Returns:
        function,collections.OrderedDict :
            MCMC function to use, dictionary of parameter names and values which can be passed to MCMC function
    '''
    
    if mcmc_type == "fixed_basis":
        expected_param = fixed_basis_expected_param()

    mcmc_type_section = "MCMC.TYPE"    
    param = config.extract_params(config_file,expected_param=expected_param,ignore_sections=[mcmc_type_section])

    for ep in expected_param:
        if not param[ep]:
            raise Exception(f"Expected parameter {param[ep]} has not been defined")

    if print_param:
        print("Input parameters: ")
        for key,value in param.items():
            print(f"{key} = {value}")

    return param    
    

if __name__=="__main__":

    acrg_path = paths.acrg
    config_file = os.path.join(acrg_path,"acrg_hbmcmc/hbmcmc_input.ini")    

    parser = argparse.ArgumentParser(description="Running hbmcmc script")
    parser.add_argument("-c","--config",help="Configuration filename",default=config_file)

    args = parser.parse_args()
    config_file = args.config or args.config_file

    mcmc_type = extract_mcmc_type(config_file)
    mcmc_function = define_mcmc_function(mcmc_type)
    print(f"Using MCMC type: {mcmc_type} - function {mcmc_function}(...)")
    
    param = hbmcmc_extract_param(config_file,mcmc_type)

    #mcmc.fixedbasisMCMC(**param)
    mcmc_function(**param)

