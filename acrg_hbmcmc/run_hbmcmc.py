'''
Wrapper script to read in parameters from a configuration file and run underlying MCMC script.

Run as:
    $ python run_hbmcmc.py [start end -c config.ini]
e.g.
    $ python run_hbmcmc.py
    $ python run_hbmcmc.py 2012-01-01 2013-01-01 -c hbmcmc_ch4_run.ini

start - Start of date range to use for MCMC inversion (YYYY-MM-DD)
end - End of date range to use for MCMC inversion (YYYY-MM-DD) (must be after start)
-c / --config - configuration file. See config/ folder for templates and examples of this input file.

If start and end are specified these will superceed the values within the configuration file, if present.
If -c option is not specified, this script will look for configuration file within the 
acrg_hbmcmc/ directory called `hbmcmc_input.ini`. 

To generate a config file from the template run this script as:
    $ python run_hbmcmc.py -r  [-c config.ini]

The MCMC run *will not be executed*. This will be named for your -c input or, if not specified, this will 
create a configuration file called `hbmcmc_input.ini` within your acrg_hbmcmc/ directory and exit. 
This file will need to be edited to add parameters for your MCMC run.
'''

import os
import sys
import argparse
from shutil import copyfile
import acrg_hbmcmc.hbmcmc as mcmc
import acrg_config.config as config
from acrg_config.paths import paths
import acrg_hbmcmc.hbmcmc_output as output

def fixed_basis_expected_param():
    '''
    Define required parameters for acrg_hbmcmc.hcmcmc.fixedbasisMCMC() function.
    
    Expected parameters currently include:
        species, sites, meas_period, start_date, end_date, domain, outputpath, outputname
    
    Returns:
        list : required parameter names
    '''
    
    #sg_expected_param = {"INPUT":["species","sites","meas_period","domain","start_date","end_date"],
    #                      "MCMC":["outputpath","outputname"]}
	#
    #expected_param = []
    #for values in sg_expected_param.values():
    #    expected_param.extend(values)

    expected_param = ["species","sites","meas_period","domain","start_date","end_date",
                      "outputpath","outputname"]

    return expected_param

def extract_mcmc_type(config_file,default="fixed_basis"):
    '''
    Find value which describes the MCMC function to use.
    Checks the input configuation file the "mcmc_type" keyword within the "MCMC.TYPE" section.
    If not present, the default is used.
    
    Args:
        config_file (str) :
            Configuration file name. Should be an .ini file.
        default (str) :
            ***
    
    Returns:
        str :
            Keyword for MCMC function to use
    '''
    mcmc_type_section = "MCMC.TYPE"
    mcmc_type_keyword = "mcmc_type"
    param_mcmc_type = config.extract_params(config_file,section=mcmc_type_section)
    
    if param_mcmc_type is not None and mcmc_type_keyword in param_mcmc_type:
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
    function_dict = {"fixed_basis":mcmc.fixedbasisMCMC}
    
    return function_dict[mcmc_type]

def hbmcmc_extract_param(config_file,mcmc_type="fixed_basis",print_param=True,**command_line):
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
        command_line :
            Any additional command line arguments to be added to the param dictionary or to superceed 
            values contained within the config file.
    
    Returns:
        function,collections.OrderedDict :
            MCMC function to use, dictionary of parameter names and values which can be passed to MCMC function
    '''
    
    if mcmc_type == "fixed_basis":
        expected_param = fixed_basis_expected_param()

    # If an expected parameter has been passed from the command line, this does not need to be within the config file
    for key,value in command_line.items():
        if key in expected_param and value is not None:
            expected_param.remove(key)

    mcmc_type_section = "MCMC.TYPE"    
    param = config.extract_params(config_file,expected_param=expected_param,ignore_sections=[mcmc_type_section])

    # Command line values added to param (or superceed inputs from the config file)
    for key,value in command_line.items():
        if value is not None:
            param[key] = value
    
    # If configuration file does not include values for the required parameters - produce an error
    for ep in expected_param:
        if not param[ep]:
            raise Exception(f"Required parameter '{ep}' has not been defined")

    if print_param:
        print("\nInput parameters: ")
        for key,value in param.items():
            print(f"{key} = {value}")

    return param    
    

if __name__=="__main__":

    acrg_path = paths.acrg
    default_config_file = os.path.join(acrg_path,"acrg_hbmcmc/hbmcmc_input.ini")
    config_file = default_config_file

    parser = argparse.ArgumentParser(description="Running Hierarchical Bayesian MCMC script")
    parser.add_argument("start", help="Start date string of the format YYYY-MM-DD",nargs="?")                  
    parser.add_argument("end", help="End date sting of the format YYYY-MM-DD",nargs="?")
    parser.add_argument("-c","--config",help="Name (including path) of configuration file",default=config_file)
    parser.add_argument("-r","--generate",action='store_true',help="Generate template config file and exit (does not run MCMC simulation)")

    args = parser.parse_args()
    
    config_file = args.config or args.config_file
    start_date = args.start
    end_date = args.end

    if args.generate == True:
        template_file = os.path.join(acrg_path,"acrg_hbmcmc/config/hbmcmc_input_template.ini")
        if os.path.exists(config_file):
            write = input(f"Config file {config_file} already exists.\nOverwrite? (y/n): ")
            if write.lower() == "y" or write.lower() == "yes":        
                copyfile(template_file,config_file)
            else:
                sys.exit(f"Previous configuration file has not been overwritten.")
        else:
            copyfile(template_file,config_file)
        sys.exit(f"New configuration file has been generated: {config_file}")

    if not os.path.exists(config_file):
        if config_file == default_config_file:
            sys.exit("No configuration file detected.\nTo generate a template configuration file run again with -r flag:\n  $ python run_tdmcmc.py -r")
        else:
            sys.exit(f"Configuration file cannot be found.\nPlease check path and filename are correct: {config_file}")

    mcmc_type = extract_mcmc_type(config_file)
    mcmc_function = define_mcmc_function(mcmc_type)
    print(f"Using MCMC type: {mcmc_type} - function {mcmc_function.__name__}(...)")
    
    param = hbmcmc_extract_param(config_file,mcmc_type,start_date=start_date,end_date=end_date)

    output.copy_config_file(config_file,param=param,start_date=start_date,end_date=end_date)

    mcmc_function(**param)

