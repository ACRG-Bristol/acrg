'''
Wrapper script to read in parameters from a configuration file and run the
acrg_hbmcmc.hbmcmc.fixedbasisMCMC() function.
'''

import os
import argparse
import acrg_hbmcmc.hbmcmc as mcmc
import acrg_config.config as config
from acrg_config.paths import paths

def hbmcmc_expected_param(section_group=None):

	sg_expected_param = {"INPUT":["species","sites","meas_period","start_date","end_date","domain"],"MCMC":["outputpath","outputname"]}
	
	if section_group is None:
		section_group = list(sg_expected_param.keys())

	expected_param = []
	for sg in section_group:
		expected_param.extend(sg_expected_param[sg])

	return expected_param


if __name__=="__main__":

    acrg_path = paths.acrg
    config_file = os.path.join(acrg_path,"acrg_hbmcmc/hbmcmc_input.ini")    

    parser = argparse.ArgumentParser(description="Running hbmcmc script")
    parser.add_argument("-c","--config",help="Configuration filename",default=config_file)

    args = parser.parse_args()
    config_file = args.config or args.config_file

    expected_param = hbmcmc_expected_param()
    param = config.extract_params(config_file,expected_param=expected_param)
    
	for ep in expected_param:
		if not param[ep]:
			raise Exception(f"Expected parameter {param[ep]} has not been defined")

    print("Input parameters: ")
    for key,value in param.items():
        print(f"{key} = {value}")

    mcmc.fixedbasisMCMC(**param)

