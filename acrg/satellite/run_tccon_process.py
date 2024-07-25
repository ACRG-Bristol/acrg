#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 14:25:59 2018

The Run_TCCON_process script takes a tccon configuration file and uses this to run the 
acrg.satellite.tccon_process(...) function.

The format of the configuration file should match acrg_config/templates/tccon_process_template.ini 
although most parameters are optional. If any input parameters for tccon_process(...) are not 
specified in the input configuration file the defaults for the function will be used.

* Please check the defaults within the documentation for the tccon_process function to ensure 
this is what is required *

Note: By default nothing will be written to file, only a merged xarray.Dataset is produced,
 and so this script will produce no output. Please set write_nc and/or write_name inputs to
 True to produce and output and specify output directories for these.

See acrg_satellite.tccon module for further details.

Run as:
    >> python run_tccon_process.py -c tccon_param.ini

@author: rt17603
"""

import os
import argparse

import acrg.satellite.tccon as tccon
import acrg.satellite.tccon_config as tccon_config

if __name__=="__main__":

    from acrg.config.paths import Paths
    acrg_path = Paths.acrg

    config_file = os.path.join(acrg_path,"acrg_satellite/tccon_process.ini")
    
    parser = argparse.ArgumentParser(description='Running tccon process script')
    parser.add_argument("-c","--config",help="Configuration filename",default=config_file)
    
    args = parser.parse_args()
    config_file = args.config or args.config_file
    
    tccon_param = tccon_config.tccon_param(config_file)
    #print("Input parameters for tccon function: ",tccon_param)
    tccon_ds = tccon.tccon_process(**tccon_param)
