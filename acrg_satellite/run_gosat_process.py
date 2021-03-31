#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 14:25:59 2018

The Run_GOSAT_process script takes a gosat configuration file and uses this to run the 
acrg_gosat.gosat_process(...) function.

The format of the configuration file should match acrg_config/templates/gosat_process_template.ini 
although most parameters are optional. If any input parameters for gosat_process(...) are not 
specified in the input configuration file the defaults for the function will be used.

* Please check the defaults within the documentation for the gosat_process function to ensure 
this is what is required *

Note: By default nothing will be written to file, only a merged xarray.Dataset is produced,
 and so this script will produce no output. Please set write_nc and/or write_name inputs to
 True to produce and output and specify output directories for these.

See acrg_satellite.gosat module for further details.

Run as:
    >> python run_gosat_process.py -c gosat_param.ini

@author: rt17603
"""

#import gosat
#import gosat_config
import acrg_satellite.gosat as gosat
import acrg_satellite.gosat_config as gosat_config
import os
import sys
import argparse

if __name__=="__main__":

    from acrg_config.paths import paths
    acrg_path = paths.acrg

    config_file = os.path.join(acrg_path,"acrg_satellite/gosat_process.ini")
    
    parser = argparse.ArgumentParser(description='Running gosat process script')
    parser.add_argument("-c","--config",help="Configuration filename",default=config_file)
    
    args = parser.parse_args()
    config_file = args.config or args.config_file
    
    gosat_param = gosat_config.gosat_param(config_file)
    #print("Input parameters for gosat function: ",gosat_param)
    gosat_ds = gosat.gosat_process(**gosat_param)
