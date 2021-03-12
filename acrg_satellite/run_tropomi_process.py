#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 17:33:47 2020

@author: rt17603
"""

import os
import sys
from pathlib import Path
import argparse
from acrg_satellite import tropomi,tropomi_config
from acrg_config.config import generate_from_template

try:
    from acrg_config.paths import paths
except ImportError:
    acrg_path = os.getenv("ACRG_PATH")
else:
    acrg_path = paths.acrg


if __name__=="__main__":
    
    default_config_file = os.path.join(acrg_path,"acrg_satellite/tropomi_process.ini")
    config_file = default_config_file
    
    parser = argparse.ArgumentParser(description='Running tropomi process script')

    parser.add_argument("start", help="Start date string yyyy-mm-dd",type=str,nargs="?")
    parser.add_argument("end", help="End date sting yyyy-mm-dd",type=str,nargs="?")
    parser.add_argument("-c","--config",help="Configuration filename",default=config_file)
    parser.add_argument("-r","--generate",action='store_true',help="Generate template config file and exit (does not run tropomi process)")

    args = parser.parse_args()
    
    config_file = args.config or args.config_file

    if args.generate is True:
        template_file = os.path.join(acrg_path,"acrg_config/templates/tropomi_process_template.ini")
        if os.path.exists(config_file):
            write = input(f"Config file {config_file} already exists.\nOverwrite? (y/n): ")
            if write.lower() == "y" or write.lower() == "yes":        
                generate_from_template(template_file,config_file)
            else:
                sys.exit(f"Previous configuration file has not been overwritten.")
        else:
            generate_from_template(template_file,config_file)
        sys.exit(f"New configuration file has been generated: {config_file}")

    if not os.path.exists(config_file):
        if config_file == default_config_file:
            sys.exit("No configuration file detected.\nTo generate a template configuration file run again with -r flag:\n  $ python run_tropomi_process.py -r")
        else:
            sys.exit(f"Configuration file cannot be found.\nPlease check path and filename are correct: {config_file}")

    command_line_args = {}
    if args.start:
        command_line_args["start_date"] = args.start
    if args.end:
        command_line_args["end_date"] = args.end

    
    tropomi_param = tropomi_config.tropomi_param(config_file,**command_line_args)
    
    print("Input parameters for tropomi function: ", tropomi_param)
    tropomi.tropomi_process(**tropomi_param)

