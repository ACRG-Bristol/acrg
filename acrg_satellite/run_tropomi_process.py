#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 17:33:47 2020

@author: rt17603
"""

#import acrg_satellite.gosat as gosat
#import acrg_satellite.gosat_config as gosat_config
#import os
#import sys
#import argparse

import os

import acrg_satellite.tropomi as tropomi
from pathlib import Path
import argparse

try:
    from acrg_config.paths import paths
except ImportError:
    acrg_path = os.getenv("ACRG_PATH")
else:
    acrg_path = paths.acrg



if __name__=="__main__":
    
    config_file = os.path.join(acrg_path,"acrg_satellite/tropomi_process.ini")
    
    parser = argparse.ArgumentParser(description='Running tropomi process script')

    parser.add_argument("start", help="Start date string yyyy-mm-dd",type=str,nargs="?")
    parser.add_argument("end", help="End date sting yyyy-mm-dd",type=str,nargs="?")
    parser.add_argument("-c","--config",help="Configuration filename",default=config_file)

    args = parser.parse_args()
    
    start_date = args.start
    end_date = args.end
    config_file = args.config or args.config_file
    
    # tropomi_param = gosat_config.gosat_param(config_file)

#    if start_date is not None:
#      tropomi_param["start_date"] = start_date
#    if end_date is not None:
#      tropomi_param["end_date"] = end_date


    ####### TO BE REMOVED ONCE CONFIG FILE IS IN PLACE #####
    start_date = "2019-07-19"
    end_date = "2019-07-22"
    
    area="GOSAT-BRAZIL"

    if area == "GOSAT-SOUTHAMERICA":
        ## Match to GOSAT-SOUTHAMERICA bounds
        lat_min = -57.0 # -56.98.071
        lat_max = +13.5 # +13.45861
        lon_min = -82.0 # -82.32861
        lon_max = -32.7 # -32.79305
    elif area == "GOSAT-BRAZIL":
        ## Match to GOSAT-BRAZIL bounds
        lat_min = -35.8 # -35.7525
        lat_max = +7.3  # +7.25138
        lon_min = -76.0 # -75.98444
        lon_max = -32.7 # -32.79305
    
    dlat = 0.234
    dlon = 0.352
#    dlat = 0.10
#    dlon = 0.10

    lat_bounds = [lat_min,lat_max]
    lon_bounds = [lon_min,lon_max]
    coord_bin = [dlat,dlon]
    
#    time_increment = "1min" # e.g. None, "10s", "30min" (always < 1 day)
    #time_increment = "30min" # e.g. None, "10s", "30min" (always < 1 day)
    time_increment = "1h" # e.g. None, "10s", "30min" (always < 1 day)
    #time_increment = None
    
    site = area.replace("GOSAT-","TROPOMI-")
    network = "TROPOMI/"+site
    
    write_name = True
    name_directory = Path("/home/rt17603/work/NAME_files/test")
    #name_directory = Path("/home/rt17603/work/NAME_files/")
#    name_directory = Path("/home/rt17603/work/NAME_files/03_bin")
#    name_directory = Path("/home/rt17603/work/NAME_files/05_bin")
#    name_directory = Path("/home/rt17603/work/NAME_files/035_1min_bin")
    #name_directory = Path("/home/rt17603/work/NAME_files/010_1min_bin")
    use_name_pressure = True
    pressure_domain = "SOUTHAMERICA"
    pressure_max_days = 1
    
    output_directory = Path("/home/rt17603/work/obs/test")
#    output_directory = Path("/home/rt17603/work/obs/")
    #output_directory = Path("/home/rt17603/work/obs/03_bin")
#    output_directory = Path("/home/rt17603/work/obs/05_bin")
#    output_directory = Path("/home/rt17603/work/obs/035_1min_bin")
    #output_directory = Path("/home/rt17603/work/obs/010_1min_bin")
    
    #max_points = None # Maximum number of points per NAME file
    max_points = 50 # Maximum number of points per NAME file
    
    #regrid_method = "conservative_normed"
    regrid_method = "conservative"

    
    tropomi.tropomi_process(start_date,end_date,lat_bounds,lon_bounds,coord_bin,
                            site=site,network=network,
                            time_increment=time_increment,write_name=write_name,
                            regrid_method=regrid_method,
                            max_points=max_points,
                            use_name_pressure=use_name_pressure,
                            pressure_domain=pressure_domain,
                            pressure_max_days=pressure_max_days,
                            name_directory=name_directory,
                            output_directory=output_directory)
    
    ################################
    
    # #print("Input parameters for tropomi function: ",tropomi_param)
    #tropomi.tropomi_process(**tropomi_param)
