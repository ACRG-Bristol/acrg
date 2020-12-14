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

import numpy as np
import acrg_satellite.tropomi as tropomi
from pathlib import Path

if __name__=="__main__":
    
    start_date = "2019-07-11"
    end_date = "2019-07-15"
    
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
    
    dlat = 0.1
    dlon = 0.1
    
    lat_bounds = [lat_min,lat_max]
    lon_bounds = [lon_min,lon_max]
    coord_bin = [dlat,dlon]
    
    time_increment = "1s" # e.g. None, "10s", "30min" (always < 1 day)
    #time_increment = "10s" # e.g. None, "10s", "30min" (always < 1 day)
    #time_increment = "30min" # e.g. None, "10s", "30min" (always < 1 day)
    #time_increment = "1h" # e.g. None, "10s", "30min" (always < 1 day)
    #time_increment = None
    
    site = area.replace("GOSAT-","TROPOMI-")
    network = "TROPOMI/"+site
    
    write_name = True
    #name_directory = Path("/home/rt17603/work/NAME_files/test")
    name_directory = Path("/home/rt17603/work/NAME_files/")
    use_name_pressure = True
    pressure_domain = "SOUTHAMERICA"
    pressure_max_days = 1
    
    #output_directory = Path("/home/rt17603/work/obs/test")
    output_directory = Path("/home/rt17603/work/obs/")
    
    #max_points = None # Maximum number of points per NAME file
    max_points = 50 # Maximum number of points per NAME file
    
    ## ** RERUN WITH THIS **
    regrid_method = "conservative_normed"
    
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
