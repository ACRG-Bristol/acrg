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
    
    start_date = "2019-07-26"
    end_date = "2019-07-27"
    
    area="GOSAT-SOUTHAMERICA"

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
    
    # This is broken for now - not sure why!!
    ##time_increment = None # e.g. None, "10s", "30min" (always < 1 day)
    
    #time_increment = "30min" # e.g. None, "10s", "30min" (always < 1 day)
    time_increment = "1h" # e.g. None, "10s", "30min" (always < 1 day)
    
    write_name = False
    name_directory = Path("/home/rt17603/work/TROPOMI_processed/SOUTHAMERICA/NAME/")
    
    tropomi.tropomi_process(start_date,end_date,lat_bounds,lon_bounds,coord_bin,
                            time_increment=time_increment,write_name=write_name,
                            name_directory=name_directory)
