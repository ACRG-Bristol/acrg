#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import os
import argparse

import sys
sys.path.insert(0,'/user/home/bq24992/acrg/')

import acrg.satellite.tccon as tccon
import acrg.satellite.config as config

if __name__=="__main__":

    from acrg.config.paths import Paths
    acrg_path = Paths.acrg

    config_file = os.path.join(acrg_path,"acrg/satellite/templates/tccon_process.ini")
    
    parser = argparse.ArgumentParser(description='Running tccon process script')
    parser.add_argument("-c","--config",help="Configuration filename",default=config_file)
    
    args = parser.parse_args()
    config_file = args.config or args.config_file
    
    tccon_param = config.param(config_file)
    #print("Input parameters for tccon function: ",tccon_param)
    tccon_ds = tccon.tccon_process(**tccon_param)