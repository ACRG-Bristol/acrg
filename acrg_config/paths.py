# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:09:47 2019

Creates a class that stores file paths required for finding data and model output
Defaults are stored in acrg_config/templates/paths_default.yaml
To create a user-specific file, copy this to acrg_config/paths.yaml, and edit

Feel free to add other paths, as required.

@author: chxmr
"""

import yaml
from pathlib import Path

_acrg_path = Path(__file__).parents[1]
_acrg_config_path = Path(__file__).parents[0]

_user_defined_data_paths = sorted(_acrg_config_path.glob("paths.y*ml"))
if len(_user_defined_data_paths) == 0:
    _data_paths_file = _acrg_config_path / "templates/paths_default.yaml"
else:
    _data_paths_file = _user_defined_data_paths[0]

with open(_data_paths_file, 'r') as f:
    _data_paths = yaml.load(f, Loader = yaml.SafeLoader)


class paths:
    '''
    Object that contains the acrg, observation and data drive paths
    ACRG path is determined from the repo directory
    Data path is populated by acrg_config/templates/paths_default.yaml
    unless, a user-defined file is present: acrg_config/paths.yaml
    To start with, copy the paths_default.yaml to acrg_config/paths.yaml
    
    All paths are pathlib.Path objects (Python >3.4)
    
    paths.acrg: path to ACRG repo
    paths.obs: path to obs folder
    path.lpdm: path to LPDM data directory
    '''
    acrg = _acrg_path
    obs = Path(_data_paths["obs_folder"])
    lpdm = Path(_data_paths["lpdm_folder"])
