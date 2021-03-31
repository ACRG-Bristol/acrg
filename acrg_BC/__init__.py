# -*- coding: utf-8 -*-
"""
Created on Thu Nov 26 17:45:56 2015

@author: chxmr

edited by vf20487
"""
from __future__ import absolute_import

import os
import pandas as pd
import numpy as np
import xarray as xr
from acrg_name.name import open_ds
from acrg_config.paths import paths
from . import climatology

data_path      = paths.data
acrg_path      = paths.acrg
cams_directory = os.path.join(data_path, 'ECMWF_CAMS', 'CAMS_inversion/')