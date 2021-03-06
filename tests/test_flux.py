#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:50:50 2019

@author: lw13938
"""
import xarray as xr
import os

from acrg.name.flux import write
from acrg.config.paths import Paths

acrg_path = Paths.acrg

def test_write(tmpdir):
    with xr.open_dataset(acrg_path / "tests/files/LPDM/emissions/EUROPE/n2o-ocean_EUROPE_2009.nc") as load:
        ds = load.load()
    del ds.attrs['author']
    del ds.attrs['date_created']
    lat = ds.lat.values
    lon = ds.lon.values
    time = ds.time.values
    flux = ds.flux.values
    species = 'n2o'
    domain = 'EUROPE'
    source= 'ocean'
    title = ds.title
    prior_info_dict = {ds.prior_file_1 : [ds.prior_file_1_version,ds.prior_file_1_raw_resolution, ds.prior_file_1_reference]}     
    write(lat, lon, time, flux, species, domain,
          source, title, prior_info_dict,
          output_directory = str(tmpdir))
    with xr.open_dataset(os.path.join(tmpdir,'EUROPE/n2o-ocean_EUROPE_2009.nc')) as load:
        testds = load.load()
    del testds.attrs['author']
    del testds.attrs['date_created']
    assert testds.equals(ds)
    assert testds.attrs == ds.attrs
    
