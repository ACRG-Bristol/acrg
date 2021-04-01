#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 15:22:18 2019

@author: lw13938
"""
import pytest
import numpy as np
from acrg_name import emissions_helperfuncs as ehf
import acrg_grid
from acrg_config.paths import paths
acrg_path = paths.acrg
import os
import xarray as xr

lon =  np.arange(-5,4,0.5)
lat =  np.arange(50,55,0.5)
year = 2012

@pytest.mark.parametrize('species', [('ch4'),('n2o')])    
def test_getedgarannualtotals(species):    
    testout = ehf.getedgarannualtotals(year, lon, lat, species=species)  
    outdim = [len(lat), len(lon)]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    
def test_getothernaturalCH4():
    testout = ehf.getothernaturalCH4(lon, lat)
    outdim = [len(lat), len(lon),12]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    
def test_getsoilsinkCH4():
    testout = ehf.getsoilsinkCH4(lon, lat)
    outdim = [len(lat), len(lon),12]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)

@pytest.mark.parametrize('modeltype', [('extended'),('full')])    
def test_getBloom2017(modeltype):
    testout = ehf.getBloom2017(year, lon, lat, modeltype=modeltype)
    outdim = [len(lat), len(lon), 12]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)

@pytest.mark.skipif(os.path.exists("/dagage2/agage/metoffice/naei") == False, reason="Files on dagage2 not available")
def test_getnaeiandedgarCH4():
    testout = ehf.getnaeiandedgarCH4(lon,lat)
    outdim = [len(lat), len(lon)]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)

@pytest.mark.skipif(os.path.exists("/dagage2/agage/metoffice/naei") == False, reason="Files on dagage2 not available")
@pytest.mark.parametrize('naei_sector, species', [("roadtrans", 'ch4'),("total", 'n2o'), ("totalexcship",'ch4')])  
def test_getNAEI(naei_sector, species):
    testout = ehf.getNAEI(year, lon, lat, species, naei_sector)
    outdim = [len(lat), len(lon)]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    
def test_getedgarv5annualsectors():
    EDGARsectorlist = ["AGS","AWB","CHE","ENE","ENF","FFF","IND","IRO","MNM",
                       "PRO_COAL","PRO_GAS","PRO_OIL","PRO","RCO","REF_TRF","SWD_INC",
                       "SWD_LDF","TNR_Aviation_CDS","TNR_Aviation_CRS",
                       "TNR_Aviation_LTO","TNR_Other","TNR_Ship","TRO","WWT"]   
    testout = ehf.getedgarv5annualsectors(year, lon, lat, EDGARsectorlist, species='CH4')
    outdim = [len(lat), len(lon)]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    
def test_getedgarv432annualsectors():
    secdict = {'powerindustry' : '1A1a', 
               'oilrefineriesandtransformationindustry' : '1A1b_1A1c_1A5b1_1B1b_1B2a5_1B2a6_1B2b5_2C1b',
               'combustionformanufacturing' : '1A2',
               'aviationclimbinganddescent' : '1A3a_CDS',
               'aviationcruise' : '1A3a_CRS',
               'aviationlandingandtakeoff' : '1A3a_LTO',
               'aviationsupersonic' : '1A3a_SPS',
               'roadtransport' : '1A3b',
               'railwayspipelinesandoffroadtransport' : '1A3c_1A3e',
               'shipping' : '1A3d_1C2',
               'energyforbuildings' : '1A4',
               'fuelexploitation' : '1B1a_1B2a1_1B2a2_1B2a3_1B2a4_1B2c',
               'nonmetallicmineralsproduction' : '2A',
               'chemicalprocesses': '2B',
               'ironandsteelproduction' : '2C1a_2C1c_2C1d_2C1e_2C1f_2C2',
               'nonferrousmetalsproduction' : '2C3_2C4_2C5',
               'nonenergyuseoffuels' : '2G',
               'solventsandproductsuse' :  '3',
               'entericfermentation' : '4A',
               'manuremanagement' : '4B',
               'agriculturalsoils' : '4C_4D',
               'indirectN2Oemissionsfromagriculture' : '4D3',
               'agriculturalwasteburning' : '4F',
               'solidwastelandfills' : '6A_6D',
               'wastewaterhandling' : '6B',
               'Solid waste incineration' : '6C',
               'fossilfuelfires' : '7A',
               'indirectemissionsfromNOxandNH3' : '7B_7C'           
    }
    sectors= secdict.keys()    
    testout = ehf.getedgarv432annualsectors(year, lon, lat, sectors, species='CH4')
    outdim = [len(lat), len(lon)]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    
def test_getedgarmonthlysectors():
    secdict = {'powerindustry' : '1A1a', 
               'oilrefineriesandtransformationindustry' : '1A1b_1A1c_1A5b1_1B1b_1B2a5_1B2a6_1B2b5_2C1b',
               'combustionformanufacturing' : '1A2',
               'aviationclimbinganddescent' : '1A3a_CDS',
               'aviationcruise' : '1A3a_CRS',
               'aviationlandingandtakeoff' : '1A3a_LTO',
               'aviationsupersonic' : '1A3a_SPS',
               'roadtransport' : '1A3b',
               'railwayspipelinesandoffroadtransport' : '1A3c_1A3e',
               'shipping' : '1A3d_1C2',
               'energyforbuildings' : '1A4',
               'fuelexploitation' : '1B1a_1B2a1_1B2a2_1B2a3_1B2a4_1B2c',
               'nonmetallicmineralsproduction' : '2A',
               'chemicalprocesses': '2B',
               'ironandsteelproduction' : '2C1a_2C1c_2C1d_2C1e_2C1f_2C2',
               'nonferrousmetalsproduction' : '2C3_2C4_2C5',
               'nonenergyuseoffuels' : '2G',
               'solventsandproductsuse' :  '3',
               'entericfermentation' : '4A',
               'manuremanagement' : '4B',
               'agriculturalsoils' : '4C_4D',
               'indirectN2Oemissionsfromagriculture' : '4D3',
               'agriculturalwasteburning' : '4F',
               'solidwastelandfills' : '6A_6D',
               'wastewaterhandling' : '6B',
               'Solid waste incineration' : '6C',
               'fossilfuelfires' : '7A',
               'indirectemissionsfromNOxandNH3' : '7B_7C'           
    }
    sectors= secdict.keys() 
    months = [5,6,7]
    testout = ehf.getedgarmonthlysectors(lon, lat, sectors, months=months,
                           species='CH4')
    outdim = [len(lat), len(lon), len(months)]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    
@pytest.mark.parametrize('scarpelli_sector', [('coal'), ('oil'), ('gas'), ('all')])  
def test_getScarpelliFossilFuelsCH4(scarpelli_sector):
    testout = ehf.getScarpelliFossilFuelsCH4(lon, lat, scarpelli_sector=scarpelli_sector)
    outdim = [len(lat), len(lon)]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)

@pytest.mark.parametrize('scale_wetlands, total_w_emission', [(True, 185e12), (True, 100e12), (False, 185e12)]) 
def test_getJULESwetlands(scale_wetlands, total_w_emission):
    testout = ehf.getJULESwetlands(year,lon,lat,scale_wetlands=scale_wetlands,
                     total_w_emission=total_w_emission)
    outdim = [len(lat), len(lon), 12]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)

@pytest.mark.parametrize('timeframe, nt', [('daily', 365), ('monthly', 12)])     
def test_getbloomwetlandsCH4(timeframe, nt):
    testout = ehf.getbloomwetlandsCH4(year, lon, lat, timeframe=timeframe)
    outdim = [len(lat), len(lon), nt]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    

@pytest.mark.parametrize('timeframe, species', [('monthly', 'ch4'), ('daily', 'ch4'), ('3hourly', 'ch4'), ('monthly', 'co2')])  
def test_getGFED(monkeypatch, timeframe, species):
    if timeframe == 'monthly':
        nt = 2
    elif timeframe == 'daily':
        nt = 62
    elif timeframe == '3hourly':
        nt = 496
    def mockreturn(narr):
        return np.random.rand(len(lat), len(lon), nt)
    monkeypatch.setattr(acrg_grid.regrid, 'regrid2d', mockreturn)
    testout = ehf.getGFED(year, lon, lat, timeframe=timeframe, months = [7,8], species=species)
    outdim = [len(lat), len(lon), nt]
    assert np.isfinite(testout).all() 
    assert np.array_equal(testout.shape, outdim)
    
def test_get_naei_public():
    filename = os.path.join(acrg_path, "tests/files/gridded_fluxes/raw_naei_mock.asc")
    lat_out = np.arange(48.0, 58.1, 2.0)
    lon_out = np.arange(-9.0,3.1, 2.0)
    output = ehf.get_naei_public(filename, lon_out, lat_out)
    output_sum = np.sum(output["data"] * acrg_grid.areagrid(lat_out, lon_out))
    
    raw = ehf.loadAscii(os.path.join(acrg_path, "tests/files/gridded_fluxes/raw_naei_mock.asc"))
    raw_sum = np.sum(raw["data"] * raw.attrs["cellsize"] * raw.attrs["cellsize"])
    
    assert np.abs( (output_sum-raw_sum)/raw_sum ) < 0.01
    
    
def test_processUKGHG():
    ukghg_dir = os.path.join(acrg_path,"tests/files/gridded_fluxes/")
    lat_out = np.arange(50.,65.,1.)
    lon_out = np.arange(-10.,5.,1.0)

    out = ehf.processUKGHG(ukghg_dir,"ch4","2015","total",lon_out,lat_out)
    out_total = np.nansum(out["flux"].values * acrg_grid.areagrid(out.lat.values,
                                                                  out.lon.values))
    
    with xr.open_dataset(os.path.join(ukghg_dir,"uk_flux_total_ch4_LonLat_0.01km_2015.nc")) as raw_f:
        raw_total = np.nansum(raw_f["ch4_flux"].values[0,:,:] * acrg_grid.areagrid(raw_f.latitude.values,
                                                                         raw_f.longitude.values))
    
    assert np.abs((out_total - raw_total)/raw_total) < 0.01
    
def test_get_US_EPA():

    out = ehf.get_US_EPA(["6A_Landfills_Municipal"])
    out_total = np.sum(out["flux"].values[:,:,0] * (acrg_grid.areagrid(out.lat.values,
                                                                       out.lon.values)))
                                                  
    with xr.open_dataset(os.path.join(acrg_path,"tests/files/gridded_fluxes/GEPA_Annual.nc")) as raw_f:
        raw_flux = raw_f["emissions_6A_Landfills_Municipal"].values / 6.02214076e23 * 1e4
        raw_total = np.sum(raw_flux * acrg_grid.areagrid(raw_f.lat.values,raw_f.lon.values))
                          
    assert np.abs((out_total - raw_total)/raw_total) < 0.01