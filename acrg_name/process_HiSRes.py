#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 12:09:42 2019

@author: al18242

Functions used to generate and process HiSRes (High spatial resolution) footprints for NAME
"""

import acrg_name as name
import xarray as xr
import numpy as np
import numpy.ma as ma
import os
import pandas as pd
from acrg_name import emissions_helperfuncs as emfuncs
from acrg_name.name import open_ds
import xesmf
import pyproj
import acrg_grid.regrid_xesmf as regrid
import acrg_grid.areagrid as areagrid

#OSGB projection for NAEI
OS_proj = pyproj.Proj(init='epsg:7405')
#lat-lon projection
ll_proj = pyproj.Proj(init='epsg:4326')

def loadAscii(filename):
    """
    Read in ASCII GIS files for NAEI raw data
    """
    data = np.loadtxt(filename, skiprows = 6)
    data = np.flip(data, axis=0)
    header = np.loadtxt(filename,usecols=(1))[:6]
    
    data[data == float(header[5])] = np.nan
    
    xValues = np.arange(header[0]) * header[4] + header[2]
    yValues = np.arange(header[1]) * header[4] + header[3]
    
    output = xr.Dataset({})
    output["lat"] = xr.DataArray(yValues,dims="lat")
    output["lon"] = xr.DataArray(xValues,dims="lon")
    output["data"] = xr.DataArray(data,dims=["lat", "lon"])
    
    return output

def getGridCC(x, y):
    """
    Create a grid centered meshgrid for pcolormesh from 1D arrays of lon lat
    """
    dx = x[2]-x[1]
    dy = y[2]-y[1]
    x = np.append(x, x[-1] + dx)
    y = np.append(y, y[-1] + dy)
    x -= dx/2.
    y -= dy/2.
    XX, YY = np.meshgrid(x, y)
    return XX, YY

def getGridLL(x, y):
    x = np.append(x, x[-1] + x[2]-x[1])
    y = np.append(y, y[-1] + y[2]-y[1])
    XX, YY = np.meshgrid(x, y)
    return XX, YY

def process_all(domain, site, height,
                force_met_empty = False,
                satellite = False,
                upper_level = None,
                processed_folder = "Processed_Fields_files",
                processed_folder_HR = "Processed_Fields_files_HR",
                base_dir = "/data/al18242/name_out/",
                **kwargs):
    """
    Process high resolution footprints, by seperately processing the low and high resolution outputs
    These are then combined together into a single dataset ready for further use
    """
    
    #process low resolution:
    name.process.process_all(domain, site, 
                         height,
                         force_met_empty = force_met_empty,
                         satellite = satellite,
                         upper_level = upper_level,
                         processed_folder = processed_folder,
                         base_dir=base_dir, **kwargs)
    #process high resolution:
    name.process.process_all(domain, site, 
                             height,
                             force_met_empty = force_met_empty,
                             satellite = satellite,
                             upper_level = upper_level,
                             fields_folder = "Fields_files_HR",
                             processed_folder = processed_folder_HR,
                             base_dir=base_dir, **kwargs)
    
    #open processed datasets to be combined
    output_base = "{}{}_{}_{}/".format(base_dir,domain,site,height)
    filenames = os.listdir("{}{}".format(output_base, processed_folder_HR))
    for i, filename in enumerate(filenames):
        combine_date(output_base, processed_folder, processed_folder_HR, filename)
        
def process(domain, site, height, year, month,
            force_met_empty = False,
            satellite = False,
            processed_folder = "/work/al18242/fp/",#"Processed_Fields_files",
            processed_folder_HR = "/work/al18242/fp_onlyHR/",#"Processed_Fields_files_HR",
            base_dir = "/work/al18242/name_out/",
            **kwargs):
    """
    Process high resolution footprints, by seperately processing the low and high resolution outputs
    These are then combined together into a single dataset ready for further use
    """
    
    #process low resolution:
    name.process.process(domain, site, height, year, month,
                         force_met_empty = force_met_empty,
                         satellite = satellite,
                         fields_folder = "Fields_files",
                         process_dir = processed_folder,
                         base_dir=base_dir,
                         **kwargs)
    #process high resolution:
    name.process.process(domain, site, height, year, month,
                             force_met_empty = force_met_empty,
                             satellite = satellite,
                             fields_folder = "Fields_files_HR",
                             process_dir = processed_folder_HR,
                             base_dir=base_dir,
                             **kwargs)
    
    #open processed datasets to be combined
    output_base = "/work/al18242/fp_hr/EUROPE/"#"{}{}_{}_{}/".format(base_dir,domain,site,height)
    filename = "{}-{}_{}_{}{}.nc".format(site,height,domain,str(year),str(month))
    combine_date(output_base, processed_folder+domain, processed_folder_HR+domain, filename)
    
def getOverlapParameters(lat_low, lon_low, lat_high, lon_high):
    """
    calculates parameters from the low and high resolution grids for flattening and computing overlap
    """
    lats, lons = np.meshgrid(lat_low, lon_low)
    lowsize = len(lon_low) * len(lat_low)
    lons_low=lons.reshape(lowsize)
    lats_low=lats.reshape(lowsize)
    
    #process high
    lats, lons = np.meshgrid(lat_high, lon_high)
    highsize = len(lon_high) * len(lat_high)
    lons_high=lons.reshape(highsize)
    lats_high=lats.reshape(highsize)
    
    #get indicies where grids overlap
    lon_first_i = np.where(lon_low > lons_high[0])[0][0]
    lon_last_i = np.where(lon_low < lons_high[-1])[0][-1]
    
    lat_first_i = np.where(lat_low > lats_high[0])[0][0]
    lat_last_i = np.where(lat_low < lats_high[-1])[0][-1]
    
    lon_indicies = np.arange(lon_first_i, lon_last_i+1,1)
    lat_indicies = np.arange(lat_first_i, lat_last_i+1,1)
    
    indicies_to_remove_lon = np.where(np.in1d(lons_low, lon_low[lon_indicies]))[0]
    indicies_to_remove_lat = np.where(np.in1d(lats_low, lat_low[lat_indicies]))[0]
    indicies_to_remove = np.intersect1d(indicies_to_remove_lon,indicies_to_remove_lat)
    
    lons_out = np.delete(lons_low, indicies_to_remove)
    lats_out = np.delete(lats_low, indicies_to_remove)
    lons_out = np.append(lons_out,lons_high)
    lats_out = np.append(lats_out,lats_high)
    
    return (lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out)
    
def merge_resolutions(fp_low, fp_high, base, dataType = "fp"):
    if dataType == "fp":
        lowName = "fp_low"
        highName = "fp_high"
        combinedName = "fp"
    elif dataType=="flux":
        lowName = "low_res"
        highName = "high_res"
        combinedName = "flux"
    else:
        raise Exception("dateType must be flux or fp")
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(fp_low.lat.values, fp_low.lon.values, fp_high.lat_high.values, fp_high.lon_high.values)
    
    #flatten footprints ready to combine
    fp_low_flat = fp_low.values.reshape(lowsize, -1, order="F")#(["index", "time"] ,fp_low.values.reshape(lowsize, -1, order="F"))
    fp_high_flat = fp_high.values.reshape(highsize, -1, order="F")#(["index", "time"] ,fp_high.values.reshape(highsize, -1, order="F"))
    
    lons_out = np.delete(lons_low, indicies_to_remove)
    lats_out = np.delete(lats_low, indicies_to_remove)
    fp_out = np.delete(fp_low_flat, indicies_to_remove, axis=0)
    
    lons_out=np.append(lons_out,lons_high)
    lats_out=np.append(lats_out,lats_high)
    fp_out=np.append(fp_out,fp_high_flat,axis=0)
    
    #create output file from existing
    output_file = base.copy(deep=True)
    output_file[lowName] = fp_low
    output_file[highName] = (["lat_high", "lon_high", "time"], fp_high.values)
    output_file["lat_high"] = (["lat_high"], fp_high.lat_high.values, {"units":"Degrees_north"})
    output_file["lon_high"] = (["lon_high"], fp_high.lon_high.values, {"units":"Degrees_east"})
    output_file[combinedName] = (["index", "time"], fp_out)
    
    output_file["index_lons"] = (["index"], lons_out)
    output_file["index_lats"] = (["index"], lats_out)
    
    return output_file

def combine_date(output_base, processed_folder, processed_folder_HR, filename):
    """
    Combine the low and high resolution outputs into a single file
    """
    fp_low = xr.open_dataset("{}/{}".format(processed_folder, filename))
    fp_high = xr.open_dataset("{}/{}".format(processed_folder_HR, filename))
    
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(fp_low.lat.values, fp_low.lon.values, fp_high.lat.values, fp_high.lon.values)
    
    #flatten footprints ready to combine
    fp_low_flat = fp_low.copy(deep=True)
    fp_low_flat["fp"] = (["index", "time"] ,fp_low.fp.values.reshape(lowsize, -1, order="F"))
    fp_high_flat = fp_high.copy(deep=True)
    fp_high_flat["fp"] = (["index", "time"] ,fp_high.fp.values.reshape(highsize, -1, order="F"))
    
    lons_out = np.delete(lons_low, indicies_to_remove)
    lats_out = np.delete(lats_low, indicies_to_remove)
    fp_out = np.delete(fp_low_flat.fp, indicies_to_remove, axis=0)
    
    lons_out=np.append(lons_out,lons_high)
    lats_out=np.append(lats_out,lats_high)
    fp_out=np.append(fp_out,fp_high_flat.fp,axis=0)
    
    #create output file from existing
    output_file = fp_low.copy(deep=True)
    output_file["fp_low"] = fp_low.fp
    output_file["fp_high"] = (["lat_high", "lon_high", "time"], fp_high.fp.values)
    output_file["lat_high"] = (["lat_high"], fp_high.lat.values, {"units":"Degrees_north"})
    output_file["lon_high"] = (["lon_high"], fp_high.lon.values, {"units":"Degrees_east"})
    output_file["fp"] = (["index", "time"], fp_out)
    
    output_file["index_lons"] = (["index"], lons_out)
    output_file["index_lats"] = (["index"], lats_out)
    
    #save output
    out_directory = output_base#"{}Processed_Fields_files_combined".format(output_base)
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)
    print("Writing file {}/{}".format(out_directory, filename))
    output_file.to_netcdf("{}/{}".format(out_directory, filename))
    
def NAEIinEDGAR(species="ch4",
                NAEI_directory = "/data/al18242/flux/NAEI_RAW",
                sector="total",
                NAEI_year = "2017",
                EDGAR_year = "2015",
                output_directory = "/data/al18242/flux_hr",
                template_file = "/data/al18242/fp_hr/EUROPE/TMB-HR_UKV_EUROPE_201805.nc"):
    """
    Create an emissions file with NAEI values embedded within EDGAR emissions
    """
    
    naei_raw = loadAscii("{}/{}{}{}.asc".format(NAEI_directory, sector, species, NAEI_year[-2:]))
    
    if species.lower() == "ch4":
        #edgar_file = open_ds("/data/shared/Gridded_fluxes/CH4/EDGAR_v4.3.2/v432_CH4_TOTALS_nc/v432_CH4_2012.0.1x0.1.nc")
        edgar_file = open_ds("/data/shared/Gridded_fluxes/CH4/EDGAR_v5.0/yearly/v50_CH4_2015.0.1x0.1.nc")
    elif species.lower() == "co2":
        edgar_file = open_ds("/data/shared/Gridded_fluxes/CO2/EDGAR_v5.0/v50_CO2_excl_short-cycle_org_C_2015.0.1x0.1.nc")
    else:
        print("Invalid species")
        return 0;
    
    hr_template = open_ds(template_file)
    
    ################################################
    #Create the input grid
    ################################################
    
    #add 500m to make NAEI cell centered
    XX, YY = np.meshgrid(naei_raw.lon.values+500,naei_raw.lat.values+500)
    XX, YY = pyproj.transform(OS_proj,ll_proj, XX, YY)
    #create boundary mesh
    XX_b, YY_b = getGridCC(naei_raw.lon.values+500,naei_raw.lat.values+500)
    XX_b, YY_b = pyproj.transform(OS_proj,ll_proj, XX_b, YY_b)
    
    input_grid = xr.Dataset({'lon': (['x', 'y'], XX),
                             'lat': (['x', 'y'], YY),
                             'lon_b': (['x_b', 'y_b'], XX_b),
                             'lat_b': (['x_b', 'y_b'], YY_b)})
    
    ################################################
    #Create the EDGAR grid
    ################################################
    XX, YY = np.meshgrid(edgar_file.lon.values, edgar_file.lat.values)
    XX_b, YY_b = getGridCC(edgar_file.lon.values, edgar_file.lat.values)
    output_grid = xr.Dataset({'lon': (['x', 'y'], XX),
                             'lat': (['x', 'y'], YY),
                             'lon_b': (['x_b', 'y_b'], XX_b),
                             'lat_b': (['x_b', 'y_b'], YY_b)})
    
    regridder = xesmf.Regridder(input_grid, output_grid, 'conservative')
    naei_raw_regridded = regridder( np.nan_to_num(naei_raw.data,copy=True))
    
    ################################################
    #Convert units and combine grids to EDGAR
    ################################################
    combined = naei_raw_regridded * 1e-3 / (3600*24*365)
    edgar_indicies = np.where(naei_raw_regridded  == 0)
    if species.lower() == "ch4":
        combined[edgar_indicies] = edgar_file.emi_ch4.values[edgar_indicies]
        mw = 16.
    elif species.lower() == "co2":
        combined[edgar_indicies] = edgar_file.emi_co2.values[edgar_indicies]
        mw = 12.
    
    ################################################
    #Create the output grid and convert data too this grid
    ################################################
    XX, YY = np.meshgrid(hr_template.lon.values, hr_template.lat.values)
    XX_b, YY_b = getGridCC(hr_template.lon.values, hr_template.lat.values)
    output_grid2 = xr.Dataset({'lon': (['x', 'y'], XX),
                             'lat': (['x', 'y'], YY),
                             'lon_b': (['x_b', 'y_b'], XX_b),
                             'lat_b': (['x_b', 'y_b'], YY_b)})
    regridder = xesmf.Regridder(output_grid, output_grid2, 'conservative')
    naei_edgar_out = regridder( np.nan_to_num(combined,copy=True))
    naei_edgar_out = np.expand_dims(naei_edgar_out,2) / (1000 * mw) *1e6
    
    ################################################
    #Create the hi-res grid and convert data too this grid
    ################################################
    XX, YY = np.meshgrid(hr_template.lon_high.values, hr_template.lat_high.values)
    XX_b, YY_b = getGridCC(hr_template.lon_high.values, hr_template.lat_high.values)
    output_grid_hr = xr.Dataset({'lon': (['x', 'y'], XX),
                             'lat': (['x', 'y'], YY),
                             'lon_b': (['x_b', 'y_b'], XX_b),
                             'lat_b': (['x_b', 'y_b'], YY_b)})
    regridder = xesmf.Regridder(input_grid, output_grid_hr, 'conservative')
    naei_edgar_out_hr = regridder( np.nan_to_num(naei_raw.data,copy=True))
    naei_edgar_out_hr = np.expand_dims(naei_edgar_out_hr,2) * 1e-3 / (3600*24*365) / (1000 * mw) *1e6
    
    ################################################
    #1D array book-keeping
    ################################################
    
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(hr_template.lat.values, hr_template.lon.values,
                             hr_template.lat_high.values, hr_template.lon_high.values)
        
    combined = naei_edgar_out[:,:,0]
    combined = combined.reshape(lowsize, order="F")
    combined = np.delete(combined, indicies_to_remove)
    combined = np.append(combined, naei_edgar_out_hr.reshape(highsize, order="F"))
    
    data_vars = {
            "low_res":(["lat", "lon", "time"], naei_edgar_out),
            "high_res":(["lat_high", "lon_high", "time"], naei_edgar_out_hr),
            "flux":(["index", "time"], np.expand_dims(combined,axis=1)),
            "index_lat":(["index"], lats_out),
            "index_lon":(["index"], lons_out)
            }
    coords = {
            "lat":np.unique(lats_low),
            "lon":np.unique(lons_low),
            "lat_high":np.unique(lats_high),
            "lon_high":np.unique(lons_high),
            "time":[np.datetime64(NAEI_year, 'Y')]
            }
    output = xr.Dataset(data_vars=data_vars, coords=coords)
    #output.attrs["Notes"] = "2017 NAEI, with 2012 EDGAR"
    output.to_netcdf("{}/EUROPE/{}_edgar{}-naei_{}.nc".format(output_directory, species, EDGAR_year, NAEI_year))

def add_inner_grid(data, template, output_file= None):
    """
    Modify an existing footprint file to contain an inner high resolution grid. This only changes the stored grid,
    it does not increase the resolution of the underlying information. 

    Parameters
    ----------
    data : DataSet or Str
        A standard 30-day integrated footprint dataset to add an inner high resolution grid to, or the filepath to such a dataset
    template : DataSet or Str
        A footprint dataset that already contains the inner high resolution grid you want to add to data, or the filepath to such a dataset
    output_file : str
        File to save the new Dataset to. Default None results in the file not being saved

    Returns
    -------
    data_with_inner_grid : Dataset
        a copy of data with the high resolution grid added

    """
    if type(data) == str:
        data = open_ds(data)
    if type(template) == str:
        template = open_ds(template)
    
    #calculate parameters to combien grids
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(template.lat.values, template.lon.values,
                             template.lat_high.values, template.lon_high.values)
        
    #cut out inner grid to improve speeds and allow checks
    inner_grid_original = data.fp.loc[dict(lat = np.unique(lats_low[indicies_to_remove]),
                                    lon= np.unique(lons_low[indicies_to_remove]))]
    
    #footprint has implicit area information - need to scale this correctly
    area_original = areagrid(inner_grid_original.lat.values, inner_grid_original.lon.values)
    area_new = areagrid(template.lat_high.values, template.lon_high.values)
    
    #regrid and check conservation:
    inner_grid = (regrid.regrid_uniform_cc(inner_grid_original.transpose('time', 'lat', 'lon').values/area_original,
                                          inner_grid_original.lat.values, inner_grid_original.lon.values,
                                          template.lat_high.values, template.lon_high.values)*area_new).transpose([1,2,0])
    
    test_equal = np.allclose(np.sum(inner_grid_original.values, axis=(0,1)), np.sum(inner_grid, axis=(0,1)), rtol=1e-3)
    if not test_equal:
        print("WARNING: Total sensitivity has not been conserved.")

    #create the new dataset
    data_with_inner_grid = data.copy(deep=True)
    data_with_inner_grid["fp_low"] = data.fp
    
    data_with_inner_grid["fp_high"] =  xr.DataArray(inner_grid,
                                  dims=["lat_high", "lon_high", "time"],
                                  coords= {"lat_high":template.lat_high,
                                          "lon_high":template.lon_high,
                                          "time":data_with_inner_grid.time})
    
    data_with_inner_grid = merge_resolutions(data_with_inner_grid.fp_low, data_with_inner_grid.fp_high, data_with_inner_grid)
    
    if output_file is not None:
        comp = dict(zlib=True, complevel=5)
        encoding = {var: comp for var in data_with_inner_grid.data_vars}
        data_with_inner_grid.to_netcdf(output_file, encoding=encoding)
    
    return data_with_inner_grid
    
class quadTreeNode:    
    
    def __init__(self, maskedArray):
        self.maskedArray = maskedArray
        if not (type(self.maskedArray) == ma.core.MaskedArray):
            self.maskedArray = ma.array(self.maskedArray, mask = np.zeros_like(self.maskedArray, dtype=bool))
        self.XX, self.YY = np.meshgrid(np.arange(maskedArray.shape[1]), np.arange(maskedArray.shape[0]))
        self.children = []
    
    def isLeaf(self):
        if self.children:
            return False
        else:
            return True
        
    def createChildrenWithInit(self, limit, init):
        for grid in init:
            self.children.append(quadTreeNode(ma.array(self.maskedArray.data, mask = grid)))
        
        leafList = []
        self.appendLeaves(leafList)
        for leaf in leafList:
            leaf.createChildren(limit)      
        
    def createChildren(self, limit):
        value = np.sum(self.maskedArray)
        if (value < limit or
            (np.sum(~self.maskedArray.mask) <=1 )):
            return
        
        weight = ~self.maskedArray.mask
        mx = np.sum(self.XX * weight) / np.sum(weight)
        my = np.sum(self.YY * weight) / np.sum(weight)
        TL = ma.array(self.maskedArray.data, mask = self.maskedArray.mask | ~((self.XX <= mx) & (self.YY <= my)))
        TR = ma.array(self.maskedArray.data, mask = self.maskedArray.mask | ~((self.XX > mx) & (self.YY <= my)))
        BR = ma.array(self.maskedArray.data, mask = self.maskedArray.mask | ~((self.XX > mx) & (self.YY > my)))
        BL = ma.array(self.maskedArray.data, mask = self.maskedArray.mask | ~((self.XX <= mx) & (self.YY > my)))
        
        self.children.append(quadTreeNode(TL))
        self.children.append(quadTreeNode(TR))
        self.children.append(quadTreeNode(BR))
        self.children.append(quadTreeNode(BL))
        
        for child in self.children:
            child.createChildren(limit)
        
    def appendLeaves(self, leafList):
        if (self.isLeaf()):
            leafList.append(self)
        else:
            for child in self.children:
                child.appendLeaves(leafList)
        
    
def quadTreeGrid(grid, limit, init = None):
    parentNode = quadTreeNode(grid)
    if init:
        parentNode.createChildrenWithInit(limit, init)
    else:
        parentNode.createChildren(limit)
    leafList = []
    parentNode.appendLeaves(leafList)
    
    outputGrid = np.zeros_like(grid.data)
    for i, leaf in enumerate(leafList):
        outputGrid[~leaf.maskedArray.mask] = i
    
    return outputGrid
        
def makeBasisFromExisting():
    #flux = name.name.flux("EUROPE", "CH4-BTT-5", flux_directory="/data/al18242/flux_HR/")
    with xr.open_dataset("/data/al18242/flux_hr/ch4-BTT-5.nc") as ds:
        flux = ds.load()
    with xr.open_dataset("/data/shared/NAME/basis_functions/EUROPE/sub-country_mask_uk-split_EUROPE_2014.nc") as ds:
        existing = ds.load()
    basis_out = existing.copy(deep=True)
    
    #do the high res section as a simple grid for now
    start_region = np.amax(existing.basis.values) + 1
    high_res_basis = np.zeros_like(flux.high_res)
    high_res_basis[:,:,0] = start_region
    xVals = np.arange(22,55)
    yVals = np.arange(22,55)
    X, Y = np.meshgrid(yVals, xVals)
    Xc = np.floor( (X-np.amin(X))/1.0)
    Yc = np.floor( (Y-np.amin(Y))/1.0)
    rows = np.amax(Yc-np.amin(Yc))-1
    high_res_basis[X.astype(int),Y.astype(int),0] = 1+start_region + Yc + Xc * rows
    
    lowsize, highsize, lons_low, lats_low, lons_high, lats_high, indicies_to_remove, lons_out, lats_out = \
        getOverlapParameters(flux.lat_low.values, flux.lon_low.values, flux.lat_high.values, flux.lon_high.values)
        
    combined = existing.basis.values.reshape(lowsize, -1,order="F")
    combined = np.delete(combined, indicies_to_remove, axis=0)
    combined = np.append(combined, high_res_basis.reshape(highsize, -1,order="F"),axis=0)
    basis_out["basis_old"] = basis_out["basis"].copy(deep=True)
    basis_out["basis"] = (["index", "time"], combined.astype(np.int32))
    basis_out["lat_high"] = flux.lat_high
    basis_out["lon_high"] = flux.lon_high
    
    basis_out.to_netcdf("/data/al18242/basis_hr/EUROPE/sub_btt-1-london_EUROPE_2000.nc")
    
    
