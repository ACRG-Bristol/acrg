# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 10:45:51 2014

"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import datetime as dt
import os
import glob
from matplotlib import ticker
import pandas as pd
import bisect
import subprocess
from progressbar import ProgressBar
import json
import acrg_agage as agage
#import acrg_regrid as regrid
import acrg_convert as unit_convert
from acrg_grid import areagrid
import xray
from os.path import split, realpath, join
from acrg_time import convert
import calendar
import pickle

acrg_path = os.getenv("ACRG_PATH")
data_path = os.getenv("DATA_PATH")

if acrg_path is None:
    acrg_path = os.getenv("HOME")
    print("Default ACRG directory is assumed to be home directory. Set path in .bashrc as \
            export ACRG_PATH=/path/to/acrg/repository/ and restart python terminal")
if data_path is None:
    data_path = "/data/shared/"
    print("Default Data directory is assumed to be /data/shared/. Set path in .bashrc as \
            export DATA_PATH=/path/to/data/directory/ and restart python terminal")

fp_directory = join(data_path, 'NAME/fp/')
flux_directory = join(data_path, 'NAME/emissions/')
basis_directory = join(data_path, 'NAME/basis_functions/')
bc_directory = join(data_path, 'NAME/bc/')
bc_basis_directory = join(data_path,'NAME/bc_basis_functions/')
fp_HiTRes_directory = join(data_path,'NAME/fp_high_time_res/')

# Get acrg_site_info file
with open(join(acrg_path, "acrg_site_info.json")) as f:
    site_info=json.load(f)

def filenames(site, domain, start, end, height = None, flux=None, basis=None, HiTRes = False):
    """
    Output a list of available footprint file names,
    for given site, domain, directory and date range.

    Doesn't work for flux or basis functions yet.
    """
    if flux is None and basis is None and HiTRes is False:
        baseDirectory = fp_directory
    else:
        if flux is not None:
            baseDirectory = flux_directory
        if basis is not None:
            baseDirectory = basis_directory
        if HiTRes is not False:
            baseDirectory = fp_HiTRes_directory
    
    # Get height
    #Get site info for heights
    if height is None:
        if not site in site_info.keys():
            print("Site code not found in arcg_site_info.json to get height information. " + \
                  "Check that site code is as intended. "+ \
                  "If so, either add new site to file or input height manually.")
            return None
        height = site_info[site]["height_name"][0]
    
    # Generate time series
    months = pd.DatetimeIndex(start = start, end = end, freq = "M").to_pydatetime()
    yearmonth = [str(d.year) + str(d.month).zfill(2) for d in months]

    files = []
    for ym in yearmonth:
        f=glob.glob(baseDirectory + \
            domain + "/" + \
            site + "*" + "-" + height + "*" + domain + "*" + ym + "*.nc")
#        f=glob.glob(baseDirectory + \
#            domain + "/" + \
#            site + "*" + height + "*" + domain + "*" + ym + "*.nc")
        if len(f) > 0:
            files += f

    files.sort()

    if len(files) == 0:
        print("Can't find file: " + baseDirectory + \
            domain + "/" + \
            site + "*" + height + "*" + domain + "*" + "*.nc")
    return files

def read_netcdfs(files, dim = "time", transform_func=None):
    '''
    Use xray to open sequential netCDF files. 
    Makes sure that file is closed after open_dataset call.
    '''
    
    def process_one_path(path):
        with xray.open_dataset(path) as ds:
            ds.load()
        return ds
    
    datasets = [process_one_path(p) for p in sorted(files)]
    combined = xray.concat(datasets, dim)
    return combined   


def footprints(sitecode_or_filename, start = "2010-01-01", end = "2016-01-01",
        domain="EUROPE", height = None, species = None, emissions_name = None, HiTRes = False):
    """
    Load a NAME footprint netCDF files into an xray dataset.
    Either specify:
    
    a) A file name:

        fp = footprints(filename)
    
    b) A site code, domain, and date range:
    
        fp = footprints("MHD", 
                    start = "2014-01-01", end = "2014-01-01",
                    domain = "EUROPE")
    
    If the HEIGHT keyword is not specified, the default height from the
    acrg_site_info.json file is assumed.
    
    If the SPECIES or EMISSIONS_NAME keyword is given, fluxes will also
    be loaded and merged into the dataset.
    
    EMISSIONS_NAME allows emissions files such as co2nee_EUROPE_2012.nc
    to be read in. In this case EMISSIONS_NAME would be 'co2nee'
    """
    
    #Chose whether we've input a site code or a file name
    #If it's a three-letter site code, assume it's been processed
    # into an annual footprint file in (mol/mol) / (mol/m2/s)
    # using acrg_name_process
    if '.nc' in sitecode_or_filename:
        if not '/' in sitecode_or_filename:
            files = [os.path.join(fp_directory, sitecode_or_filename)]
        else:
            files=[sitecode_or_filename]
    else:
        site=sitecode_or_filename[:]
        files = filenames(site, domain, start, end, height = height)

    if len(files) == 0:
        print("Can't find files, " + sitecode_or_filename)
        return None

    else:
        fp=read_netcdfs(files)
        
        # If a species is specified, also get flux and vmr at domain edges
        if emissions_name is not None:
            flux_ds = flux(domain, emissions_name)
            if flux_ds is not None:
                fp = combine_datasets(fp, flux_ds)
        elif species is not None:
            flux_ds = flux(domain, species)
            if flux_ds is not None:
                fp = combine_datasets(fp, flux_ds)
        
        if species is not None:
            bc_ds = boundary_conditions(domain, species)
            if bc_ds is not None:
                fp = combine_datasets(fp, bc_ds)
        
        if HiTRes == True:
            HiTRes_files = filenames(site, domain, start, end, height = height, HiTRes=True)
            HiTRes_ds = read_netcdfs(HiTRes_files)
            fp = combine_datasets(fp, HiTRes_ds)

        return fp


def flux(domain, species):
    """
    Read in a flux dataset.
    To be consistent with the footprints, fluxes should be in
    mol/m2/s. 
    Note that at present ALL flux data is read in per species per domain
    this may get slow for very large flux datasets, and we may want to subset.    
    """

    files = sorted(glob.glob(flux_directory + domain + "/" + 
                   species.lower() + "_" + "*.nc"))
    if len(files) == 0:
        print("Can't find flux: " + domain + " " + species)
        return None
    
    flux_ds = read_netcdfs(files)
    
    # Check that time coordinate is present
    if not "time" in flux_ds.coords.keys():
        print("ERROR: No 'time' coordinate " + \
              "in flux dataset for " + domain + ", " + species)
        return None

    # Check for level coordinate. If one level, assume surface and drop
    if "lev" in flux_ds.coords.keys():
        print("WARNING: Can't support multi-level fluxes. Trying to remove 'lev' coordinate " + \
              "from flux dataset for " + domain + ", " + species)
        if len(flux_ds.lev) > 1:
            print("ERROR: More than one flux level")
        else:
            return flux_ds.drop("lev")

    return flux_ds


def boundary_conditions(domain, species):
    """
    Read in the files with the global model vmrs at the domain edges to give
    the boundary conditions.
    """
    
    files = sorted(glob.glob(bc_directory + domain + "/" + 
                   species.lower() + "_" + "*.nc"))
    if len(files) == 0:
        print("Can't find boundary condition files: " + domain + " " + species)
        return None

    bc_ds = read_netcdfs(files)

    return bc_ds


def basis(domain, basis_case = 'voronoi'):
    """
    Read in a basis function file.
    """
    
    files = sorted(glob.glob(basis_directory + domain + "/" +
                    basis_case + "*.nc"))
    if len(files) == 0:
        print("Can't find basis functions: " + domain + " " + basis_case)
        return None

    basis_ds = read_netcdfs(files)

    return basis_ds


def basis_boundary_conditions(domain, basis_case = 'NESW'):
    
    files = sorted(glob.glob(bc_basis_directory + domain + "/" +
                    basis_case + "_*.nc"))
    if len(files) == 0:
        print("Can't find boundary condition basis functions: " + domain + " " + basis_case)
        return None
    
    basis_ds = read_netcdfs(files)

    return basis_ds


def combine_datasets(dsa, dsb, method = "nearest", tolerance = None):
    """
    Merge two datasets. Assumes that you want to 
    re-index to the FIRST dataset.
    
    Example:
    
    ds = combine_datasets(dsa, dsb)

    ds will have the index of dsa    
    """
    # merge the two datasets within a tolerance and remove times that are NaN (i.e. when FPs don't exist)
    ds_temp = dsa.merge(dsb.reindex_like(dsa, method, tolerance = tolerance))
    if 'fp' in ds_temp.keys():
        flag = np.where(np.isfinite(ds_temp.fp.mean(dim=["lat","lon"]).values))
        ds_temp = ds_temp[dict(time = flag[0])]
    return ds_temp

def timeseries(ds):
    """
    Compute flux * footprint time series.
    All that is required is that you input an xray
    dataset with both the flux and footprint fields present    
    
    Example:
    
        ts = timeseries(dataset)
    
    There are almost certainly much more efficient ways of doing this.
    """

    if "flux" in ds.keys():
        return (ds.fp*ds.flux).sum(["lat", "lon"])
    else:
        print("Can't calculate time series " + \
              "no fluxes. Check flux file.")
        return None

def timeseries_HiTRes(fp_HiTRes_ds, domain, HiTRes_flux_name, Resid_flux_name,
                      output_TS = True, output_fpXflux = True):
    flux_HiTRes = flux(domain, HiTRes_flux_name)
    flux_resid = flux(domain, Resid_flux_name)
    
    fp_HiTRes = fp_HiTRes_ds.fp_HiTRes.to_dataset()
    fpXflux = np.zeros((len(fp_HiTRes.lat), len(fp_HiTRes.lon), len(fp_HiTRes.time)))
    
    for ti, time in enumerate(fp_HiTRes.time):
        fp = fp_HiTRes.sel(time=time).fp_HiTRes.to_dataset()
        time_back = [fp.time.values - np.timedelta64(i,'h') for i in fp.H_back.values]
        fp = fp.update({'H_back':time_back})
        fp = fp.drop('time')
        fp = fp.rename({'H_back':'time'})
        #To make  Hour Back' time go forward
        fp= fp.update({'fp_HiTRes' : fp.fp_HiTRes[:,:,::-1], 'time' : fp.time[::-1]})
        em = flux_HiTRes.reindex_like(fp, method='ffill')
        #Use end of hours back as closest point for finding the emissions file
        emend = flux_resid.sel(time = fp.time[0], method = 'nearest')
        em.flux[:,:,0] = emend.flux
        fpXflux[:,:,ti] = (fp.fp_HiTRes*em.flux).sum(["time"])
        
    timeseries= np.sum(fpXflux, axis = (0,1))
    
    if output_fpXflux == True and output_TS ==True:
        return timeseries, fpXflux
    
    elif output_fpXflux == False and output_TS ==True:
        return timeseries
        
    elif output_fpXflux == True and output_TS ==False:
        return fpXflux       

def timeseries_boundary_conditions(ds):
    """
    Compute particle location * global model edges time series.
    All that is required is that you input an xray
    dataset with both the particle locations and vmr at domain edge fields present    
    """ 

    return (ds.particle_locations_n*ds.vmr_n).sum(["height", "lon"]) + \
           (ds.particle_locations_e*ds.vmr_e).sum(["height", "lat"]) + \
           (ds.particle_locations_s*ds.vmr_s).sum(["height", "lon"]) + \
           (ds.particle_locations_w*ds.vmr_w).sum(["height", "lat"])

    
def footprints_data_merge(data, domain = "EUROPE", species = "CH4", load_flux = True,
                          calc_timeseries = True, calc_bc = True, HiTRes = False,
                          average = None, site_modifier = {}, height = None,
                          emissions_name = None, 
                          perturbed=False, fp_dir_pert=None, pert_year=None, pert_month=None):
    """
    Output a dictionary of xray footprint datasets, that correspond to a given
    dictionary of Pandas dataframes, containing mole fraction time series.
    
    Example:
    
    Input dictionary contains time series at Mace Head and Tacolneston:
        
        data = {"MHD": MHD_dataframe, "TAC": TAC_dataframe}

    The dataset must be labeled with "time" index, "mf" and "dmf" columns.
    To combine this with corresponding NAME footprints:
            if not site in site_info.keys():
        dataset = footprints_data_merge(data)
        
    An optional site modifier dictionary is used that maps the site name in the
    obs file to the site name in the footprint file, if they are different. This
    is useful for example if the same site FPs are run with a different met and 
    there named slightly differently from the obs file.
        
    Output dataset will contain a dictionary of merged data and footprints:
        
        dataset = {"MHD": MHD_xray_dataset, "TAC": TAC_xray_dataset}
    """

    sites = [key for key in data.keys() if key[0] != '.']
    attributes = [key for key in data.keys() if key[0] == '.']
    
    if average is not None:
        if type(average) is not list:
            average = [average]
        if len(average) != len(sites):
            print("WARNING: average list must be the same length as " + \
                  "number of sites. Ignoring. Output dataset will not be resampled.")
            average = [None for i in sites]
    else:
        average = [None for i in sites]

    # If not given, check if species is defined in data dictionary:
    if species is None:
        if ".species" in data.keys():
            species = data[".species"]
        else:
            print "Species is not specified and can't be found in data dictionary."

    # Output array
    fp_and_data = {}
    
    for si, site in enumerate(sites):

        # Dataframe for this site            
        site_df = data[site] 
            
        # Get time range
        df_start = min(site_df.index).to_pydatetime()
        start = dt.datetime(df_start.year, df_start.month, 1, 0, 0)
        
        df_end = max(site_df.index).to_pydatetime()
        month_days = calendar.monthrange(df_end.year, df_end.month)[1]
        end = dt.datetime(df_end.year, df_end.month, 1, 0, 0) + \
                dt.timedelta(days = month_days)
      
        # Convert to dataset
        site_ds = xray.Dataset.from_dataframe(site_df)
        
        if site in site_modifier.keys():
            site_modifier_fp = site_modifier[site]
        else:    
            site_modifier_fp = site
         
        if height is not None:
            
            if type(height) is not dict:
                print("Height input needs to be a dictionary with sitename:height")
                return None
                
            height_site = height[site] 
        else:
            height_site = height
        
        # Get footprints
        if perturbed:
            fp_dir_pert2=fp_dir_pert[site]
            site_modifier_fp = fp_dir_pert2 + str(site) + '-' + str(height_site) + 'magl_EUROPE_' + str(pert_year) + str(pert_month) + '.nc'
            
            site_fp = footprints(site_modifier_fp, start = start, end = end,
                             domain = domain,
                             species = [species if load_flux == True or \
                                         calc_timeseries == True or \
                                         calc_bc == True \
                                         else None][0], \
                             height = height_site,
                             emissions_name = [emissions_name if load_flux == True or \
                                         calc_timeseries == True \
                                         else None][0],
                             HiTRes = HiTRes)

        else:
            site_fp = footprints(site_modifier_fp, start = start, end = end,
                             domain = domain,
                             species = [species if load_flux == True or \
                                         calc_timeseries == True or \
                                         calc_bc == True \
                                         else None][0], \
                             height = height_site,
                             emissions_name = [emissions_name if load_flux == True or \
                                         calc_timeseries == True \
                                         else None][0],
                             HiTRes = HiTRes)
                        
        if site_fp is not None:
                        
        # If satellite data, check that the max_level in the obs and the max level in the processed FPs are the same
        # Set tolerance tin time to merge footprints and data            
            if "GOSAT" in site.upper():
                ml_obs = site_df.max_level
                ml_fp = site_fp.max_level
                tolerance = 60e9 # footprints must match data with this tolerance in [ns]
                if ml_obs != ml_fp:
                    print "ERROR: MAX LEVEL OF SAT OBS DOES NOT EQUAL MAX LEVEL IN FP"
                    print "max_level_fp =",ml_fp
                    print "max_level_obs =",ml_obs
                    return None
            elif "GAUGE-FERRY" in site.upper():
                tolerance = '5min'
            else:
                tolerance = None
                
            site_ds = combine_datasets(site_ds, site_fp,
                                       method = "nearest",
                                       tolerance = tolerance)

            
            # If units are specified, multiply by scaling factor
            if ".units" in attributes:
                site_ds.update({'fp' : (site_ds.fp.dims, site_ds.fp / data[".units"])})
                if calc_bc:
                    for key in site_ds.keys():
                        if "vmr" in key:
                            site_ds.update({key :
                                            (site_ds[key].dims, site_ds[key] / \
                                            data[".units"])})

                if HiTRes:
                    site_ds.update({'fp_HiTRes' : (site_ds.fp_HiTRes.dims, site_ds.fp_HiTRes / data[".units"])})

               
            # Calculate model time series, if required
            if calc_timeseries:
                site_ds["mf_mod"] = timeseries(site_ds)
            
            # Calculate boundary conditions, if required         
            if calc_bc:
                site_ds["bc"] = timeseries_boundary_conditions(site_ds)  
            
            # Resample, if required
            if average[si] is not None:
                site_ds = site_ds.resample(average[si], dim = "time")
            
            fp_and_data[site] = site_ds
        
    for a in attributes:
        fp_and_data[a] = data[a]

    return fp_and_data


def fp_sensitivity(fp_and_data, domain = 'EUROPE', basis_case = 'voronoi',
                   HiTRes_flux_name = None, Resid_flux_name=None):
    """
    Adds a sensitivity matrix, H, to each site xray dataframe in fp_and_data.
    
    Basis function data in an array: lat, lon, no. regions. In each 'region'
    element of array there is a lt lon grid with 1 in region and 0 outside region.
    
    Region numbering must start from 1
    """    
    
    sites = [key for key in fp_and_data.keys() if key[0] != '.']
    attributes = [key for key in fp_and_data.keys() if key[0] == '.']
    
    basis_func = basis(domain = domain, basis_case = basis_case)
    
    for site in sites:

#        site_bf = combine_datasets(fp_and_data[site]["fp", "flux", "mf_mod"],
#                                   basis_func)

        if 'fp_HiTRes' in fp_and_data[site].keys():
            site_bf_temp = xray.Dataset({"fp":fp_and_data[site]["fp"],
                                         "fp_HiTRes":fp_and_data[site]["fp_HiTRes"],
                                         "flux":fp_and_data[site]["flux"]})
        else:
            site_bf_temp = xray.Dataset({"fp":fp_and_data[site]["fp"],
                                         "flux":fp_and_data[site]["flux"]})
        
        site_bf = combine_datasets(site_bf_temp,basis_func)

        basis_scale = xray.Dataset({'basis_scale': (['lat','lon','time'],
                                                    np.zeros(np.shape(site_bf.basis)))},
                                   coords = site_bf.coords)
        site_bf = site_bf.merge(basis_scale)

        H = np.zeros((int(np.max(site_bf.basis)),len(site_bf.time)))
        
        if 'fp_HiTRes' in site_bf.keys():
            H_all_arr=timeseries_HiTRes(site_bf, domain, HiTRes_flux_name, Resid_flux_name, output_TS = False, output_fpXflux = True)
            H_all = xray.DataArray(H_all_arr, coords=[site_bf.lat, site_bf.lon, site_bf.time], dims = ['lat','lon','time'])
        else:
            H_all=site_bf.fp*site_bf.flux 
            H_all_arr = H_all.values
            
        H_all_v=np.zeros((len(site_bf.lat)*len(site_bf.lon),len(site_bf.time)))
        for ti in range(len(site_bf.time)):
            H_all_v[:,ti]=np.ravel(H_all[:,:,ti])
            
        base_temp = np.ravel(site_bf.basis.values[:,:,0])
        for i in range(int(np.max(site_bf.basis))):
            wh_ri = np.where(base_temp == i+1)
             
            for ti in range(len(site_bf.time)):
                 H[i,ti]=np.sum(H_all_v[wh_ri,ti])        
        
#        for i in range(int(np.max(site_bf.basis))):
#            site_bf.basis_scale.values = np.zeros(np.shape(site_bf.basis_scale))
#            site_bf.basis_scale.values[np.where(site_bf.basis == i+1)] = 1
#            fpalign, fluxalign, scalealign = xray.align(site_bf.fp,
#                                                        site_bf.flux,
#                                                        site_bf.basis_scale)
#            H[i,:] = (fpalign*fluxalign*scalealign).sum(["lat", "lon"])
        
        sensitivity = xray.Dataset({'H': (['region','time'], H)},
                                    coords = {'region' : range(1,np.max(site_bf.basis)+1),
                                              'time' : (fp_and_data[site].coords['time'])})

        fp_and_data[site] = fp_and_data[site].merge(sensitivity)
        
        if any([word in basis_case for word in ['transd','test', 'alcompare', 'sense', 'mcf']]):
            sub_fp_temp = site_bf.fp.sel(lon=slice(min(site_bf.sub_lon),max(site_bf.sub_lon)), 
                                    lat=slice(min(site_bf.sub_lat),max(site_bf.sub_lat))) 
            sub_fp = xray.Dataset({'sub_fp': (['sub_lat','sub_lon','time'], sub_fp_temp)},
                               coords = {'sub_lat': (site_bf.coords['sub_lat']),
                                         'sub_lon': (site_bf.coords['sub_lon']),
                                'time' : (fp_and_data[site].coords['time'])})
                                
            sub_H_temp = H_all.sel(lon=slice(min(site_bf.sub_lon),max(site_bf.sub_lon)), 
                                    lat=slice(min(site_bf.sub_lat),max(site_bf.sub_lat)))                             
            sub_H = xray.Dataset({'sub_H': (['sub_lat','sub_lon','time'], sub_H_temp)},
                               coords = {'sub_lat': (site_bf.coords['sub_lat']),
                                         'sub_lon': (site_bf.coords['sub_lon']),
                                'time' : (fp_and_data[site].coords['time'])})
                                
            fp_and_data[site] = fp_and_data[site].merge(sub_fp)
            fp_and_data[site] = fp_and_data[site].merge(sub_H)
            
    return fp_and_data


def bc_sensitivity(fp_and_data, domain = 'EUROPE', basis_case = 'NESW'):
    
    sites = [key for key in fp_and_data.keys() if key[0] != '.']
#    attributes = [key for key in fp_and_data.keys() if key[0] == '.']
    basis_func = basis_boundary_conditions(domain = domain,
                                           basis_case = basis_case)
# sort basis_func into time order    
    ind = basis_func.time.argsort()                                        
    timenew = basis_func.time[ind]
    basis_func = basis_func.reindex({"time":timenew})
    
    for site in sites:

        # stitch together the particle locations, vmrs at domain edges and
        #boundary condition basis functions
#        DS = combine_datasets(fp_and_data[site]["particle_locations_n",
#                                                "particle_locations_e",
#                                                "particle_locations_s",
#                                                "particle_locations_w",
#                                                "vmr_n",
#                                                "vmr_e",
#                                                "vmr_s",
#                                                "vmr_w",
#                                                "bc"],
#                                                basis_func)
        DS_temp = xray.Dataset({"particle_locations_n":fp_and_data[site]["particle_locations_n"],
                                "particle_locations_e":fp_and_data[site]["particle_locations_e"],
                                "particle_locations_s":fp_and_data[site]["particle_locations_s"],
                                "particle_locations_w":fp_and_data[site]["particle_locations_w"],
                                "vmr_n":fp_and_data[site]["vmr_n"],
                                "vmr_e":fp_and_data[site]["vmr_e"],
                                "vmr_s":fp_and_data[site]["vmr_s"],
                                "vmr_w":fp_and_data[site]["vmr_w"],
                                "bc":fp_and_data[site]["bc"]})
                        
        DS = combine_datasets(DS_temp, basis_func)                                    

        part_loc = np.hstack([DS.particle_locations_n,
                                DS.particle_locations_e,
                                DS.particle_locations_s,
                                DS.particle_locations_w])
        # Added by ML but not yet uploaded to repository:
        for ii in range(len(DS.time)):
            DS.vmr_n[:,:,ii]=DS.vmr_n[:,:,0].values
            DS.vmr_e[:,:,ii]=DS.vmr_e[:,:,0].values
            DS.vmr_s[:,:,ii]=DS.vmr_s[:,:,0].values
            DS.vmr_w[:,:,ii]=DS.vmr_w[:,:,0].values
            #########################################
        vmr_ed = np.hstack([DS.vmr_n,
                           DS.vmr_e,
                           DS.vmr_s,
                           DS.vmr_w])
        
        bf = np.hstack([DS.bc_basis_n,
                        DS.bc_basis_e,
                        DS.bc_basis_s,
                        DS.bc_basis_w])
        
        H_bc = np.zeros((len(DS.coords['region']),len(DS.bc)))
        
        for i in range(len(DS.coords['region'])):
            reg = bf[:,:,i,:]
            H_bc[i,:] = np.sum((part_loc*vmr_ed*reg), axis=(0,1))
        
        sensitivity = xray.Dataset({'H_bc': (['region_bc','time'], H_bc)},
                                    coords = {'region_bc': (DS.coords['region'].values),
                                              'time' : (DS.coords['time'])})

        fp_and_data[site] = fp_and_data[site].merge(sensitivity)
    
    return fp_and_data


def merge_sensitivity(fp_data_H,
                      out_filename = None,
                      remove_nan = True):
    """
    Outputs y, y_site, y_time in a single array for all sites
    (as opposed to a dictionary) and H and H_bc if present in dataset.
    Assumes my default that NaN values in y will be removed.
    """

    y = []
    y_error = []
    y_site = []
    y_time = []
    H = []
    H_bc = []
    
    sites = [key for key in fp_data_H.keys() if key[0] != '.']
    
    for si, site in enumerate(sites):

        y.append(fp_data_H[site].mf.values)
        
        # Approximate y_error
        if "vmf" in fp_data_H[site].keys():
            y_error.append(fp_data_H[site].vmf.values)
        elif "dmf" in fp_data_H[site].keys():
            y_error.append(fp_data_H[site].dmf.values)
        
        y_site.append([site for i in range(len(fp_data_H[site].coords['time']))])

        y_time.append(fp_data_H[site].coords['time'].values)

        if 'H' in fp_data_H[site].data_vars:
            # Make sure H matrices are aligned in the correct dimensions
            if fp_data_H[site].H.dims[0] == "time":
                H.append(fp_data_H[site].H.values)
            else:
                H.append(fp_data_H[site].H.values.T)
        if 'H_bc' in fp_data_H[site].data_vars:
            if fp_data_H[site].H_bc.dims[0] == "time":
                H_bc.append(fp_data_H[site].H_bc.values)
            else:
                H_bc.append(fp_data_H[site].H_bc.values.T)

    y = np.hstack(y)
    y_error = np.hstack(y_error)
    y_site = np.hstack(y_site)
    y_time = np.hstack(y_time)
    
    if remove_nan:
        why = np.isfinite(y)
        y = y[why]
        y_error = y_error[why]
        y_site = y_site[why]
        y_time = y_time[why]
    
    if len(H) > 0:
        H = np.vstack(H)
        if remove_nan:
            H = H[why, :]
            
    if len(H_bc) > 0:
        H_bc = np.vstack(H_bc)
        if remove_nan:
            H_bc = H_bc[why, :]
    
    out_variables = y, y_error, y_site, y_time
    
    if len(H_bc) > 0:
        out_variables += (H_bc,)
    else:
        out_variables += (None,)
        
    if len(H) > 0:
        out_variables += (H,)
    else:
        out_variables += (None,)

    # Save or return y, y_error, y_site, y_time, H_bc, H
    if out_filename is None:
        return out_variables
    else:
        with open(out_filename, "w") as outfile:
            pickle.dump(out_variables, outfile)
        print("Written " + out_filename)
        return out_variables


def filtering(datasets_in, filters, keep_missing=False):
    """
    Apply filtering (in time dimension) to entire dataset.
    
    Filters supplied in a list and then applied in order. So, if you want
    a daily, daytime average, you could do this:
    
    datasets_dictionary = filtering(datasets_dictionary, 
                                    ["daytime", "daily_median"])
    
    The first filter "daytime" selects data between 1000 and 1500 UTC,
    the second "daily_median" calculates the daily median. Obviously in this
    case, you need to do the first filter before the second.    
    """

    if type(filters) is not list:
        filters = [filters]

    datasets = datasets_in.copy()

    def ferry_loc(dataset, site,keep_missing=False):
        # Subset during daytime hours
        if site == 'GAUGE-FERRY':
                    
            lats = dataset.release_lat
            lons = dataset.release_lon
            ti = [i for i, h in enumerate(lons) if h > 0. and lats[i] > 54.4]
            ti2=np.arange(len(lons))
            ti3=np.delete(ti2,ti)
            
            if keep_missing:
                dataset_temp = dataset[dict(time = ti3)]   
                dataset_out = dataset_temp.reindex_like(dataset)
                return dataset_out
            else:
                return dataset[dict(time = ti3)]
        else:
            return dataset

    # Filter functions
    def daily_median(dataset, site,keep_missing=False):
        # Calculate daily median
        return dataset.resample("1D", "time", how = "median")
        
    def six_hr_mean(dataset, site,keep_missing=False):
        # Calculate daily median
        return dataset.resample("6H", "time", how = "mean")
    

    def daytime(dataset, site,keep_missing=False):
        # Subset during daytime hours
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h >= 11 and h <= 15]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def nighttime(dataset, site,keep_missing=False):
        # Subset during daytime hours
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h >= 23 or h <= 3]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def noon(dataset, site,keep_missing=False):
        # Subset during daytime hours
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h >= 12 and h < 13]
        
        if keep_missing:
            dataset_temp = dataset[dict(time = ti)]   
            dataset_out = dataset_temp.reindex_like(dataset)
            return dataset_out
        else:
            return dataset[dict(time = ti)] 
        


    def pblh_gt_500(dataset,site, keep_missing=False):
        # Subset for times when boundary layer height is > 500m
        ti = [i for i, pblh in enumerate(dataset.PBLH) if pblh > 500.]
        
        if keep_missing:
            mf_data_array = dataset.mf            
            dataset_temp = dataset.drop('mf')
            
            dataarray_temp = mf_data_array[dict(time = ti)]   
            #dataarray_temp2 = dataarray_temp.reindex_like(dataset)
            
            mf_ds = xray.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})
            
            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def pblh_gt_250(dataset,site, keep_missing=False):
        # Subset for times when boundary layer height is > 500m
        ti = [i for i, pblh in enumerate(dataset.PBLH) if pblh > 250.]
        
        if keep_missing:
            mf_data_array = dataset.mf            
            dataset_temp = dataset.drop('mf')
            
            dataarray_temp = mf_data_array[dict(time = ti)]   
            #dataarray_temp2 = dataarray_temp.reindex_like(dataset)
            
            mf_ds = xray.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})
            
            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
            
    def blxws_gt_5000(dataset,site, keep_missing=False):
        # Subset for times when boundary layer height is > 500m
        blxws = dataset.PBLH*dataset.wind_speed
        ti = [i for i, dum in enumerate(blxws) if dum > 5000.]
        
        if keep_missing:
            mf_data_array = dataset.mf            
            dataset_temp = dataset.drop('mf')
            
            dataarray_temp = mf_data_array[dict(time = ti)]   
            #dataarray_temp2 = dataarray_temp.reindex_like(dataset)
            
            mf_ds = xray.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})
            
            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]   
            
    def ferry_mf(dataset,site, keep_missing=False):
        # Subset for times when ferry mole fraction < 2400 ppb
        # Filter should only apply to ferry
        if site == 'GAUGE-FERRY':
            ti = [i for i, mf in enumerate(dataset.mf) if mf < 2400.]
            
            if keep_missing:
                mf_data_array = dataset.mf            
                dataset_temp = dataset.drop('mf')
                
                dataarray_temp = mf_data_array[dict(time = ti)]   
                #dataarray_temp2 = dataarray_temp.reindex_like(dataset)
                
                mf_ds = xray.Dataset({'mf': (['time'], dataarray_temp)}, 
                                      coords = {'time' : (dataarray_temp.coords['time'])})
                
                dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
                return dataset_out
            else:
                return dataset[dict(time = ti)]
        else:
            return dataset
            
    def ferry_fp_zero(dataset,site, keep_missing=False):
        # For some reason some ferry fps are zero - so filter them out
        # Filter should only apply to ferry
        if site == 'GAUGE-FERRY':
            fp_sum = dataset.fp.sum(["lat","lon"])
            ti = [i for i, mf in enumerate(fp_sum) if mf > 0.]
            
            if keep_missing:
                mf_data_array = dataset.mf            
                dataset_temp = dataset.drop('mf')
                
                dataarray_temp = mf_data_array[dict(time = ti)]   
                #dataarray_temp2 = dataarray_temp.reindex_like(dataset)
                
                mf_ds = xray.Dataset({'mf': (['time'], dataarray_temp)}, 
                                      coords = {'time' : (dataarray_temp.coords['time'])})
                
                dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
                return dataset_out
            else:
                return dataset[dict(time = ti)]
        else:
            return dataset
            
    def local_influence(dataset,site, keep_missing=False):
        
        #lr = dataset.local_ratio - np.min(dataset.local_ratio)
        #lr = dataset.local_ratio/np.percentile(dataset.local_ratio,5)
        lr = dataset.local_ratio
        pc = np.percentile(lr,5)*2
        ti = [i for i, local_ratio in enumerate(lr) if local_ratio <= pc]
        #ti = [i for i, local_ratio in enumerate(lr) if local_ratio < 0.04]
        #ti = [i for i, local_ratio in enumerate(lr) if local_ratio < 0.1]
        #ti = [i for i, local_ratio in enumerate(lr) if local_ratio < 0.05]
        if keep_missing is True: 
            mf_data_array = dataset.mf            
            dataset_temp = dataset.drop('mf')
            
            dataarray_temp = mf_data_array[dict(time = ti)]   
            #dataarray_temp2 = dataarray_temp.reindex_like(dataset)
            
            mf_ds = xray.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})
            
            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
       
        
    filtering_functions={"daily_median":daily_median,
                         "daytime":daytime,
                         "nighttime":nighttime,
                         "noon":noon,
                         "pblh_gt_500": pblh_gt_500,
                         "pblh_gt_250": pblh_gt_250,
                         "local_influence":local_influence,
                         "six_hr_mean":six_hr_mean,
                         "ferry_loc":ferry_loc,
                         "ferry_mf":ferry_mf,
                         "ferry_fp_zero":ferry_fp_zero,
                         "blxws_gt_5000":blxws_gt_5000}

    # Get list of sites
    sites = [key for key in datasets.keys() if key[0] != '.']
    
    # Do filtering
    for site in sites:
    
            for filt in filters:
                datasets[site] = filtering_functions[filt](datasets[site], site, keep_missing=keep_missing)

    return datasets


def baseline(y, y_time, y_site, baseline_error = 100., days_to_average = 5.):
    """
    Prepares an add-on to the sensitivity (H) matrix, that allows the baseline
    to be solved within the inversion.
    
    Specify how many days of data you want the baseline to be solved for at a time.
    If this number will leave one or two days at the end of the month, these days
    will be added to the previous bucket. For example, a days_to_average value of 5
    might leave a 6 day long bucket at the end of the month.
    
    Also specify the error on x for the inversion.
    """
    keys = np.unique(y_site)

    n = float(days_to_average)
    pos = np.zeros(len(y))
    for site in keys:
        val = np.max(pos)
        wh = np.where(y_site == site)
        ts = pd.Series(1, y_time[wh])
        fiveday = np.clip((ts.index.day-1) // n, 0, n)
        months = (ts.index.month - ts.index.month[0])
        pos[wh] = (val + 1) + fiveday + months*(n+1)
    
    HB = np.zeros((len(y), len(np.unique(pos))))
    xerror = np.zeros(len(np.unique(pos)))
    col = 0
    for i in range(len(np.unique(pos))):
        wh = np.where(pos == i+1)
        HB[wh, col] = 1
        xerror[col] = float(baseline_error)
        col += 1
                
    return HB, xerror


def gauss_inversion(data, prior, model, data_error, prior_error):
    """
    Simple gaussian inversion. Inputs are required to be arrays,
    already in the correct dimensions for the inversion:
    data : (m,)
    prior : (n,)
    model : (m,n)
    data_error : (m,m)
    prior_error : (n,n)
    """
    if type(data) is not np.matrix:
        data = np.matrix(data)
    if type(prior) is not np.matrix:
        prior = np.matrix(prior)
    if type(model) is not np.matrix:
        model = np.matrix(model)
    if type(data_error) is not np.matrix:
        data_error = np.matrix(data_error)
    if type(prior_error) is not np.matrix:
        prior_error= np.matrix(prior_error)

    y = data.T
    xa = prior.T
    H = model
    R1 = data_error.I
    P1 = prior_error.I
    
    P = ((H.T*R1)*H + P1).I
    x = xa + P*H.T*R1*(y - H*xa)
        
    return x, P


def prior_flux(species, domain, basis_case, av_date, emissions_name = None):
    """
    Calculates an area weighted flux for each basis region using the prior emissions.
    'av_date' is a date representative of the timeseries being looked at. It is used
    to find the most appropriate flux and basis files. Will only work for a flux that
    doesn't change during the peirod of the inversion.
    
    WARNING! This has a bug at the moment, only finds the nearest flux or basis
    file ON OR BEFORE the inputted av_date.
    """

    if emissions_name == None:
        emissions_name = species

    flux_data= flux(domain, emissions_name)
    basis_data = basis(domain, basis_case)
    
    print av_date
    
#    fi = np.argmin(np.abs(flux_data.time.values.to_pydatetime() - av_date))  
#    bi = np.argmin(np.abs(basis_data.time.values.to_pydatetime() - av_date)) 
#    
#    print fi, bi
#    
#    flux_timestamp = flux_data.time.values[fi]
#    basis_timestamp = basis_data.time.values[bi]
    
#    print flux_timestamp
#    print basis_timestamp
    
    flux_timestamp = pd.DatetimeIndex(flux_data.time.values).asof(av_date)
    basis_timestamp = pd.DatetimeIndex(basis_data.time.values).asof(av_date)
    
    flux0 = flux_data.flux.sel(time=flux_timestamp).values
    basis0 = basis_data.basis.sel(time=basis_timestamp).values
    
    basis_nos = range(int(np.min(basis0)), int(np.max(basis0))+1)
    
    area = areagrid(flux_data.lat.values, flux_data.lon.values)
            
    awflux = flux0*area
    basisflux = np.zeros(np.shape(basis_nos))
        
    for i in basis_nos:
        basisflux[i-1] = np.sum(awflux[basis0[:,:]==i])

    prior_x = unit_convert.mol2kg(basisflux,species)*(3600*24*365)
    
    return prior_x
    
    
def scaling_to_post_flux(prior_x, x, P):
    """
    For Gaussian inversion will convert outputs (emission scaling factors and associated error)
    to emission estimates using the prior flux calculated using the function 'prior_flux'.
    """

    if type(x) is not np.array:
        x = np.array(x)
    x2 = x.reshape(np.shape(prior_x))

    post_x = np.array(x2)*prior_x

    qmatrix = np.zeros((len(post_x), len(post_x)))
    for i in range(len(post_x)):
        qmatrix[:,i] = post_x[:]*post_x[i]

    if type(P) is not np.array:
        P = np.array(P)
    V = np.array(P)*qmatrix
    uncert = sum(sum(V))**0.5
    
    return post_x, uncert


class analytical_inversion:
    def __init__(self, out_var_file, species, domain, basis_case='voronoi',
                 emissions_name = None, prior_error=1., baseline_error = 100., baseline_days = 5.):
        """
        Using the output file from merge_sensitivity will calculate emissions estimates
        using a Gaussian analyical inversion.
        
        x_error is percentage error for the emissions scaling factor (1 = 100% error).
        
        If H_bc exists (the boundary conditions have been found using vmrs from a global model)
        then there is no need to specify baseline_days, however baseline_error is still needed.
        
        species_key is used if the species name gives more information than just the species
        (e.g. ch4-fossil-fire : species is 'ch4-fossil-fire', species_key is 'ch4')
        """
        
        y, y_error, y_site, y_time, H_bc, H = pickle.load(open(out_var_file))
        
        x0 = np.ones(len(H[0,:]))
        
#       Solve for baseline
        if H_bc is not None:
            H = np.append(H, H_bc, axis=1)
            xerror_bl = np.ones(len(H_bc[0,:]))*float(baseline_error)
        else:
            H_bc, xerror_bl = baseline(y, y_time, y_site, baseline_error = baseline_error, days_to_average = baseline_days)
            H = np.append(H, H_bc, axis = 1)
        
        prior_bl = np.dot(H_bc,np.ones(len(H_bc[0,:])))
        
#       Inversion
        xa = np.append(x0,np.ones(len(xerror_bl)))
        xerror = np.ones(len(x0))*float(prior_error)
        P = np.diagflat(np.append(xerror**2, xerror_bl**2))
        if y_error == None:
            y_error = y*0.1
        
        R = np.diagflat(y_error**2)
        
        x, P = gauss_inversion(y, xa, H, R, P)      
        
#       Find middle time in y_time to get appropriate flux and basis year in scaling_to_emissions
        av_date = y_time[len(y_time)/2]
                
#       Find real emissions values
        prior = prior_flux(species, domain, basis_case, av_date, emissions_name = emissions_name)
        posterior, uncertainty = scaling_to_post_flux(prior, x[:len(x0)], P[:len(x0),:len(x0)])   
        
#       Find baseline solution
        BL = H[:,len(x0):]*x[len(x0):] #H[:,len(H_bc[0,:]):]*x[len(H_bc[0,:]):]
    
        self.prior_scal = xa
        self.model = H
        self.obs = y
        self.time = y_time
        self.site_code = y_site
        self.obs_error = y_error
        self.post_scal = x
        self.prior_emi = prior
        self.post_emi = posterior
        self.uncert = uncertainty
        self.prior_bl = prior_bl
        self.post_bl = BL



#    
#    time, H0 = sensitivity_single_site(sites[0], species,
#                                       years=years, flux_years=flux_years,
#                                       domain=domain, basis_case=basis_case,
#                                       filt=filt)
#    H_shape=H0.shape
#    H=np.zeros((H_shape[0]*len(sites), H_shape[1]))
#    H[0:H_shape[0],:]=H0
#    
#    H_site_name = []
#    H_site_name += [sites[0] for dummy in range(H_shape[0])]
#    
#    H_time = time
#    print("... done " + sites[0])
#    
#    if len(sites) > 1:
#        for si, site in enumerate(sites[1:]):
#            time_site, H_site = sensitivity_single_site(site, species,
#                                       years=years, flux_years=flux_years,
#                                       domain=domain, basis_case=basis_case,
#                                       filt=filt)
#            H[H_shape[0]*(si+1):H_shape[0]*(si+2), :] = H_site
#            H_site_name += [site for dummy in range(H_shape[0])]
#            H_time += time_site
#            print("... done " + site)
#            # Add capability to incorporate footprints of different times
#
#    return H_time, H_site_name, H


class plot_map_setup:
    def __init__(self, fp_data, 
                 lon_range = None, lat_range = None,
                 bottom_left = False,
                 map_resolution = "l"):

        if lon_range is None:
            lon_range = (min(fp_data.lon.values),
                         max(fp_data.lon.values))
        if lat_range is None:
            lat_range = (min(fp_data.lat.values),
                         max(fp_data.lat.values))
        
        m = Basemap(projection='gall',
            llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
            resolution = map_resolution)

        if bottom_left == False:
            lons, lats = np.meshgrid(fp_data.lon.values,
                                     fp_data.lat.values)
        else:
            dlon = fp_data.lon.values[1] - fp_data.lon.values[0]
            dlat = fp_data.lat.values[1] - fp_data.lat.values[0]            
            lons, lats = np.meshgrid(fp_data.lon.values - dlon,
                                     fp_data.lat.values - dlat)
        
        x, y = m(lons, lats)
        
        self.x = x
        self.y = y
        self.m = m


def plot_default_colors(site):
    
    cmap = {"GAUGE-FERRY": plt.cm.Blues,
            "GAUGE-FAAM": plt.cm.Reds}
            
    if site in cmap.keys():
        return cmap[site]
    else:
        return (plt.cm.BuPu + plt.cm.Reds)/2.


def plot_map_zoom(fp_data):
    
    sites = [key for key in fp_data.keys() if key[0] != '.']
    
    dlon = max(fp_data[sites[0]].lon.values) - \
            min(fp_data[sites[0]].lon.values)
    lon_range = [min(fp_data[sites[0]].lon.values) + 0.6*dlon,
                 max(fp_data[sites[0]].lon.values) - 0.25*dlon]
    dlat = max(fp_data[sites[0]].lat.values) - \
            min(fp_data[sites[0]].lat.values)
    lat_range = [min(fp_data[sites[0]].lat.values) + 0.53*dlat,
                 max(fp_data[sites[0]].lat.values) - 0.25*dlat]

    return lat_range, lon_range


def plot(fp_data, date, out_filename=None, 
         lon_range=None, lat_range=None, log_range = [5., 9.],
         map_data = None, zoom = False,
         map_resolution = "l", 
         map_background = "countryborders",
         colormap = plt.cm.YlGnBu,
         tolerance = None,
         interpolate = False,
         dpi = None):
    """
    Plot footprint using pcolormesh.
    
    Arguments:
    fp_data: Dictionary of xray datasets containing footprints
        and other variables
    date: Almost any time format should work (datetime object, string, etc).
        Footprints from all sites in dictionary that have time indices at that time
        will be plotted.
        
    Keywords:
    lon_range: list of min and max longitudes [min, max]
    lat_range: list of min and max latitudes [min, max]
    log_range: list of min and max LOG10(footprints) for 
        color scale
    map_data: contains plot_map_setup class (useful if animating, 
        so that entire map does not need to be re-calculated each time)
    zoom: shortcut to zoom in to subset of domain (needs work)
    map_resolution: resolution of map
    map_background: choice of "countryborders", "shadedrelief"
    colormap: color map to use for contours
    tolerance: Xray doc "Maximum distance between original and new
        labels for inexact matches."
    """
    
    def fp_nearest(fp, tolerance = None):
        return fp.reindex_like( \
                            xray.Dataset(coords = {"time": [date]}),
                            method = "nearest",
                            tolerance = tolerance)

    
    # Looks for nearest time point aviable in footprint   
    date = np.datetime64(convert.reftime(date))

    # Get sites
    sites = [key for key in fp_data.keys() if key[0] != '.']

    # Zoom in. Assumes release point is to the East of centre
    if zoom:
        lat_range, lon_range = plot_map_zoom(fp_data)
    
    # Get map data
    if map_data is None:
        map_data = plot_map_setup(fp_data[sites[0]],
                                  lon_range = lon_range,
                                  lat_range = lat_range,
                                  bottom_left=True,
                                  map_resolution = map_resolution)

    # Open plot
    fig = plt.figure(figsize=(8,8))
    fig.add_axes([0.1,0.1,0.8,0.8])

    if map_background == "shadedrelief":
        map_data.m.shadedrelief()
    elif map_background == "countryborders":
        map_data.m.drawcoastlines()
        map_data.m.fillcontinents(color='grey',lake_color=None, alpha = 0.2)
        map_data.m.drawcountries()

    #Calculate color levels
    cmap = colormap
    rp_color = {"SURFACE": "black",
                "SHIP": "purple",
                "AIRCRAFT": "red",
                "SATELLITE": "green"}
            
    levels = MaxNLocator(nbins=256).tick_values(log_range[0], log_range[1])
    #levels=np.arange(log_range[0],log_range[1], 0.2)
    norm = BoundaryNorm(levels,
                        ncolors=cmap.N,
                        clip=True)

    # Create dictionaries and arrays    
    release_lon = {}
    release_lat = {}

    data = {}

    data = np.zeros(np.shape(
            fp_data[sites[0]].fp[dict(time = [0])].values.squeeze()))
    
    time = None

    # Generate plot data
    for site in sites:
    
        time = fp_data[site].time.values.squeeze().copy()
        
        if tolerance is None:
            tol = np.median(time[1:] - time[:-1])
        
        if interpolate:

            ti = bisect.bisect(time, date)
            
            if (ti > 0) and (ti < len(time)):

                if np.float(date - time[ti-1]) < np.float(tol):

                    dt = np.float(date - time[ti-1])/np.float(time[ti] - time[ti-1])
                    fp_ti_0 = fp_data[site][dict(time = ti-1)].fp.values.squeeze()
                    fp_ti_1 = fp_data[site][dict(time = ti)].fp.values.squeeze()
                    fp_ti = fp_ti_0 + (fp_ti_1 - fp_ti_0)*dt
                    
                    data += np.nan_to_num(fp_ti)
    
                    # Store release location to overplot later
                    if "release_lat" in fp_data[site].keys():
                        lon_0 = fp_data[site][dict(time = ti-1)].release_lon.values.squeeze()
                        lon_1 = fp_data[site][dict(time = ti)].release_lon.values.squeeze()                
                        lat_0 = fp_data[site][dict(time = ti-1)].release_lat.values.squeeze()
                        lat_1 = fp_data[site][dict(time = ti)].release_lat.values.squeeze()
                        release_lon[site] = lon_0 + (lon_1 - lon_0)*dt
                        release_lat[site] = lat_0 + (lat_1 - lat_0)*dt
                    
            else:

                fp_data_ti = fp_nearest(fp_data[site], tolerance = tolerance)
                data += np.nan_to_num(fp_data_ti.fp.values.squeeze())
                # Store release location to overplot later
                if "release_lat" in dir(fp_data_ti):
                    release_lon[site] = fp_data_ti.release_lon.values
                    release_lat[site] = fp_data_ti.release_lat.values

        else:

            fp_data_ti = fp_nearest(fp_data[site], tolerance = tolerance)
            data += np.nan_to_num(fp_data_ti.fp.values.squeeze())

            # Store release location to overplot later
            if "release_lat" in dir(fp_data_ti):
                release_lon[site] = fp_data_ti.release_lon.values
                release_lat[site] = fp_data_ti.release_lat.values

    #Set very small elements to zero
    data = np.log10(data)
    data[np.where(data <  log_range[0])]=np.nan
    
    #Plot footprint
    cs = map_data.m.pcolormesh(map_data.x, map_data.y,
                               np.ma.masked_where(np.isnan(data), data),
                               cmap = cmap, norm = norm)

    # over-plot release location
    if len(release_lon) > 0:
        for site in sites:
            if site in release_lon:
                rplons, rplats = np.meshgrid(release_lon[site],
                                             release_lat[site])
                rpx, rpy = map_data.m(rplons, rplats)
                if "platform" in site_info[site]:
                    color = rp_color[site_info[site]["platform"].upper()]
                else:
                    color = rp_color["SURFACE"]
                rp = map_data.m.scatter(rpx, rpy, 100, color = color)
    
    plt.title(str(pd.to_datetime(str(date))), fontsize=20)

    cb = map_data.m.colorbar(cs, location='bottom', pad="5%")
    
    tick_locator = ticker.MaxNLocator(nbins=10)
    cb.locator = tick_locator
    cb.update_ticks()
 
    cb.set_label('log$_{10}$( (nmol/mol) / (mol/m$^2$/s) )', 
                 fontsize=15)
    cb.ax.tick_params(labelsize=13) 
    
    if out_filename is not None:
        plt.savefig(out_filename, dpi = dpi)
        plt.close()
    else:
        plt.show()


def time_unique(fp_data, time_regular = False):
    
    sites = [key for key in fp_data.keys() if key[0] != '.']
    
    time_array = fp_data[sites[0]].time
    time_array.name = "times"
    time = time_array.to_dataset()
    if len(sites) > 1:
        for site in sites[1:]:
            time.merge(fp_data[site].time.to_dataset(), inplace = True)
    
    time = time.time.to_pandas().index.to_pydatetime()
    if time_regular is not False:
        time = pd.date_range(min(time), max(time), freq = time_regular)
    
    return time


def plot3d(fp_data, date, out_filename=None, 
           log_range = [5., 9.]):

    """date as string "d/m/y H:M" or datetime object 
    datetime.datetime(yyyy,mm,dd,hh,mm)
    """
    
    #looks for nearest time point aviable in footprint   
    if isinstance(date, str):
        date=dt.datetime.strptime(date, '%d/%m/%Y %H:%M')

    time_index = bisect.bisect_left(fp_data.time, date)

    data = np.log10(fp_data.fp.values[:,:,time_index].squeeze())
    lon_range = (fp_data.lon.min().values, fp_data.lon.max().values)
    lat_range = (fp_data.lat.min().values, fp_data.lat.max().values)

    #Set very small elements to zero
    data[np.where(data <  log_range[0])]=np.nan

    fig = plt.figure(figsize=(16,12))

    fig.text(0.1, 0.2, str(date), fontsize = 20)
    
    ax = Axes3D(fig)
    
    ax.set_ylim(lat_range)
    ax.set_xlim(lon_range)
    ax.set_zlim((min(fp_data.height), max(fp_data.height)))

    fpX, fpY = np.meshgrid(fp_data.lon, fp_data.lat)
    
    levels = np.arange(log_range[0], log_range[1], 0.05)

    plfp = ax.contourf(fpX, fpY, data, levels, offset = 0.)
    plnX, plnY = np.meshgrid(fp_data.lon.values.squeeze(), fp_data.height.values.squeeze())
    plwX, plwY = np.meshgrid(fp_data.lat.values.squeeze(), fp_data.height.values.squeeze())
    pllevs = np.arange(0., 0.0031, 0.0001)
    
    plnvals = fp_data.particle_locations_n.values[:,:,time_index].squeeze()
    plnvals[np.where(plnvals == 0.)]=np.nan
    plwvals = fp_data.particle_locations_w.values[:,:,time_index].squeeze()
    plwvals[np.where(plwvals == 0.)]=np.nan
    plpln = ax.contourf(plnX, plnvals, plnY,
                        zdir = 'y', offset = max(fp_data.lat.values), levels = pllevs)
    plplw = ax.contourf(plwvals, plwX, plwY,
                        zdir = 'x', offset = min(fp_data.lon.values), levels = pllevs)    
    ax.view_init(50)

    cb = plt.colorbar(plfp, location='bottom', shrink=0.8)
    tick_locator = ticker.MaxNLocator(nbins=7)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_label('log$_{10}$( (mol/mol) / (mol/m$^2$/s))', 
             fontsize=15)
    cb.ax.tick_params(labelsize=13) 
    
    if out_filename is not None:
        plt.savefig(out_filename, dpi = 300)
        plt.close()
    else:
        plt.show()


def animate(fp_data, output_directory, 
            lon_range = None, lat_range=None, zoom = False,
            log_range = [5., 9.],
            overwrite=True, file_label = 'fp',
            framerate=10, delete_png=False,
            video_os="mac", ffmpeg_only = False,
            plot_function = "plot", time_regular = "1H",
            frame_limit = None,
            dpi = 150, interpolate = True):
    
    sites = [key for key in fp_data.keys() if key[0] != '.']
    
    if ffmpeg_only is False:

        # Set up map        
        if zoom:
            lat_range, lon_range = plot_map_zoom(fp_data)
        
        map_data = plot_map_setup(fp_data[sites[0]], 
                                  lon_range = lon_range, 
                                  lat_range = lat_range, bottom_left = True)
        
        # Find unique times
        times = time_unique(fp_data,
                            time_regular = time_regular)
        
        if frame_limit:
            times = times[:min([frame_limit, len(times)])]

        # Start progress bar
        pbar = ProgressBar(maxval=len(times)).start()

        # Plot each timestep
        for ti, t in enumerate(times):
            
            fname=os.path.join(output_directory, 
                               file_label + '_' + str(ti).zfill(5) + '.png')
            
            if len(glob.glob(fname)) == 0 or overwrite == True:
                if plot_function == "plot":
                    plot(fp_data, t, out_filename = fname, 
                         lon_range = lon_range, lat_range= lat_range,
                         log_range = log_range, map_data = map_data,
                         interpolate = interpolate,
                         dpi = dpi)
                elif plot_function == "plot3d":
                    plot3d(fp_data[sites[0]], t, out_filename = fname,
                           log_range = log_range)
                     
            pbar.update(ti)
        pbar.finish()
    
    print("")
    print("... running ffmpeg")

    if video_os.lower() == "mac":
        ffmpeg_status = subprocess.call("ffmpeg -r " + str(framerate) + \
            " -i '" + os.path.join(output_directory, file_label) + "_%05d.png' " + \
            "-f mp4 -vcodec libx264 " + \
            "-pix_fmt yuv420p -intra -qscale 0 -y " + \
            os.path.join(output_directory, file_label) + ".mp4", shell=True)
    elif video_os.lower() == "pc":
#        os.remove(os.path.join(output_directory, file_label) + ".wmv")
        ffmpeg_status = subprocess.call("ffmpeg -r " + str(framerate) + \
            " -i '" + os.path.join(output_directory, file_label) + "_%05d.png' " + \
            " -b 5000k -f asf -vcodec wmv2 -acodec wmav2 " + \
            os.path.join(output_directory, file_label) + ".wmv", shell=True)
    else:
        print("ERROR: video_os must be mac or pc")
        return None
    
    print("... done with status " + str(ffmpeg_status))

    if delete_png:
        filelist = glob.glob(os.path.join(output_directory, "*.png"))
        for f in filelist:
            os.remove(f)


class get_country:
  def __init__(self, domain, ocean=False, al=False, uk_split=False):

        if ocean is False:

            countryDirectory=data_path +'NAME/countries/'
            filename=glob.glob(countryDirectory + \
                 "/" + "country_" \
                 + domain + ".nc")
             
        else:
            if uk_split is True:
                countryDirectory=data_path +'NAME/countries/'
                filename=glob.glob(countryDirectory + \
                     "/" + "country_"\
                     + domain + "_uk_split" + ".nc")
            else:
                if al is False:
                    countryDirectory=data_path +'NAME/countries/'
                    filename=glob.glob(countryDirectory + \
                         "/" + "country_ocean_"\
                         + domain + ".nc")
                else:
                    countryDirectory=data_path +'NAME/countries/'
                    filename=glob.glob(countryDirectory + \
                         "/" + "al_countries_"\
                         + domain + ".nc")
        
        f = nc.Dataset(filename[0], 'r')
    
        lon = f.variables['lon'][:]
        lat = f.variables['lat'][:]
    
        #Get country indices and names
        country = f.variables['country'][:, :]
        if al is True:
            name_temp = f.variables['name'][:]  
            f.close()
            name=np.asarray(name_temp)
        elif uk_split is True:
            name_temp = f.variables['name'][:]  
            f.close()
            name=np.asarray(name_temp)
            
        else:
        #if al is False:
            name_temp = f.variables['name'][:,:]
            f.close()
    
            name=[]
            for ii in range(len(name_temp[:,0])):
                name.append(''.join(name_temp[ii,:]))
            name=np.asarray(name)
#        else:
#            name_temp = f.variables['name'][:]  
#            f.close()
#            name=np.asarray(name_temp)
    
    
        self.lon = lon
        self.lat = lat
        self.lonmax = np.max(lon)
        self.lonmin = np.min(lon)
        self.latmax = np.max(lat)
        self.latmin = np.min(lat)
        self.country = np.asarray(country)
        self.name = name 
