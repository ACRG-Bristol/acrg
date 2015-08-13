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
from acrg_grid import areagrid
import xray
from os.path import split, realpath
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

fp_directory = data_path + 'NAME/fp/'
flux_directory = data_path +'NAME/emissions/'
basis_directory = data_path + 'NAME/basis_functions/'
bc_directory = data_path +'NAME/bc/'
bc_basis_directory = data_path +'NAME/bc_basis_functions/'

# Get acrg_site_info file
acrg_path=split(realpath(__file__))
with open(acrg_path[0] + "/acrg_site_info.json") as f:
    site_info=json.load(f)

def filenames(site, domain, start, end, height = None, flux=None, basis=None):
    """
    Output a list of available footprint file names,
    for given site, domain, directory and date range.

    Doesn't work for flux or basis functions yet.
    """
    if flux is None and basis is None:
        baseDirectory = fp_directory
    else:
        if flux is not None:
            baseDirectory = flux_directory
        if basis is not None:
            baseDirectory = basis_directory
    
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
            site + "*" + height + "*" + domain + "*" + ym + "*.nc")
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
        domain="EUROPE", height = None, species = None):
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
    
    If the SPECIES keyword is given, fluxes will also be loaded and
    merged into the dataset.    
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
        if species is not None:

            flux_ds = flux(domain, species)

            if flux_ds is not None:
                fp = combine_datasets(fp, flux_ds)

            bc_ds = boundary_conditions(domain, species)
            if bc_ds is not None:
                fp = combine_datasets(fp, bc_ds)

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
                    basis_case + "*.nc"))
    if len(files) == 0:
        print("Can't find boundary condition basis functions: " + domain + " " + basis_case)
        return None
    
    basis_ds = read_netcdfs(files)

    return basis_ds


def combine_datasets(dsa, dsb, method = "nearest"):
    """
    Merge two datasets. Assumes that you want to 
    re-index to the FIRST dataset.
    
    Example:
    
    ds = combine_datasets(dsa, dsb)

    ds will have the index of dsa    
    """
    
    return dsa.merge(dsb.reindex_like(dsa, method))


def timeseries(ds):
    """
    Compute flux * footprint time series.
    All that is required is that you input an xray
    dataset with both the flux and footprint fields present    
    
    Example:
    
        ts = timeseries(dataset)
    
    There are almost certainly much more efficient ways of doing this.
    """

    return (ds.fp*ds.flux).sum(["lat", "lon"])


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

    
def footprints_data_merge(data, domain = "EUROPE", species = "CH4",
                          calc_timeseries = True, calc_bc = True,
                          average = None, site_modifier = {}, height = {}):
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

    # Check if species is defined in data dictionary
    if ".species" in data.keys():
        species = data[".species"]

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
        
        # Get footprints
        site_fp = footprints(site_modifier_fp, start = start, end = end,
                             domain = domain,
                             species = [species if calc_timeseries == True or \
                                         calc_bc == True \
                                         else None][0], \
                             height = height[site])
        
        if site_fp is not None:
            
            # Merge datasets
            site_ds = combine_datasets(site_ds, site_fp, method = "nearest")
            
            # If units are specified, multiply by scaling factor
            if ".units" in attributes:
                site_ds.update({'fp' : (site_ds.fp.dims, site_ds.fp / data[".units"])})
                if calc_bc:
                    for key in site_ds.keys():
                        if "vmr" in key:
                            site_ds.update({key :
                                            (site_ds[key].dims, site_ds[key] / \
                                            data[".units"])})

               

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


def fp_sensitivity(fp_and_data, domain = 'EUROPE', basis_case = 'voronoi'):
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

        site_bf = combine_datasets(fp_and_data[site]["fp", "flux", "mf_mod"],
                                   basis_func)

        basis_scale = xray.Dataset({'basis_scale': (['lat','lon','time'],
                                                    np.zeros(np.shape(site_bf.basis)))},
                                   coords = site_bf.coords)
        site_bf = site_bf.merge(basis_scale)

        H = np.zeros((int(np.max(site_bf.basis)),len(site_bf.mf_mod)))
        
        fp=site_bf.fp.values
        flux=site_bf.flux.values
        H_all=fp*flux 
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
        
        if basis_case == 'transd':
            
            sub_fp_temp = site_bf.fp.sel(lon=slice(min(site_bf.sub_lon),max(site_bf.sub_lon)), 
                                    lat=slice(min(site_bf.sub_lat),max(site_bf.sub_lat)))   
            sub_fp = xray.Dataset({'sub_fp': (['sub_lat','sub_lon','time'], sub_fp_temp)},
                               coords = {'sub_lat': (site_bf.coords['sub_lat']),
                                         'sub_lon': (site_bf.coords['sub_lon']),
                                'time' : (fp_and_data[site].coords['time'])})
            
            fp_and_data[site] = fp_and_data[site].merge(sub_fp)
    
    return fp_and_data


def bc_sensitivity(fp_and_data, domain = 'EUROPE', basis_case = 'NESW'):
    
    sites = [key for key in fp_and_data.keys() if key[0] != '.']
#    attributes = [key for key in fp_and_data.keys() if key[0] == '.']
    basis_func = basis_boundary_conditions(domain = domain,
                                           basis_case = basis_case)
    
    for site in sites:

        # stitch together the particle locations, vmrs at domain edges and
        #boundary condition basis functions
        DS = combine_datasets(fp_and_data[site]["particle_locations_n",
                                                "particle_locations_e",
                                                "particle_locations_s",
                                                "particle_locations_w",
                                                "vmr_n",
                                                "vmr_e",
                                                "vmr_s",
                                                "vmr_w",
                                                "bc"],
                                                basis_func)

        part_loc = np.hstack([DS.particle_locations_n,
                                DS.particle_locations_e,
                                DS.particle_locations_s,
                                DS.particle_locations_w])
        

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


def filtering(datasets_in, filters):
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

    # Filter functions
    def daily_median(dataset):
        # Calculate daily median
        return dataset.resample("1D", "time", how = "median")
    
    def daytime(dataset):
        # Subset during daytime hours
        hours = dataset.time.to_pandas().index.hour
        ti = [i for i, h in enumerate(hours) if h >= 10 and h <= 15]
        return dataset[dict(time = ti)]

    def pblh_gt_500(dataset):
        # Subset for times when boundary layer height is > 500m
        ti = [i for i, pblh in enumerate(dataset.PBLH) if pblh > 500.]
        return dataset[dict(time = ti)]
        
    filtering_functions={"daily_median":daily_median,
                         "daytime":daytime,
                         "pblh_gt_500": pblh_gt_500}

    # Get list of sites
    sites = [key for key in datasets.keys() if key[0] != '.']
    
    # Do filtering
    for site in sites:
        for filt in filters:
            datasets[site] = filtering_functions[filt](datasets[site])

    return datasets


def baseline(y, y_time, y_site, x_error = 10000, days_to_average = 5):
    
    keys = np.unique(y_site)

    n = days_to_average
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
        xerror[col] = x_error
        col += 1
                
    return HB, xerror


class analytical_inversion:
    def __init__(self, obs, species, years=[2012], flux_years=None,
                domain="small", basis_case='voronoi', filt=None,
                species_key = None, baseline_days = 5, alt_fp_filename = None):
    
        y_time, y_site, y, H = sensitivity(obs, species, years, flux_years=flux_years,
                domain=domain, basis_case=basis_case, filt=filt, alt_fp_filename = alt_fp_filename)
        
        if species_key == None:
            species_key = species
            
        acrg_path=os.path.split(os.path.realpath(__file__))
        with open(acrg_path[0] + "/acrg_species_info.json") as f:
            species_info=json.load(f)
            
            
        if type(species_key) is not str:
            species_key = str(species_key) 
        species_key = agage.synonyms(species_key, species_info)
    
        mol_mass = float(species_info[species_key]['mol_mass'])
        u = species_info[species_key]['units']
        
        units = {"ppm" : 1e6,
                 "ppb" : 1e9,
                 "ppt" : 1e12,
                 "ppq" : 1e15}
        
        H0 = H*units[u]
        x0 = np.ones(len(H[0,:]))
        
#       Solve for baseline    
        HB, xerror = baseline(y, y_time, y_site,obs, days_to_average = baseline_days)
        H = np.append(H0, HB, axis = 1)
        
#       Inversion
        xa = np.append(x0,np.zeros(len(xerror)))
        P = np.diagflat(np.append(x0**2, xerror))
        R = np.diagflat((y*0.1)**2)
        
        H = np.matrix(H)
        xa = np.matrix(xa).T
        y = np.matrix(y)
        P1 = np.matrix(P).I
        R1= np.matrix(R).I
            
        P = ((H.T*R1)*H + P1).I
        x = xa + P*H.T*R1*(y - H*xa)
        
#       Find real emissions values
        flux_data=flux(species, flux_years, domain=domain)
        basis_data = basis_function(basis_case, years=years, domain=domain)
        area = areagrid(flux_data.lat, flux_data.lon)
            
        awflux = flux_data.flux[:,:,0]*area
        fl = np.zeros(np.shape(x))
        
        for i in range(int(np.min(basis_data.basis)), int(np.max(basis_data.basis))+1):
            fl[i-1] = np.sum(awflux[basis_data.basis[:,:,0]==i])
        
        prior = sum(fl*(3600*24*365)*mol_mass*1e-3)
        E = np.array(x)*fl*(3600*24*365)*mol_mass*1e-3
        E_tot = sum(E)

        qmatrix = np.zeros((len(E), len(E)))
        for i in range(len(E)):
            qmatrix[:,i] = E[:,0]*E[i]

        V = np.array(P)*qmatrix
        sigma = sum(sum(V))**0.5
        
#       Find baseline solution
        BL = H[:,len(H0[0,:]):]*x[len(H0[0,:]):]
    
        self.prior_scal = xa
        self.model = H
        self.obs = y
        self.post_scal = x
        self.prior_emi = prior
        self.post_emi = float(E_tot)
        self.post_emi_allregions = E
        self.uncert = sigma
        self.baseline = BL

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
         colormap = plt.cm.YlGnBu):
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
    """
    
    # Looks for nearest time point aviable in footprint   
    date = convert.reftime(date)

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
                                  bottom_left = True,
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
    rp_color = {"SURFACE": "blue",
                "SHIP": "purple",
                "AIRCRAFT": "red",
                "SATELLITE": "green"}
            
    levels = MaxNLocator(nbins=100).tick_values(log_range[0], log_range[1])

    norm = BoundaryNorm(levels,
                        ncolors=cmap.N,
                        clip=True)

    # Create dictionaries and arrays    
    release_lon = {}
    release_lat = {}

    data = {}

    data = np.zeros(np.shape(
            fp_data[sites[0]].fp[dict(time = [0])].values.squeeze()))

    # Generate plot data
    for site in sites:
    
        tdelta = fp_data[site].time - date
        if np.min(np.abs(tdelta)) < 2*3600.*1e9:

            fp_data_ti = fp_data[site].reindex_like( \
                            xray.Dataset(coords = {"time": [date]}),
                            method = "nearest")

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
                rp = map_data.m.scatter(rpx, rpy, 40, color = color)
    
    plt.title(str(pd.to_datetime(str(date))), fontsize=20)

    cb = map_data.m.colorbar(cs, location='bottom', pad="5%")
    
    tick_locator = ticker.MaxNLocator(nbins=7)
    cb.locator = tick_locator
    cb.update_ticks()
 
    cb.set_label('log$_{10}$( (nmol/mol) / (mol/m$^2$/s) )', 
                 fontsize=15)
    cb.ax.tick_params(labelsize=13) 
    
    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()


def time_unique(fp_data):
    
    sites = [key for key in fp_data.keys() if key[0] != '.']
    
    time = fp_data[sites[0]].time.to_dataset()
    if len(sites) > 1:
        for site in sites[1:]:
            print(site)
            time.merge(fp_data[site].time.to_dataset(), inplace = True)
    
    return time


def plot3d(fp_data, date, out_filename=None, 
         cutoff = -3.5):

    """date as string "d/m/y H:M" or datetime object 
    datetime.datetime(yyyy,mm,dd,hh,mm)
    """
    
    #looks for nearest time point aviable in footprint   
    if isinstance(date, str):
        date=dt.datetime.strptime(date, '%d/%m/%Y %H:%M')

    time_index = bisect.bisect_left(fp_data.time, date)

    data = np.log10(fp_data.fp[:,:,time_index])
    lon_range = (fp_data.lonmin, fp_data.lonmax)
    lat_range = (fp_data.latmin, fp_data.latmax)

    #Set very small elements to zero
    data[np.where(data <  cutoff)]=np.nan

    fig = plt.figure(figsize=(16,12))

    fig.text(0.1, 0.2, str(date), fontsize = 20)
    
    ax = Axes3D(fig)
    
    ax.set_ylim(lat_range)
    ax.set_xlim(lon_range)
    ax.set_zlim((min(fp_data.particle_height), max(fp_data.particle_height)))

    fpX, fpY = np.meshgrid(fp_data.lon, fp_data.lat)
    
    levels = np.arange(cutoff, 0., 0.05)

    plfp = ax.contourf(fpX, fpY, data, levels, offset = 0.)
    plnX, plnY = np.meshgrid(fp_data.lon, fp_data.particle_height)
    plwX, plwY = np.meshgrid(fp_data.lat, fp_data.particle_height)
    pllevs = np.arange(0., 0.0031, 0.0001)
    
    plnvals = fp_data.particle_locations["N"][:,:,time_index]
    plnvals[np.where(plnvals == 0.)]=np.nan
    plwvals = fp_data.particle_locations["W"][:,:,time_index]
    plwvals[np.where(plwvals == 0.)]=np.nan
    plpln = ax.contourf(plnX, plnvals, plnY,
                        zdir = 'y', offset = max(fp_data.lat), levels = pllevs)
    plplw = ax.contourf(plwvals, plwX, plwY,
                        zdir = 'x', offset = min(fp_data.lon), levels = pllevs)    
    ax.view_init(50)

    cb = plt.colorbar(plfp, location='bottom', shrink=0.8)
    tick_locator = ticker.MaxNLocator(nbins=7)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.set_label('log$_{10}$( (mol/mol) / (mol/m$^2$/s))', 
             fontsize=15)
    cb.ax.tick_params(labelsize=13) 
    
    plt.tight_layout()

    if out_filename is not None:
        plt.savefig(out_filename, dpi = 600)
        plt.close()
    else:
        plt.show()


def animate(fp_data, output_directory, 
            lon_range = None, lat_range=None, zoom = False,
            cutoff = -3.,
            overwrite=True, file_label = 'fp', 
            framerate=10, delete_png=False,
            video_os="mac", ffmpeg_only = False):
    
    sites = [key for key in fp_data.keys() if key[0] != '.']
    
    if ffmpeg_only is False:

        # Set up map        
        if zoom:
            lat_range, lon_range = plot_map_zoom(fp_data)
            
        map_data = plot_map_setup(fp_data[sites[0]], 
                                  lon_range = lon_range, 
                                  lat_range= lat_range, bottom_left = True)
        
        # Find unique times
        times = time_unique(fp_data)

        # Start progress bar
        pbar=ProgressBar(maxval=len(times.time.values)).start()

        # Plot each timestep
        for ti, t in enumerate(times.time.values):
            
            fname=os.path.join(output_directory, 
                               file_label + '_' + str(ti).zfill(5) + '.png')
                               
            if len(glob.glob(fname)) == 0 or overwrite == True:            
                plot_scatter(fp_data, t, out_filename = fname, 
                     lon_range = lon_range, lat_range= lat_range,
                     cutoff=cutoff, map_data = map_data)
                     
            pbar.update(ti)
        pbar.finish()
    
    print ""
    print "... running ffmpeg"

    if video_os.lower() == "mac":
        ffmpeg_status = subprocess.call("ffmpeg -r " + str(framerate) + \
            " -i '" + os.path.join(output_directory, file_label) + "_%05d.png' " + \
            "-f mp4 -vcodec libx264 " + \
            "-pix_fmt yuv420p -intra -qscale 0 -y " + \
            os.path.join(output_directory, file_label) + ".mp4", shell=True)
    elif video_os.lower() == "pc":
        ffmpeg_status = subprocess.call("ffmpeg -r " + str(framerate) + \
            " -i '" + os.path.join(output_directory, file_label) + "_%05d.png' " + \
            " -b 5000k -f asf -vcodec wmv2 -acodec wmav2 " + \
            os.path.join(output_directory, file_label) + ".wmv", shell=True)
    else:
        print("ERROR: video_os must be mac or pc")
        
    print(ffmpeg_status)
    
    if delete_png:
        filelist = glob.glob(os.path.join(output_directory, "*.png"))
        for f in filelist:
            os.remove(f)


class get_country:
  def __init__(self, domain, ocean=False):

        if ocean is False:

            countryDirectory=data_path +'NAME/countries/'
            filename=glob.glob(countryDirectory + \
                 "/" + "country_" \
                 + domain + ".nc")
             
        else:
            countryDirectory=data_path +'NAME/countries/'
            filename=glob.glob(countryDirectory + \
                 "/" + "country_ocean_"\
                 + domain + ".nc")
        
        f = nc.Dataset(filename[0], 'r')
    
        lon = f.variables['lon'][:]
        lat = f.variables['lat'][:]
    
        #Get country indices and names
        country = f.variables['country'][:, :]
        name_temp = f.variables['name'][:,:]
        f.close()
    
        name=[]
        for ii in range(len(name_temp[:,0])):
            name.append(''.join(name_temp[ii,:]))
            
        name=np.asarray(name)
    
    
        self.lon = lon
        self.lat = lat
        self.lonmax = np.max(lon)
        self.lonmin = np.min(lon)
        self.latmax = np.max(lat)
        self.latmin = np.min(lat)
        self.country = np.asarray(country)
        self.name = name 
