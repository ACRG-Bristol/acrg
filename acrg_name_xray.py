# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 10:45:51 2014

Script to create footprint maps using netCDF processed NAME output.

site:
'MHD' Mace Head
'RGL' Ridgehill
'TAC' Tacolneston
'BSD' Bilsdale
'HFD' Heathfield

date in datetime format:
import datetime as dt
yourdate = dt.datetime(year,month,day,hour)
Footprint data available 2 hourly from midnight.

fpread:
adds footprints together if more than one.
r = fpread('MHD',2012)
r.fp 
array([[[  0.00000000e+00,   5.02990042e-06,   1.00698571e-05, ...,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   1.50677261e-05,   0.00000000e+00, ...,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   5.50260556e-05,   0.00000000e+00, ...,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        ..., 
        [  0.00000000e+00,   6.43811800e-05,   0.00000000e+00, ...,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   2.29741690e-05,   0.00000000e+00, ...,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
        [  0.00000000e+00,   7.58096758e-06,   0.00000000e+00, ...,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00]],
       ..., 
       [[  3.01224500e-05,   7.02261314e-05,   2.00294908e-05, ...,
           5.01100658e-05,   0.00000000e+00,   3.01291821e-05],
        [  1.00021070e-05,   1.99628830e-05,   4.00702156e-05, ...,
           5.99743616e-05,   0.00000000e+00,   1.99514634e-05],
        [  1.99381211e-05,   3.99558521e-05,   1.99148908e-05, ...,
           5.98267434e-05,   9.96659401e-06,   1.99097249e-05],
        ..., 
        [  5.11917424e-05,   1.12853151e-04,   5.12245788e-05, ...,
           1.02568538e-05,   0.00000000e+00,   0.00000000e+00],
        [  1.01825599e-05,   3.05091289e-05,   7.61858028e-05, ...,
           1.52862031e-05,   0.00000000e+00,   1.52607008e-05],
        [  1.00746984e-05,   1.05599909e-04,   4.53294342e-05, ...,
           0.00000000e+00,   0.00000000e+00,   0.00000000e+00]]], dtype=float32)

r.timedate
 [datetime.datetime(2012, 1, 1, 2, 0),
 datetime.datetime(2012, 1, 1, 4, 0),
 datetime.datetime(2012, 1, 1, 6, 0),
 ...]

plot:
read in footprint first with read function
then use plot function with read output and date to plot a still.
--> outputs a figure you can then save
if you use the plot function with addfp output:
--> outputs one figure with footprints for all stations files in addfp added together

animate:
first use 'addfp' to add footprints together and/or to resample footprints on the right timeline.
You can either just select two dates, a start and end, to make the timeline or select a different
frequency to the 2 hour frequency of the footprints.

fpanim(fpaddOutput, filedomain = 'If you want to create a separate file')
--> outputs still files in a directory called FP_Animation (or FP_Animation+domain) within the directory that the program is run
this can be turned into a movie using the following in the command line:
MP4: ffmpeg -r 16 -i ?footprints_%05d.png' -f mp4 -vcodec libx264 -pix_fmt yuv420p -intra -qscale 0 footprints.mp4
WMV: ffmpeg -r 16 -i ?footprints_%05d.png' -b 5000k -f asf -vcodec wmv2 -acodec wmav2 footprints.wmv
(change the number 16 to change the speed of the animation (files per second))

The colourbar for the animation stills is fixed between -4 and 0, this works quite well for multiple sites added
together but you might prefer a larger range for single sites.

@author: ew14860
"""

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D
import datetime as dt
import os
import glob
from matplotlib import ticker
from acrg_time.convert import sec2time
import pandas as pd
import bisect
#from numbapro import autojit, vectorize
#from timeit import default_timer as timer
import numexpr
import subprocess
from progressbar import ProgressBar
import json
from acrg_grid import areagrid
import acrg_agage as agage
import xray
from copy import deepcopy
from os.path import split, realpath
from acrg_time import convert
import calendar

fp_directory = '/data/shared/NAME/fp_netcdf/'
flux_directory = '/data/shared/NAME/emissions/'
basis_directory = '/data/shared/NAME/basis_functions/'

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
    
    return  files


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
        print("Can't find files, exiting")
        return None
    else:
        files.sort()
        fp = []
        for f in files:
            fp.append(xray.open_dataset(f))
            
        fp = xray.concat(fp, dim = 'time')

        # If a species is specified, also get flux            
        if species is not None:
            flux_ds = flux(domain, species)
            if flux is not None:
                fp = combine_datasets(fp, flux_ds)
        
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
                   species.lower() + "*.nc"))
    if len(files) == 0:
        print("Can't find flux: " + domain + " " + species)
        return None
        
    flux_ds = []
    for f in files:
        flux_ds.append(xray.open_dataset(f))
    flux_ds = xray.concat(flux_ds, dim = "time")

    return flux_ds


def basis(domain, basis_case = 'voronoi'):
    """
    Read in a basis function file.
    """
    
    files = sorted(glob.glob(basis_directory + domain + "/" +
                    basis_case + "*.nc"))
    if len(files) == 0:
        print("Can't find basis functions: " + domain + " " + basis_case)
        return None
        
    basis_ds = []
    for f in files:
        basis_ds.append(xray.open_dataset(f))
    basis_ds = xray.concat(basis_ds, dim = "time")

    return basis_ds


def combine_datasets(dsa, dsb, method = "ffill"):
    """
    Merge two datasets. Assumes that you want to 
    re-index to the FIRST dataset.
    
    Example:
    
    ds = combine_datasets(dsa, dsb)

    ds will have the index of dsa    
    """
    ds = dsa.merge(dsb.reindex_like(dsa, method))
    return ds


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


def footprints_data_merge(data, domain = "EUROPE", species = "CH4",
                          calc_timeseries = True, average = None):
    """
    Output a dictionary of xray footprint datasets, that correspond to a given
    dictionary of Pandas dataframes, containing mole fraction time series.
    
    Example:
    
    Input dictionary contains time series at Mace Head and Tacolneston:
        
        data = {"MHD": MHD_dataframe, "TAC": TAC_dataframe}

    The dataset must be labeled with "time" index, "mf" and "dmf" columns.
    To combine this with corresponding NAME footprints:
    
        dataset = footprints_data_merge(data)
        
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

        # Get footprints
        site_fp = footprints(site, start = start, end = end,
                             domain = domain, species = species)
        
        if site_fp is not None:
            
            # Merge datasets
            site_ds = combine_datasets(site_ds, site_fp, method = "bfill")
            
            # If units are specified, multiply by scaling factor
            if ".units" in attributes:
                site_ds.fp = site_ds.fp / data[".units"]
            
            # Calculate model time series, if required
            if calc_timeseries:
                site_ds["mf_mod"] = timeseries(site_ds)
    
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
    """    
    
    sites = [key for key in fp_and_data.keys() if key[0] != '.']
    attributes = [key for key in fp_and_data.keys() if key[0] == '.']
    basis_func = basis(domain = domain, basis_case = basis_case)
    fp_data_H = {}
    
    for site in sites:
        site_ds = fp_and_data[site]
        if ".units" in attributes:
            site_ds.fp = site_ds.fp / fp_and_data[".units"]
        site_bf = combine_datasets(site_ds, basis_func)
        
        reference = site_bf.mf_mod
        
        H = np.zeros((len(site_bf.coords['region']),len(reference)))
        
        for i in range(len(site_bf.coords['region'])):
            reg = site_bf.basis.sel(region=i)
            flux_scale = reg + 1
            perturbed = (site_bf.fp*site_bf.flux*flux_scale).sum(["lat", "lon"])
            H[i,:] = perturbed - reference
        
        sensitivity = xray.Dataset({'H': (['region','time'], H)},
                                    coords = {'region': (site_bf.coords['region']),
                                              'time' : (site_ds.coords['time'])})
        site_ds = combine_datasets(site_ds, sensitivity)
        fp_data_H[site] = site_ds
    
    return fp_data_H
    
    
def merge_sensitivity(fp_data_H):
#    outputs y, y_site, y_time, H
    y = []
    y_error = []
    y_site = []
    y_time = []
    H = []
    
    sites = [key for key in fp_data_H.keys() if key[0] != '.']
    for si, site in enumerate(sites):
        y.append(fp_data_H[site].mf.values)
        y_error.append(fp_data_H[site].dmf.values)
        y_site.append([site for i in range(len(fp_data_H[site].coords['time']))])
        y_time.append(fp_data_H[site].coords['time'].values)
        H.append(fp_data_H[site].H.values)
    
    y = np.hstack(y)
    y_error = np.hstack(y_error)
    y_site = np.hstack(y_site)
    y_time = np.hstack(y_time)
    H = np.hstack(H)
    
    return y, y_error, y_site, y_time, H.T
    

def filtering(time, mf, filt):

    def midday(time, mf):
        df=pandas.DataFrame(mf, index=time, columns=['mf'])
        dfpm=df[(df.index.hour>=10) * (df.index.hour<=15)]
        dfr=dfpm.resample("1D", how="median")
        return [t.to_pydatetime() + dt.timedelta(0.5) for t in dfr.index], \
                np.array(dfr['mf'])
    
    def daytime2hr(time,mf):
        df=pandas.DataFrame(mf, index=time, columns=['mf'])
        dfpm=df[(df.index.hour>=10) * (df.index.hour<=16)]
        dfr = dfpm.resample("2H", how="mean")
        dfn = dfr.dropna()
        return [t.to_pydatetime() for t in dfn.index], \
                np.array(dfn['mf'])

    filters={"midday":midday,
             "daytime2hr":daytime2hr}
            
    return filters[filt](time, mf)


def baseline(y, y_time, y_site, x_error = 10000, days_to_average = 5):
    
    keys = np.unique(y_site)

    n = days_to_average
    pos = np.zeros(len(y))
    for site in keys:
        val = np.max(pos)
        wh = np.where(y_site == site)
        ts = pandas.Series(1, y_time[wh])
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
                 lon_range = None, lat_range = None):

        if lon_range is None:
            lon_range = (min(fp_data.lon.values),
                         max(fp_data.lon.values))
        if lat_range is None:
            lat_range = (min(fp_data.lat.values),
                         max(fp_data.lat.values))
        
        m = Basemap(projection='gall',
            llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
            resolution='l')

        lons, lats = np.meshgrid(fp_data.lon.values,
                                 fp_data.lat.values)
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
        return plt.cm.BuPu

def plot_map_zoom(fp_data):
    
    sites = fp_data.keys()
    
    dlon = max(fp_data[sites[0]].lon.values) - \
            min(fp_data[sites[0]].lon.values)
    lon_range = [min(fp_data[sites[0]].lon.values) + 0.6*dlon,
                 max(fp_data[sites[0]].lon.values) - 0.2*dlon]
    dlat = max(fp_data[sites[0]].lat.values) - \
            min(fp_data[sites[0]].lat.values)
    lat_range = [min(fp_data[sites[0]].lat.values) + 0.5*dlat,
                 max(fp_data[sites[0]].lat.values) - 0.25*dlat]

    return lat_range, lon_range

def plot(fp_data, date, out_filename=None, 
         lon_range=None, lat_range=None, cutoff = -3.5,
         map_data = None, zoom = False):
    
    """date as string "d/m/y H:M" or datetime object 
    datetime.datetime(yyyy,mm,dd,hh,mm)
    """
    
    # Looks for nearest time point aviable in footprint   
    date = convert.reftime(date)

    # Get sites
    sites = fp_data.keys()

    # Zoom in. Assumes release point is to the East of centre
    if zoom:
        lat_range, lon_range = plot_map_zoom(fp_data)
        
    # Get map data
    if map_data is None:
        map_data = plot_map_setup(fp_data[sites[0]],
                                  lon_range = lon_range,
                                  lat_range = lat_range)
    
    # Open plot
    fig = plt.figure(figsize=(8,8))
    fig.add_axes([0.1,0.1,0.8,0.8])

    map_data.m.drawcoastlines()
    map_data.m.drawcountries()

    levels = np.arange(cutoff, 0., 0.05)
    
    release_lon = {}
    release_lat = {}

#    data = np.zeros(np.shape(
#                    fp_data[sites[0]].fp.values.squeeze())) - cutoff

    for site in sites:
    
        fp_data_ti = fp_data[site].reindex_like( \
                        xray.Dataset(coords = {"time": [date]}), method = "pad")
        
        data = np.log10(fp_data_ti.fp.values.squeeze())
        
        #Set very small elements to zero
        data[np.where(data <  cutoff)]=np.nan
        
        #Plot map
        cs = map_data.m.contourf(map_data.x, map_data.y, data,
                                 levels, cmap = plot_default_colors(site),
                                 alpha = 0.8,
                                 antialiased = True)

        # Store release location to overplot later
        if "release_lat" in dir(fp_data_ti):
            release_lon[site] = fp_data_ti.release_lon.values
            release_lat[site] = fp_data_ti.release_lat.values

    # over-plot release location
    if len(release_lon) > 0:
        for site in sites:            
            rplons, rplats = np.meshgrid(release_lon[site],
                                         release_lat[site])
            rpx, rpy = map_data.m(rplons, rplats)
            rp = map_data.m.scatter(rpx, rpy, 40, color = 'black')
    
    plt.title(str(date), fontsize=20)

    cb = map_data.m.colorbar(cs, location='bottom', pad="5%")
    
    tick_locator = ticker.MaxNLocator(nbins=7)
    cb.locator = tick_locator
    cb.update_ticks()
 
    cb.set_label('log$_{10}$( (mol/mol) / (mol/m$^2$/s))', 
                 fontsize=15)
    cb.ax.tick_params(labelsize=13) 
    
    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()


def time_unique(fp_data):
    
    sites = fp_data.keys()
    
    time = fp_data[sites[0]].time.to_dataset()
    if len(sites) > 1:
        for site in sites[1:]:
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
            cutoff = -3.5,
            overwrite=True, file_label = 'fp', 
            framerate=10, delete_png=False,
            video_os="mac", ffmpeg_only = False):
    
    sites = fp_data.keys()
    
    if ffmpeg_only is False:

        # Set up map        
        if zoom:
            lat_range, lon_range = plot_map_zoom(fp_data)
            
        map_data = plot_map_setup(fp_data[sites[0]], 
                                  lon_range = lon_range, 
                                  lat_range= lat_range)
        
        # Find unique times
        times = time_unique(fp_data)

        # Start progress bar
        pbar=ProgressBar(maxval=len(times.time.values)).start()

        # Plot each timestep
        for ti, t in enumerate(times.time.values):
            
            fname=os.path.join(output_directory, 
                               file_label + '_' + str(ti).zfill(5) + '.png')
                               
            if len(glob.glob(fname)) == 0 or overwrite == True:            
                plot(fp_data, t, out_filename = fname, 
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
  def __init__(self, domain, ocean=None):

        if ocean is None:

            countryDirectory='/data/shared/NAME/countries/'
            filename=glob.glob(countryDirectory + \
                 "/" + "country_" \
                 + domain + ".nc")
             
        else:
            countryDirectory='/data/shared/NAME/countries/'
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
