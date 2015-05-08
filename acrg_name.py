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
import pandas
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

fp_directory = '/data/shared/NAME/fp_netcdf/'
flux_directory = '/data/shared/NAME/emissions/'
basis_directory = '/data/shared/NAME/basis_functions/'

def filename(site, domain, years, flux=None, basis=None):
        
    if flux is None and basis is None:
        baseDirectory = fp_directory
    else:
        if flux is not None:
            baseDirectory = flux_directory
        if basis is not None:
            baseDirectory = basis_directory

    filename=[]
    fileyear=[]
    
    for year in years:
        files=glob.glob(baseDirectory + \
            domain + "/" + \
            site + "_" + domain + "_" + str(year) + ".nc")
        if len(files) == 0:
            print("Can't find file " + baseDirectory + \
                domain + "/" + \
                site + "_" + domain + "_" + str(year) + ".nc")
        else:
            filename.append(files[0])
            fileyear.append(year)
        
    return  filename, fileyear


def ncread(filenames, varname):

    f = nc.Dataset(filenames[0], 'r')
    
    #Get grid and time first
    lon = f.variables['lon'][:]
    lat = f.variables['lat'][:]
    
    #Get footprints and time
    var = f.variables[varname][:, :, :]
    time = sec2time(f.variables['time'][:], \
        str(f.variables['time'].units)[14:])
    f.close()
    
    if len(filenames) > 1:
        for fpname in filenames[1:]:
            f=nc.Dataset(fpname, 'r')
            varfile = f.variables[varname][:, :, :]
            var = np.append(var, varfile, axis=2)
            time = time + \
                sec2time(f.variables['time'][:], \
                str(f.variables['time'].units)[14:])
            f.close()

    return lat, lon, time, var
    

def read_particle_locations(filenames):
    
    f = nc.Dataset(filenames[0], 'r')
    
    #Get grid and time first
    pl_height = f.variables['height'][:]
    pl_n = [f.variables['particle_locations_n'][:,:,:]]
    pl_e = [f.variables['particle_locations_e'][:,:,:]]
    pl_s = [f.variables['particle_locations_s'][:,:,:]]
    pl_w = [f.variables['particle_locations_w'][:,:,:]]
    f.close()
    
    for f in filenames[1:]:
        pl_n.append(f.variables['particle_locations_n'][:,:,:])
        pl_e.append(f.variables['particle_locations_e'][:,:,:])
        pl_s.append(f.variables['particle_locations_s'][:,:,:])
        pl_w.append(f.variables['particle_locations_w'][:,:,:])
        
    pl = {"N": np.dstack(pl_n),
          "E": np.dstack(pl_e),
          "S": np.dstack(pl_s),
          "W": np.dstack(pl_w)}

    return pl_height, pl

class read:
    def __init__(self, sitecode_or_filename, years=[2012], domain="small"):

        #Chose whether we've input a site code or a file name
        #If it's a three-letter site code, assume it's been processed
        # into an annual footprint file in (mol/mol) / (mol/m2/s)
        # using acrg_name_process
        if '.nc' in sitecode_or_filename:
            if not '/' in sitecode_or_filename:
                filenames = [os.path.join(fp_directory, sitecode_or_filename)]
            else:
                filenames=[sitecode_or_filename]
            fileyears=[0]
        else:
            site=sitecode_or_filename[:]
            if type(years) is not list:
                years = [years]
            filenames, fileyears = filename(site, domain, years)
        
        #Get footprints
        lat, lon, time, fp = ncread(filenames, 'fp')
            
        self.lon = lon
        self.lat = lat
        self.lonmax = np.max(lon)
        self.lonmin = np.min(lon)
        self.latmax = np.max(lat)
        self.latmin = np.min(lat)
        self.fp = np.asarray(fp)
        self.time = time

        if domain is "NWEU":
            #Get particle locations (if available)
            pl_height, pl = read_particle_locations(filenames)

            self.particle_locations = pl
            self.particle_height = pl_height
        

class flux:
    def __init__(self, species_or_filename, years=[2010], domain="small"):

        #Chose whether we've input a species or a file name
        if '/' in species_or_filename and \
            '.nc'in species_or_filename:
            filenames=[species_or_filename]
            fileyears=[0]
        else:
            species=species_or_filename[:]
            if type(years) is not list:
                years = [years]
            filenames, fileyears = \
                filename(species, domain, years, flux=True)
                
        lat, lon, time, flux = ncread(filenames, 'flux')
        
        self.lon = lon
        self.lat = lat
        self.flux = np.asarray(flux)
        self.time = time

class basis_function:    
    def __init__(self, case_or_filename, years=[2012], domain='small'):
    
        if '/' in case_or_filename and \
            '.nc'in case_or_filename:
            filenames=[case_or_filename]
        else:
            case=case_or_filename[:]
            if type(years) is not list:
                years = [years]
            filenames, fileyears = \
                filename(case, domain, years, basis=True)
    
        lat, lon, time, basis = ncread(filenames, 'basis')
        
        self.lat = lat
        self.lon = lon
        self.time = time
        self.basis = basis


def timeseries(site, species, years, flux_years=None, domain='small',
               filt=None):
  
    if flux_years is None:
        flux_years=years
    
    fp_data=read(site, years, domain=domain)
    flux_data=flux(species, flux_years, domain=domain)

    return footprint_x_flux(fp_data, flux_data, filt=filt)


def footprint_x_flux(fp_data, flux_data, basis=None, basis_scale=1.,
                     filt=None):

    #Calculate flux timestep
    if len(flux_data.time) == 1:
        flux_time_delta=dt.timedelta(0) 
    else:
        flux_time_delta=(flux_data.time[1] - flux_data.time[0])/2

    #Calculate scaling to apply to flux field
    #CONSTANT SCALING IN TIME AT THE MOMENT!
    flux_scale=np.ones((flux_data.flux[:,:,0]).shape)
    if basis is not None:
        for si, scale in enumerate(basis_scale):
            wh=np.where(basis == si + 1) # Basis numbers should start at 1
            flux_scale[wh] = flux_scale[wh]*scale
    
    #Calculate time series
    mf=np.zeros(len(fp_data.time))
    for ti, t in enumerate(fp_data.time):
        #find nearest point in emissions dataset
        flux_ti=min([len(flux_data.time)-1,
                     bisect.bisect(flux_data.time, t - flux_time_delta)])
        fp=fp_data.fp[:, :, ti]
        flux=flux_data.flux[:, :, flux_ti]
#        mf[ti]=0.
        mf[ti]=numexpr.evaluate("sum(fp*flux*flux_scale)")
    
    time=fp_data.time
    
    if filt is not None:
        time, mf = filtering(time, mf, filt)

    return time, mf


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


def sensitivity_single_site(site, species,
                years=[2012], flux_years=None,
                domain="small", basis_case='voronoi', filt=None):
    
    if flux_years is None:
        flux_years=years
    
    fp_data=read(site, years, domain=domain)
    flux_data=flux(species, flux_years, domain=domain)    
    basis_data = basis_function(basis_case, years=years, domain=domain)

    basis_scale=np.ones(np.max(basis_data.basis))
    
    time, reference = footprint_x_flux(fp_data, flux_data, 
                    basis=basis_data.basis[:,:,0], basis_scale=basis_scale, 
                    filt=filt)

    sensitivity=np.zeros((len(reference), len(basis_scale)))
    
    
    for xi, scale in enumerate(basis_scale):
        basis_scale_perturbed=basis_scale.copy()
        basis_scale_perturbed[xi] += 1.
        time, perturbed = footprint_x_flux(fp_data, flux_data, 
                        basis=basis_data.basis[:,:,0], 
                        basis_scale=basis_scale_perturbed, 
                        filt=filt)
        sensitivity[:, xi] = perturbed - reference

    return time, sensitivity


def sensitivity(obs, species, years=[2012], flux_years=None,
                domain="small", basis_case='voronoi', filt=None, alt_fp_filename = None):
#Using alt_fp_filename only works for one site at the moment
    H=[]
    y_time=[]
    y_site=[]
    y=[]

    for site in sorted(obs.iterkeys()):
        if alt_fp_filename is not None:
            ts, Hs = sensitivity_single_site(alt_fp_filename, species, years, 
                                flux_years=flux_years, domain=domain,
                                basis_case = basis_case, filt=filt)
        else:
            ts, Hs = sensitivity_single_site(site, species, years, 
                                flux_years=flux_years, domain=domain, 
                                basis_case = basis_case, filt=filt)
        df_site = pandas.DataFrame(Hs, index=ts)
        obsdf = obs[site].dropna()
        Hdf = df_site.reindex(obsdf.index)
        Hdf2 = Hdf.dropna()
        obsdf2 = obsdf.reindex(Hdf2.index)
        y_time.append(Hdf2.index.to_pydatetime())
        y_site.append([site for i in range(len(Hdf2.index))])
        y.append(obsdf2.values)
        H.append(Hdf2.values)
        
    H=np.vstack(H)
    y_time=np.hstack(y_time)
    y_site=np.hstack(y_site)
    y=np.vstack(y)

    return y_time, y_site, y, H


def baseline(y, y_time, y_site, days_to_average = 5):
    
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
        xerror[col] = 10000
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

class addfp:
    def __init__(self, filenames, start, end, frequency = 2):
        
        timeslice = pandas.date_range(start, end, freq = '%iH' %frequency)
        filenames.sort()    
    
        lat, lon, time, fp = ncread(filenames, 'fp')
        time = np.array(time)
    
        fpadd = np.zeros((len(lat),len(lon),len(timeslice)))
       
        for i in range(len(timeslice)):
            index0 = np.where(time==timeslice[i])
            index = list(index0[0])
            if len(index)==0:
                fpadd[:,:,i]=np.empty(np.shape(fpadd[:,:,0]))
                fpadd[:,:,i].fill('NaN')
            elif len(index)==1:
                fpadd[:,:,i]=fp[:,:,int(index0[0])]
            elif len(index)>1:
                fpadd[:,:,i]=np.sum(fp[:,:,index], axis=2)
            
        self.lat = lat
        self.lon = lon
        self.lonmax = np.max(lon)
        self.lonmin = np.min(lon)
        self.latmax = np.max(lat)
        self.latmin = np.min(lat)
        self.time = list(timeslice)
        self.fp = fpadd


class plot_map_setup:
    def __init__(self, fp_data, 
                 lon_range = None, lat_range = None):

        if lon_range is None:
            lon_range = (fp_data.lonmin, fp_data.lonmax)
        if lat_range is None:
            lat_range = (fp_data.latmin, fp_data.latmax)
        
        m = Basemap(projection='gall',
            llcrnrlat=lat_range[0], urcrnrlat=lat_range[1],
            llcrnrlon=lon_range[0], urcrnrlon=lon_range[1],
            resolution='l')

        lons, lats = np.meshgrid(fp_data.lon, fp_data.lat)
        x, y = m(lons, lats)
        
        self.x = x
        self.y = y
        self.m = m

def plot(fp_data, date, out_filename=None, 
         lon_range=None, lat_range=None, cutoff = -3.5,
         map_data = None):

    """date as string "d/m/y H:M" or datetime object 
    datetime.datetime(yyyy,mm,dd,hh,mm)
    """

    if map_data is None:
        map_data = plot_map_setup(fp_data, 
                                  lon_range = lon_range,
                                  lat_range = lat_range)
    
    #looks for nearest time point aviable in footprint   
    if isinstance(date, str):
        date=dt.datetime.strptime(date, '%d/%m/%Y %H:%M')

    time_index = bisect.bisect_left(fp_data.time, date)

    data = np.log10(fp_data.fp[:,:,time_index])

    #Set very small elements to zero
    data[np.where(data <  cutoff)]=str("nan")

    fig = plt.figure(figsize=(8,8))
    fig.add_axes([0.1,0.1,0.8,0.8])

    map_data.m.drawcoastlines()
    map_data.m.drawstates()
    map_data.m.drawcountries()
    
    levels = np.arange(cutoff, 0., 0.05)
    
#    cm = plt.cm.Purples
    cs = map_data.m.contourf(map_data.x,map_data.y,data,
                             levels)
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

    fig = plt.figure(figsize=(8,7))

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

    if out_filename is not None:
        plt.savefig(out_filename)
        plt.close()
    else:
        plt.show()


def animate(allfpdata, output_directory, 
            lon_range = None, lat_range=None, 
            cutoff = -3.5,
            overwrite=True, file_label = 'fp', 
            framerate=10, delete_png=False,
            video_os="mac", ffmpeg_only = False):
    
    if ffmpeg_only is False:
        map_data = plot_map_setup(allfpdata, 
                                  lon_range = lon_range, 
                                  lat_range= lat_range)
        
        pbar=ProgressBar(maxval=len(allfpdata.time)).start()
        for ti, t in enumerate(allfpdata.time):
            
            fname=os.path.join(output_directory, 
                               file_label + '_' + str(ti).zfill(5) + '.png')
                               
            if len(glob.glob(fname)) == 0 or overwrite == True:            
                plot(allfpdata, t, out_filename = fname, 
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
            
def fp_resample(fp,lat,lon,time, av_period=None, dimension='time', av_how='mean',
                startM=None, endM=None):
                    
    da = xray.DataArray(fp, [('lat', lat),('lon', lon), 
                                 ('time', time)], name = 'fp')
                                 
    if av_period is not None:
        da = da.resample(av_period, dimension, how=av_how, base=1)
    
    if startM is not None:
        if endM is None:
                print 'Need to specify an end month as well'
                return
        else:
            if endM == startM:
                print 'End month needs to be different from start month'
                return
                
            year=time[1].timetuple().tm_year
            da = da.sel(time=slice(dt.datetime(year,startM,01), dt.datetime(year,endM,01)))
    
    fp_out=da
    time_out = da.time    
    
    return fp_out, time_out
   
