#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 13:29:24 2019

@author: ew14860
"""

import iris
import xarray as xr
import gzip
import numpy as np
import pandas as pd
import os
import math
import json
import collections as c

acrg_path = os.getenv("ACRG_PATH")

if acrg_path is None:
    acrg_path = os.getenv("HOME")
    print("Default ACRG directory is assumed to be home directory. Set path in .bashrc as \
            export ACRG_PATH=/path/to/acrg/repository/ and restart python terminal")

# Function to get wind speed and wind direction from x_wind and y_wind

def xy_wind_to_speed_direction(x_wind, y_wind):
    
    """
    Input x_wind and y_wind from met files as arrays
    Returns arrays of wind speed and wind direction
    """
    
    # Wind speed is hypotenuse of x and y wind speeds
    wind_speed = np.sqrt(x_wind**2+y_wind**2)
    
    # Wind direction is angle from North of x_wind, y_wind vector
    # Need to rotate x and y wind 90 degrees clockwise to use math.atan2 function
    # Then find angle from 0-360 by removing answer from 360 (and modulus to one rotation of 360 degrees)
    
    x_rot = y_wind
    y_rot = x_wind*-1
    
    wind_direction = np.rad2deg(2*math.pi - np.arctan2(y_rot,x_rot))%360
    
    return wind_speed, wind_direction
    

# Function to work out correct MK
    
def get_met_Mark(date):

    """
    Input date as string of format "YYYY-MM-DD HH:MM:SS"
    Returns Mark number for met
    """
    
    if date >= "2009-10-06 00:00" and date < "2010-03-09 00:00":
        Mk = 5
    elif date >= "2010-03-09 00:00" and date < "2013-04-30 00:00":
        Mk = 6
    elif date >= "2013-04-30 00:00" and date < "2014-07-15 00:00":
        Mk = 7
    elif date >= "2014-07-15 00:00" and date < "2015-08-25 00:00":
        Mk = 8
    elif date >= "2015-08-25 00:00" and date < "2017-07-11 00:00":
        Mk = 9
    elif date >= "2017-07-11 00:00" and date <= "2018-05-31 21:00":
        Mk = 10

    return Mk


# Function to work out correct PT

def get_met_Part(lat, lon):
    
    """
    Input lat (-90 - 90) and lon (-180 - 180) - defaults from acrg_site_info.json
    Returns Part number of met

    PT 1: lat 79.921875 - 89.921875, lon 0.1171875 - 359.8828
    PT 2: lat 24.921875 - 80.078125, lon 314.8828 - 405.1172
    PT 3: lat 24.921875 - 80.078125, lon 44.882812 - 135.11719
    PT 4: lat 24.921875 - 80.078125, lon 134.88281 - 225.11719
    PT 5: lat 24.921875 - 80.078125, lon 224.88281 - 315.1172
    PT 6: lat -25.078125 - 25.078125, lon 314.8828 - 405.1172
    PT 7: lat -25.078125 - 25.078125, lon 44.882812 - 135.11719
    PT 8: lat -25.078125 - 25.078125, lon 134.88281 - 225.11719
    PT 9: lat -25.078125 - 25.078125, lon 224.88281 - 315.1172
    PT 10: lat -80.078125 - -24.921875, lon 314.8828 - 405.1172
    PT 11: lat -80.078125 - -24.921875, lon 44.882812 - 135.11719
    PT 12: lat -80.078125 - -24.921875, lon 134.88281 - 225.11719
    PT 13: lat -80.078125 - -24.921875, lon 224.88281 - 315.1172
    PT 14: lat -89.921875 - 79.921875, lon 0.1171875 - 359.8828
    """
    
    #Convert lat lon to convention for met ( lat: -90 to 90, lon: 0 - 360)
    
    lon = lon + 180.

    if lat > 79.921875:
        PT = 1
    elif lat > 24.921875 and lat < 80.078125:
        if lon > 314.8828:
            PT = 2
        elif lon < 45.1172:
            PT = 2
        elif lon > 44.882812 and lon < 135.11719:
            PT = 3
        elif lon > 134.88281 and lon < 225.11719:
            PT = 4
        elif lon > 224.88281 and lon < 315.1172:
            PT = 5
    elif lat > -25.078125 and lat < 25.078125:
        if lon > 314.8828:
            PT = 6
        elif lon < 45.1172:
            PT = 6
        elif lon > 44.882812 and lon < 135.11719:
            PT = 7
        elif lon > 134.88281 and lon < 225.11719:
            PT = 8
        elif lon > 224.88281 and lon < 315.1172:
            PT = 9
    elif lat > -80.078125 and lat < -24.921875:
        if lon > 314.8828:
            PT = 10
        elif lon < 45.1172:
            PT = 10
        elif lon > 44.882812 and lon < 135.11719:
            PT = 11
        elif lon > 134.88281 and lon < 225.11719:
            PT = 12
        elif lon > 224.88281 and lon < 315.1172:
            PT = 13
    elif lat < -79.921875:
        PT = 14

    return PT


def get_vertical_profile(site, start_date, end_date, output_vars, temp_dir, output_dir, met_dir = "/mnt/storage/private/acrg/met_archive/NAME/Met/Global/"):
    
    """
    Function extracts a vertical meteorological profile at the location requested (the nearest grid cell in the met files) for the required variables
    
    site: site code for the location of the vertical profile to be extracted
    start_date: start of vertical profile timeseries
    end_date: end of vertical profile timeseries
    output_vars: options are air_pressure, air_temperature, specific_humidity, wind_speed, wind_direction
                            (wind_speed and wind_direction calculated automatically from x_wind and y_wind)
    temp_dir: pathname to directory where met files are stored as they are unzipped (they are deleted once used)
    output_dir: where to save the netcdf file of the vertical profile
    met_dir: where the met is stored (default is storage location on BC4)
    """
    
    # Get lat lon from site_name

    with open(acrg_path + "/acrg_site_info.json") as f:
            site_info=json.load(f)
    lat_coord = site_info[site]['latitude']
    lon_coord = site_info[site]['longitude']

    # If wind_speed and wind_direction required, add x_wind and y_wind to output_vars

    if "wind_speed" in output_vars:
        output_vars = output_vars + ["x_wind","y_wind"]
    elif "wind_direction" in output_vars:
        output_vars = output_vars + ["x_wind","y_wind"]
    
    # Find correct filename - get right MK and right PT based on dates and coordinates
    # Then need to find correct files from date range

    date_range = pd.date_range(start=start_date, end = end_date, freq = '3H')

    for di, d in enumerate(date_range):
    
        file_date = "%04d%02d%02d%02d%02d" %(d.year, d.month, d.day, d.hour, d.minute)
    
        Mark = get_met_Mark(str(d))
        Part = get_met_Part(lat_coord, lon_coord)

        fname = met_dir + "UMG_Mk" + str(Mark) + "PT/" "MO" + file_date + ".UMG_Mk" + str(Mark) + "_I_L59PT" + str(Part) + ".pp.gz"
        temp_dest = temp_dir + fname.split('/')[-1][:-3]

        # Unzip file

        with gzip.open(fname, 'rb') as infile:
            with open(temp_dest, 'wb') as outfile:
                for line in infile:
                    outfile.write(line)

        # Load iris cube

        met = iris.load(temp_dest)

        # Get the cube corresponding to the desired variable

        variables = [var.name() for var in met]

        i = -1
        for vi, v in enumerate(variables):
            if v in output_vars:
                print v
                cube = met[vi]
                if len(cube.shape) == 3:
                    i+=1
            
                    # Find lat and lon closest to desired point
            
                    lon = cube.coord('longitude').points
                    lat = cube.coord('latitude').points
                
                    #Make sure lon_coord is in the right format (originally -180 - 180, change to 45.1172 - 405.1172)
                    if lon_coord < 45.1172:
                        lon_coord = lon_coord +360
            
                    wh_lon = (np.abs(lon - lon_coord)).argmin()
                    wh_lat = (np.abs(lat - lat_coord)).argmin()
            
                    # Extract vertical profile
            
                    constraint = iris.Constraint(latitude=lat[wh_lat], longitude = lon[wh_lon])
                    vertical_profile = cube.extract(constraint).data
            
                    #Get datetime from time
                    time = cube.coord('time')
                    dates = time.units.num2date(time.points)
            
                    da = xr.DataArray(vertical_profile[None,:], coords=[dates,cube.coord('model_level_number').points], dims=['time','level'], name = v,
                                      attrs = {'units': str(cube.units), 'lat': str(lon[wh_lon]), 'lon': str(lat[wh_lat])})
            
                    if i == 0:
                        ds = da.to_dataset()
                    else:
                        ds = xr.merge([ds,da])
        
                else:
                    pass
            else:
                pass

        if di == 0:
            DS = ds
        else:
            DS = xr.concat((DS,ds), dim = 'time')
        
        # Remove temp met file
        os.remove(temp_dest)

    # If x and y wind present in DS: Get wind speed and direction

    if "x_wind" in DS.keys():
        if "y_wind" in DS.keys():
            wind_speed, wind_direction = xy_wind_to_speed_direction(DS.x_wind.values, DS.y_wind.values)
        
            if "wind_speed" in output_vars:
                ws = xr.DataArray(wind_speed, coords=[DS.time,DS.level], dims=['time','level'], name = "wind_speed", attrs = {'units': DS.x_wind.units})
                DS = xr.merge([DS,ws])
            if "wind_direction" in output_vars:
                wd = xr.DataArray(wind_direction, coords=[DS.time,DS.level], dims=['time','level'], name = "wind_direction", attrs = {'units': "Degrees from North"})
                DS = xr.merge([DS,wd])

    # Convert levels to heights in m asl
    #(taken from /mnt/storage/private/acrg/NAME/Model/NAMEIII_v7_2_particlelocation/Resources/Defns/MetDefnUMG_Mk9_L59PTpp.txt)
            
    heights = np.array([10.000,36.667,76.667,130.000,196.667,276.667,370.000,476.667,596.667,730.000,876.667,1036.667,
                        1210.000,1396.667,1596.667,1810.000,2036.667,2276.667,2530.000,2796.667,3076.667,3370.001,
                        3676.667,3996.667,4330.001,4676.668,5036.667,5410.001,5796.668,6196.667,6610.001,7036.672,
                        7476.690,7930.085,8396.920,8877.306,9371.434,9879.598,10402.239,10939.989,11493.719,12064.604,
                        12654.205,13264.532,13898.152,14558.290,15248.947,15975.034,16742.518,17558.574,18431.783,
                        19372.307,20392.113,21505.199,22727.850,24078.898,25580.037,27256.115,29135.490])

    DS = DS.update({'level': heights})
    DS = DS.rename({'level': "height"})
    DS.height.attrs = c.OrderedDict([("units", "masl")])
    
    # Give DS attributes
    glob_attrs = c.OrderedDict([("title" , "Vertical meteorology profiles (masl)"),
                                ("site" , site),
                                ("site_height_asl", site_info[site]['height_station_masl']),
                                ("comments" , "Vertical profiles taken from nearest grid cell to site in raw meterology files from Met Office (.pp files)")])

    DS.attrs = c.OrderedDict(glob_attrs)
        
    # Compress and save as netcdf
    
    for key in DS.keys():
            DS[key].encoding['zlib'] = True                    
    DS.to_netcdf(path=output_dir+"vertical_profile_" + site.upper() + ".nc", mode='w')



