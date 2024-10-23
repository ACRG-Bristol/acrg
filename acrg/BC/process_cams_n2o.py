# *****************************************************************************
# Created: 14 Feb. 2023
# Author: Eric Saboya, School of Geographical Sciences, University of Bristol
# Contact: eric.saboya@bristol.ac.uk
# *****************************************************************************
# Updated: 02 Oct. 2024
# Author: Helene De Longueville, School of Chemistry, University of Bristol
# Contact: helene.delongueville@bristol.ac.uk
# *****************************************************************************
# About:
#   Script for processing CAMS Inversion N2O files to match
#   format of CH4 CAMS files for calculating BC conditions
# *****************************************************************************

import os
import numpy as np
import xarray as xr
import datetime as dt
import argparse

def pressure_to_altitude(p_level, p_surf):
    """
    Function for converting atmospheric pressure levels
    to altitudes from using scale height and a simple
    atmospheric pressure model. 
    Uses the average scale height for lowest 0-70 km 
    part of the atmosphere

    More info: https://en.wikipedia.org/wiki/Scale_height
    
    Args:
    -----------------------------------
      p_level (float/ array): 
        Pressure value at each model level
      p_surf (float/ array):
        Surface pressure values 
    
    Returns:
      Correpsonding altitude in metres
    -----------------------------------
    """
    scale_height = 7.64e3 #in metres
    z = -scale_height * np.log(p_level/p_surf)
    return z

def calculate_altitude_levels(p_surf, a, b, nt, nlat, nlon, nlvl):
    """
    Function to calculate the corresponding altitudes
    for the CAMS model levels
    -----------------------------------
    Args:
      p_surf (array)
        Surface pressure values
      a (array)
        Coefficient 'a' as used for calculating 
        the vertical pressure values
      b (array)
        Coefficient 'b' as used for calculting 
        the vertical pressure values 
      nt (int)
        Size of time dimension
      nlat (int)
        Size of latitude dimension
      nlon (int)
        Size of longitude dimension 
      nlvl (int)
        Size of vertical dimension (i.e. no. model levels)
    
    Returns:
      Array of altitudes and pressures for different model levels
    -----------------------------------
    """
    # Calculate pressure half-levels
    p_levels = np.zeros((nt, nlvl, nlat, nlon))
    for i in range(0, nlvl):
        p_levels[:,i,:,:] = a[i] + b[i]*p_surf[:,:,:]
    
    # Calculate altitudes
    z_levels = np.zeros((nt, nlvl, nlat, nlon))
    for i in range(0, nlvl):
        z_levels[:,i,:,:] = pressure_to_altitude(p_levels[:,i,:,:], p_surf)

    return z_levels, p_levels

def daily_means(df):
    """
    Function to calculate daily means from CAMS N2O files
    Calculates daily mean variables for:
      fluxes, 
      surface pressure, 
      pressure levels, 
      altitude,
    -----------------------------------
    Args:
      df (loaded CAMS N2O file)

    Returns
      xarray dataset object with daily means and altitudes
    -----------------------------------
    """
    # Extract relevant data
    time=df['time'].values
    lon=df['longitude'].values
    lat=df['latitude'].values
    n2o=df['N2O'].values
    psurf=df['Psurf'].values
    a=df['ap'].values
    b=df['bp'].values 
    level=df['level'].values
    hlevel=df['hlevel'].values

    nt=len(time)
    nlat=len(lat)
    nlon=len(lon)
    nlvl=len(level)
    
    # Calculate altitudes and pressure levels
    altitudes, plvl=calculate_altitude_levels(psurf, a, b, nt, nlat, nlon, nlvl)
    
    year0=str(df['time'].values[0])[0:4]
    month0=str(df['time'].values[0])[5:7]
    day0=str(df['time'].values[0])[8:10]
    dayn=int(str(df['time'].values[-2])[8:10])
    basetime=dt.datetime(int(year0), int(month0), int(day0), 12, 0, 0)

    dailyt=[]
    for i in range(0, dayn):
        dailyt.append(basetime + dt.timedelta(days=i))

    n2o_dm=np.zeros(shape=(dayn, nlvl, nlat, nlon))
    psurf_dm=np.zeros(shape=(dayn, nlat, nlon))
    plvls_dm=np.zeros(shape=(dayn, nlvl, nlat, nlon))
    alt_dm=np.zeros(shape=(dayn, nlvl, nlat, nlon))

    # Calculate daily means
    for i in range(len(dailyt)):
        tdelta_sec=(df['time'].values - np.datetime64(dailyt[i]))/1e9
        #dailyinds=np.where(np.abs(tdelta_sec.astype(float))<86400)[0]
        dailyinds= np.intersect1d(np.where(tdelta_sec.astype(float)<43200), np.where(tdelta_sec.astype(float)>=-43200))

        n2o_dm[i,:,:,:]=np.mean(n2o[dailyinds,:,:,:], axis=0)
        psurf_dm[i,:,:]=np.mean(psurf[dailyinds,:,:], axis=0)     
        plvls_dm[i,:,:,:]=np.mean(plvl[dailyinds,:,:,:], axis=0)
        alt_dm[i,:,:,:]=np.mean(altitudes[dailyinds,:,:,:], axis=0)
    
    # Create new xarray object with data
    data_vars = {
        'N2O' : (['time','level','latitude','longitude'], n2o_dm,
        {'units': '1e-9 mol mol-1', 
         'long_name': 'N2O dry mole fraction'}),

        'ps' : (['time','latitude','longitude'], psurf_dm,
        {'units': 'Pa',
         'long_name': 'Surface pressure'}),
        
        'altitude' : (['time','level','latitude','longitude'], alt_dm,
        {'units': 'm',
         'long_name':'Altitude'}),
        
        'pressure' : (['time','level','latitude','longitude'], plvls_dm,
        {'units': 'Pa',
        'long_name': 'Atmospheric pressure at corresponding model levels'}),

        }

    coords = {'time': (['time'], np.array(dailyt)),
              'level': (['level'], level),
              'latitude': (['latitude'], lat),
              'longitude': (['longitude'], lon),
              'hlevel': (['hlevel'], hlevel)}
    
    attrs = {
        **df.attrs,
        'altitude_notes': 'Calculated using pressure values and average scale height',
        'author': os.getenv('USER'),  
        'date_created': dt.datetime.today().strftime('%Y-%m-%d')
    }

    dm = xr.Dataset(data_vars=data_vars, 
                    coords=coords, 
                    attrs=attrs)
    
    return dm 
    
# def main():
#   outdir = '/user/work/qq24644/ECMWF_CAMS/CAMS_inversion_processed/'
#   cams_inv_dir = '/group/chem/acrg/ECMWF_CAMS/N2O/'
#   cams_files = [f for f in os.listdir(cams_inv_dir) if not f.endswith('.zip')]
#   # outdir = '/user/work/wz22079/ECMWF_CAMS/CAMS_inversion_processed/'
#   # cams_inv_dir = '/user/work/wz22079/ECMWF_CAMS/CAMS_inversion/N2O/'
#   # cams_files = os.listdir(cams_inv_dir)
  
#   for fname in cams_files:
#     print(fname)
#     df = xr.open_dataset(os.path.join(cams_inv_dir, fname))
#     # fname_save = fname[0:30]+'dm_'+fname[35::]
#     fname_save = fname[0:30]+'_alt'+fname[35::]
#     dm = daily_means(df)
#     dm.to_netcdf(path = os.path.join(outdir, fname_save), mode='w')
   

# if __name__ == '__main__':
#     main() 

def main(cams_inv_dir, outdir):
  cams_files = [f for f in os.listdir(cams_inv_dir) if f.endswith('.nc')]

  for fname in cams_files:
    print(fname)
    df = xr.open_dataset(os.path.join(cams_inv_dir, fname))
    fname_save = fname[0:30]+'alt_'+fname[35::]
    dm = daily_means(df)
    dm.to_netcdf(path = os.path.join(outdir, fname_save), mode='w')
   

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process CAMS inversion data.")
    parser.add_argument('--cams_inv_dir', default='/group/chem/acrg/ECMWF_CAMS/N2O/v22r1/', help='Directory containing CAMS inversion files')
    parser.add_argument('--outdir', required=True, help='Output directory for results')

    args = parser.parse_args()

    # Call the main function with parsed arguments
    main(args.cams_inv_dir, args.outdir)
