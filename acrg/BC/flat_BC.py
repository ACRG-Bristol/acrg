from openghg.retrieve import get_obs_surface
import xarray as xr
import numpy as np

def create_flat_bc_prior(species,start_date,end_date,bc_dir=None,save_name=None):
    '''
    Takes the mean MHD obs over the time period and uses 
    this value to create a uniform bc curtain prior file, based on
    the format of files used by the ACRG and OpenGHG.
    '''

    if species == 'ch4':
        units_scale = 1e-9
    elif 'hfc' in species:
        units_scale = 1e-12

    if not bc_dir:
        bc_dir = '/group/chemistry/acrg/LPDM/bc/EUROPE'

    bc_in = xr.open_dataset(f'{bc_dir}/EUROPE/ch4_EUROPE_201901.nc')
    bc_out = bc_in.copy()
    
    mhd_o = get_obs_surface(site='MHD',
                                        species=species.lower(),
                                        inlet='10m',
                                        start_date=start_date,
                                        end_date=end_date,
                                        average='4H',
                                        instrument=None,
                                        store='shared_obs')
        
    note = None
    
    if mhd_o == None:
        
        mhd_o = get_obs_surface(site='JFJ',
                                        species=species.lower(),
                                        inlet=None,
                                        start_date=start_date,
                                        end_date=end_date,
                                        average='4H',
                                        instrument=None,
                                        store='shared_obs')
        note = 'JFJ used for background as no MHD obs available'
        print('JFJ!')
        
    mhd_obs = mhd_o.data['mf'].values
    mhd_times = mhd_o.data.time.values

    #mhd_obs,mhd_times = extract_acrg_obs('MHD',species,start_date,end_date)

    mean_mhd = np.round(np.nanmean(mhd_obs),3)
    mhd_25perc = np.percentile(mhd_obs,25)

    bc_n = np.ones(bc_out['vmr_n'].values.shape) * mhd_25perc * units_scale #mean_mhd
    bc_e = np.ones(bc_out['vmr_e'].values.shape) * mhd_25perc * units_scale #mean_mhd
    bc_s = np.ones(bc_out['vmr_s'].values.shape) * mhd_25perc * units_scale #mean_mhd
    bc_w = np.ones(bc_out['vmr_w'].values.shape) * mhd_25perc * units_scale #mean_mhd

    bc_out['vmr_n'].values = bc_n
    bc_out['vmr_e'].values = bc_e
    bc_out['vmr_s'].values = bc_s
    bc_out['vmr_w'].values = bc_w
    bc_out['time'] = np.datetime64(start_date)

    bc_out.attrs['title'] = f'Uniform boundary conditions for {species} using the monthly 25th percentile of MHD obs: {mean_mhd}.'
    del(bc_out.attrs['CAMS_resolution'])
    bc_out.attrs['author'] = 'aramsden'
    bc_out.attrs['date_created'] = str(np.datetime64('now'))
    
    if note is not None:
        bc_out.attrs['note'] = note

    print(f'\nBC file created for date {start_date} with value {mean_mhd}\n')

    if save_name:
        save_path = f'{bc_dir}/EUROPE/{save_name}.nc'
        print(f'Saving to {save_path}\n')

        bc_out.to_netcdf(save_path)

    #fig,ax = plt.subplots(1,1,figsize=(8,8))
    #ax.plot(np.arange(mhd_obs.shape[0]),np.ones(mhd_obs.shape[0])*mhd_25perc)
    #ax.plot(mhd_obs)

    return bc_out#,fig