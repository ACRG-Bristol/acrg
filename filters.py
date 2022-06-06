from acrg.convert import convert_lons_0360
import numba
import numpy as np
import xarray as xr

def filtering(fp_data_H,sites,heights,species,filtering_types,
              network=None,secondary_heights=None,
              start_date=None,end_date=None,average=None):
    """
    Filters observations in fp_data_H dataset, to account for local influence
    and PBLH thresholds.
    Based on Mark Lunt's code and advice from Alistair Manning.
    Input variables after 'filtering_types' only required for height filtering.
    
    filtering_types = ['localness','pblh','height']
    
    Written by Alice, adapted for Arctic domain by Rebecca.
    """
    
    @numba.jit()
    def assign_localness(fp_data_H_site_fp,wh_rlat,wh_rlon):
        local_sum=np.zeros(fp_data_H_site_fp.shape[2])
        for ti in range(fp_data_H_site_fp.shape[2]):
            fp_data_H_site_local = fp_data_H_site_fp[wh_rlat[0]-2:wh_rlat[0]+3,
                                                     wh_rlon[0]-2:wh_rlon[0]+3,ti]
            local_sum[ti] = np.sum(fp_data_H_site_local)/np.sum(fp_data_H_site_fp[:,:,ti])
        return local_sum
    
    def define_localness(fp_data_H,site):
        """
        Define the localness of each time point for each site.
        Sum up the 25 grid boxes surrounding the site.
        """
        release_lon=fp_data_H[site].release_lon[0].values
        if any(fp_data_H[site].lon.values < 0) & (release_lon < min(fp_data_H[site].lon.values)):
            release_lon = convert_lons_0360(release_lon)
        release_lat=fp_data_H[site].release_lat[0].values
        dlon=fp_data_H[site].lon[1].values-fp_data_H[site].lon[0].values
        dlat=fp_data_H[site].lat[1].values-fp_data_H[site].lat[0].values
        wh_rlon = np.where(abs(fp_data_H[site].lon.values-release_lon) < dlon/2.)[0]
        wh_rlat = np.where(abs(fp_data_H[site].lat.values-release_lat) < dlat/2.)[0]
        fp_data_H_site_fp = fp_data_H[site].fp.values
        local_sum = assign_localness(fp_data_H_site_fp,wh_rlat,wh_rlon)
        local_ds = xr.Dataset({'local_ratio': (['time'], local_sum)},
                                coords = {'time' : (fp_data_H[site].coords['time'])})
        fp_data_H[site] = fp_data_H[site].merge(local_ds)
        return fp_data_H
    
    def localness_filter(dataset,site,keep_missing=False):
        """
        Subset for times when local influence is below threshold.       
        Local influence expressed as a fraction of the sum of entire footprint domain.
        """
        pc = 0.04  #localness 'ratio' limit - originally 0.04
        lr = dataset.local_ratio
        ti = [i for i, local_ratio in enumerate(lr) if local_ratio <= pc]
        if keep_missing is True: 
            mf_data_array = dataset.mf
            dataset_temp = dataset.drop('mf')
            dataarray_temp = mf_data_array[dict(time = ti)]   
            mf_ds = xr.Dataset({'mf': (['time'], dataarray_temp)}, 
                                  coords = {'time' : (dataarray_temp.coords['time'])})
            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]
        
    def pblh_filter(dataset,site,keep_missing=False):
        """
        Subset for times when boundary layer height > threshold.
        Threshold needs to be set in dataset as pblh_threshold.
        Only works for sites as height is variable for aicraft.
        """
        threshold = dataset.pblh_threshold.values
        ti = [i for i, pblh in enumerate(dataset.PBLH) if pblh > threshold]
        if keep_missing:
            mf_data_array = dataset.mf
            dataset_temp = dataset.drop('mf')
            dataarray_temp = mf_data_array[dict(time = ti)]
            mf_ds = xr.Dataset({'mf': (['time'], dataarray_temp)},
                                   coords = {'time' : (dataarray_temp.coords['time'])})
            dataset_out = combine_datasets(dataset_temp, mf_ds, method=None)
            return dataset_out
        else:
            return dataset[dict(time = ti)]   
        
    def height_filter(fp_data_H,site,network,heights,secondary_heights,
                     start_date,end_date,average,species):
        """
        Filters obs by comparing mf observed from multiple inlets at the same height.
        If mole fractions observed at all heights are within 20 ppm (ppb for C2H6) 
        then the observation is kept.
        No filtering applied to sites with only one height.
        """
        data_path = '/user/home/cv18710/work/'
        fp_data_H_out = fp_data_H.copy()
        if secondary_heights is not None:
            obs1 = fp_data_H.mf.values
            times1 = fp_data_H.time.values
            with io.capture_output() as captured:
                data = read.get_obs(sites=[site],species=species,start_date=start_date,end_date=end_date,
                                average=average,network=[network],inlet=[secondary_heights],keep_missing=False,
                                data_directory=data_path+'obs/')
            obs2 = data[site][0].mf.values
            times2 = data[site][0].time.values
            keep_index = []
            for t,time in enumerate(times1):
                for t2,time2 in enumerate(times2):
                    if time == time2:
                        if np.abs(obs1[t] - obs2[t2]) < 20.:
                            keep_index.append(t)
            fp_data_H_out = fp_data_H[dict(time=keep_index)]  
        else:
            print(f'{site} obs not filtered with height comparison as no second inlet height provided.')
        return fp_data_H_out
    
    n_obs = []
    n_obs_filtered = []
    for i,site in enumerate(sites):
        n_obs.append(fp_data_H[site].mf.values.shape[0])
        if len(filtering_types) == len(sites):
            filter_site = filtering_types[i]
        else:
            filter_site = filtering_types
        if 'localness' in filter_site:
            print('Applying localness filtering')
            fp_data_H = define_localness(fp_data_H,site)
            fp_data_H[site] = localness_filter(fp_data_H[site],site)
        if 'pblh_site' in filter_site:
            print('Applying PBLH filtering')
            fp_data_H[site]["pblh_threshold"] = max(int(heights[i][:-1])+100,250)
            fp_data_H[site] = pblh_filter(fp_data_H[site],site)
        if 'height' in filter_site:
            print('Appling height filtering based on obs from other heights')
            fp_data_H[site] = height_filter(fp_data_H[site],site,network[i],heights[i],secondary_heights[i],
                                      start_date,end_date,average,species)
        n_obs_filtered.append(fp_data_H[site].mf.values.shape[0])   
    perc_filtered = np.round((np.array(n_obs) - np.array(n_obs_filtered))/np.array(n_obs)*100,2)
    print(f'% of {species} filtered: {perc_filtered}') 
    
    return fp_data_H,perc_filtered